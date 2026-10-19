#!/usr/bin/env python3
"""バッチ間のバッチエフェクト診断 (サイト単位, 2 バッチ以上)

各サイトについてバッチ間のアレル頻度差・欠測率差・バッチ内 HWE / ExcHet を計算し、
頻度帯ごとの λGC・|ΔAF| と、filter_gwas.py と同じ基準で QC 不合格になるサイトの割合をまとめる。
バッチ間の AF 差 (flag_batch_af) は診断のためだけに出力する。GWAS のサイトフィルタには使わないこと
(表現型や祖先と関連する本物のシグナルまで除外してしまう)。フィルタには filter_gwas.py を使う。

usage:
  batch_af_test_vcf.py -i joint.bcf -S groups.txt -o out/batchQC --batches BatchA,BatchB[,BatchC...] \\
      [--hwe-groups BatchA_ctrl,BatchB_ctrl] [--sex sex.txt] [--no-gt-filter] [--jobs 8]
  入力の渡し方は filter_gwas.py と同じ: -i (index 付きなら分割して並列) / -i --no-chunk / -I chunks.list
    -I では一覧の順序は任意。隣接 chunk の境界で重複した行は後ろ側から除去する。

前処理 (既定で有効、filter_gwas.py と同じ): 低品質の遺伝型を欠測 (./.) にする。欠測率の診断はこれを前提にしている。
  --gt-filter 'FMT/GQ<20 | FMT/DP<10 | FMT/DP>65' / --het-binom 1e-5 / --no-gt-filter (マスク済みの入力向け)
バッチ間の検定: 2 バッチは 2x2 アレル数表のカイ二乗 (期待度数<5 は Fisher)、3 バッチ以上は 2xk カイ二乗 (自由度 k-1)。
chrX: --sex を指定すると対象に含める。PAR 以外では HWE / ExcHet を女性のみで判定し、男性ヘテロ割合を出す。
対象: 常染色体 (+ chrX)、2 アレルのサイト (多アレルは事前に分割)。

出力:
  <prefix>.tsv.gz       全サイトの統計と判定フラグ (flag_* と qc_fail は filter_gwas.py のラベルに対応)
  <prefix>.summary.txt  領域 x 変異の種類 x 頻度帯ごとのサマリ (標準エラーにも出力)
"""
import argparse
import collections
import concurrent.futures as cf
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from gwasqc_common import Groups, Mask, Regions, add_chunk_args, add_common_args, list_inputs, make_log, q, run  # noqa: E402

BINS = ['<0.1%', '0.1-1%', '1-5%', '>=5%']
FLAGS = ['flag_mono', 'flag_miss', 'flag_dmiss', 'flag_hwe', 'flag_exchet', 'flag_malehet', 'qc_fail', 'flag_batch_af']


def batch_p(AC, AN):
    """AC, AN: (サイト数, k) 配列。2 群は 2x2 カイ二乗 (期待度数<5 かつ p<1e-3 は Fisher)、k 群は 2xk カイ二乗"""
    k = AC.shape[1]
    REF = AN - AC
    tot_a, tot_r, n = AC.sum(1), REF.sum(1), AN.sum(1)
    with np.errstate(divide='ignore', invalid='ignore'):
        ea = tot_a[:, None] * AN / n[:, None]
        er = tot_r[:, None] * AN / n[:, None]
        valid = AN > 0
        x = np.where(valid, (AC - ea) ** 2 / ea + (REF - er) ** 2 / er, 0).sum(1)
        df = valid.sum(1) - 1
        p = chi2.sf(x, np.maximum(df, 1))
    p[(tot_a == 0) | (tot_r == 0)] = 1.0          # 単型
    p[df < 1] = np.nan                             # 比べられるバッチが 2 つ未満
    if k == 2:
        emin = np.minimum(np.minimum(ea[:, 0], ea[:, 1]), np.minimum(er[:, 0], er[:, 1]))
        for i in np.where((emin < 5) & np.isfinite(p) & (p < 1e-3))[0]:
            p[i] = fisher_exact([[AC[i, 0], REF[i, 0]], [AC[i, 1], REF[i, 1]]])[1]
    return p


def process_chunk(job):
    """1 chunk: bcftools で統計を出し、列と判定フラグを計算して part ファイルに書く"""
    i, cmd, cfg = job
    B, Hg, XH, XF, XM = cfg['batches'], cfg['hwe_groups'], cfg['x_hwe'], cfg['x_f'], cfg['x_m']
    cols = cfg['cols']
    out_path = os.path.join(cfg['tmp'], f'part{i:05d}.tsv')
    proc = subprocess.Popen(['bash', '-o', 'pipefail', '-c', cmd], stdout=subprocess.PIPE, text=True)
    first = last = None
    n = 0
    with open(out_path, 'w') as fo:
        try:
            reader = pd.read_csv(proc.stdout, sep='\t', names=cols, na_values='.', chunksize=cfg['chunksize'],
                                 dtype={'CHROM': str, 'REF': str, 'ALT': str})
            for d in reader:
                if not len(d):
                    continue
                num = [c for c in cols if c not in ('CHROM', 'POS', 'REF', 'ALT')]
                d[num] = d[num].apply(pd.to_numeric, errors='coerce')
                AC = d[[f'AC_{g}' for g in B]].fillna(0).values
                AN = d[[f'AN_{g}' for g in B]].fillna(0).values
                NS = d[[f'NS_{g}' for g in B]].fillna(0).values
                o = pd.DataFrame({'CHROM': d.CHROM, 'POS': d.POS, 'REF': d.REF, 'ALT': d.ALT})
                o['type'] = np.where((d.REF.str.len() == 1) & (d.ALT.str.len() == 1), 'SNV', 'INDEL')
                xnp = np.zeros(len(d), bool)
                if cfg['x']:
                    pos = d.POS.values
                    inpar = np.zeros(len(d), bool)
                    for s, e in cfg['par']:
                        inpar |= (pos >= s) & (pos <= e)
                    isx = d.CHROM.values == cfg['x']
                    xnp = isx & ~inpar
                    o['region'] = np.where(xnp, 'X', np.where(isx, 'X_PAR', 'auto'))
                else:
                    o['region'] = 'auto'
                with np.errstate(divide='ignore', invalid='ignore'):
                    AF = np.where(AN > 0, AC / AN, np.nan)
                    MISS = 1 - NS / np.array([cfg['size'][g] for g in B])   # 遺伝型が呼ばれた人の割合から (倍数性に依らない)
                for j, g in enumerate(B):
                    o[f'AC_{g}'], o[f'AN_{g}'], o[f'NS_{g}'] = AC[:, j], AN[:, j], NS[:, j]
                    o[f'AF_{g}'], o[f'MISS_{g}'] = AF[:, j], MISS[:, j]
                if len(B) == 2:
                    o['dAF'], o['dMISS'] = AF[:, 0] - AF[:, 1], MISS[:, 0] - MISS[:, 1]
                    o['dAF_range'], o['dMISS_range'] = np.abs(o.dAF), np.abs(o.dMISS)
                else:
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', RuntimeWarning)   # 全バッチ AN=0 の行
                        o['dAF_range'] = np.nanmax(AF, 1) - np.nanmin(AF, 1)
                    o['dMISS_range'] = MISS.max(1) - MISS.min(1)
                act, ant = AC.sum(1), AN.sum(1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    pool = act / ant
                o['MAF_pooled'] = np.minimum(pool, 1 - pool)
                o['P_batch'] = batch_p(AC, AN)
                for g in Hg:
                    o[f'HWE_{g}'] = d[f'HWE_{g}']
                for g in B:
                    o[f'ExcHet_{g}'] = d[f'ExcHet_{g}']
                o['HWE_all'] = d['HWE_all']
                for g in XH:
                    o[f'HWE_{g}'] = d[f'HWE_{g}']
                for g in XF:
                    o[f'ExcHet_{g}'] = d[f'ExcHet_{g}']
                for g in XM:
                    with np.errstate(divide='ignore', invalid='ignore'):
                        o[f'MaleHet_{g}'] = d[f'AC_Het_{g}'] / d[f'NS_{g}']

                # filter_gwas.py と同じ基準 (AF 差は使わない)
                anyb = lambda cs, thr: np.logical_or.reduce([(o[c] < thr).values for c in cs]) if cs else np.zeros(len(o), bool)
                o['flag_mono'] = (act == 0) | (act == ant)
                o['flag_miss'] = (MISS >= cfg['max_miss']).any(1)
                o['flag_dmiss'] = (o.dMISS_range >= cfg['max_dmiss']).values
                o['flag_hwe'] = np.where(xnp, anyb([f'HWE_{g}' for g in XH], cfg['hwe']), anyb([f'HWE_{g}' for g in Hg], cfg['hwe']))
                o['flag_exchet'] = np.where(xnp, anyb([f'ExcHet_{g}' for g in XF], cfg['exchet']),
                                            anyb([f'ExcHet_{g}' for g in B], cfg['exchet']))
                mh = np.logical_or.reduce([(o[f'MaleHet_{g}'] > cfg['max_male_het']).values for g in XM]) if XM else False
                o['flag_malehet'] = xnp & mh
                o['qc_fail'] = o[FLAGS[:6]].any(axis=1)
                o['flag_batch_af'] = o.P_batch < cfg['p_batch']          # 診断用のみ
                o.to_csv(fo, sep='\t', index=False, header=(n == 0), float_format='%.4g')
                if first is None:
                    first = (o.CHROM.iat[0], int(o.POS.iat[0]))
                last = (o.CHROM.iat[-1], int(o.POS.iat[-1]))
                n += len(o)
        except pd.errors.EmptyDataError:
            pass
    if proc.wait() != 0:
        raise RuntimeError(f'bcftools pipeline failed (exit {proc.returncode}): {cmd[:300]}')
    return i, out_path, first, last, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_chunk_args(ap)
    ap.add_argument('-o', '--prefix', required=True)
    add_common_args(ap, '比べるバッチのグループ名 (カンマ区切り, 2 つ以上)', need_two=True)
    t = ap.add_argument_group('thresholds (filter_gwas.py と同じ既定値)')
    t.add_argument('--max-miss', type=float, default=0.02, help='flag_miss: いずれかのバッチで欠測率 >= この値 (既定 0.02)')
    t.add_argument('--max-dmiss', type=float, default=0.01, help='flag_dmiss: バッチ間の欠測率の差 (最大-最小) >= この値 (既定 0.01)')
    t.add_argument('--hwe', type=float, default=1e-6, help='flag_hwe: HWE p < この値 (既定 1e-6)')
    t.add_argument('--exchet', type=float, default=1e-6, help='flag_exchet: ExcHet p < この値 (既定 1e-6)')
    t.add_argument('--max-male-het', type=float, default=0.01, help='flag_malehet: chrX 非 PAR の男性ヘテロ割合 > この値 (既定 0.01)')
    t.add_argument('--p-batch', type=float, default=1e-8, help='flag_batch_af (診断用): バッチ間 AF 差の検定 p < この値 (既定 1e-8)')
    ap.add_argument('--chunksize', type=int, default=500_000, help='1 chunk 内で一度に処理する行数')
    ap.add_argument('--tmpdir', help='一時ディレクトリ (既定: <prefix> と同じ場所)')
    a = ap.parse_args()
    log = make_log('batch_af_test')
    if not shutil.which('bcftools'):
        sys.exit('bcftools not found in PATH')

    inputs = list_inputs(a)
    hdr = run(['bcftools', 'view', '-h', inputs[0]])
    if '##FORMAT=<ID=GT,' not in hdr:
        sys.exit('FORMAT/GT がありません')
    samples = run(['bcftools', 'query', '-l', inputs[0]]).split()
    os.makedirs(os.path.dirname(os.path.abspath(a.prefix)), exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='batch_af_test.', dir=a.tmpdir or os.path.dirname(os.path.abspath(a.prefix)))
    G = Groups(a, samples, tmp, log, min_batches=2)
    R = Regions(a, hdr, log)
    M = Mask(a, hdr, inputs[0], R, log)
    B, Hg = G.batches, G.hwe_groups
    XH, XF, XM = (G.x_hwe_groups, G.x_batches_f, G.x_batches_m) if R.with_x else ([], [], [])

    # --- 既存タグを使えるか (マスクしたら必ず計算し直す) ---
    need = ([f'{t}_{g}' for g in B for t in ('AC', 'AN', 'NS', 'ExcHet')] + [f'HWE_{g}' for g in Hg + XH]
            + [f'ExcHet_{g}' for g in XF] + [f'{t}_{g}' for g in XM for t in ('AC_Het', 'NS')])
    tagged = all(f'##INFO=<ID={t},' in hdr for t in need)
    recompute = M.active or not tagged
    if not recompute:
        log('入力の既存のバッチ別タグを使います')
    hwe_all = '%HWE' if (recompute or '##INFO=<ID=HWE,' in hdr) else '.'
    fields = ([f'%{t}_{g}' for g in B for t in ('AC', 'AN', 'NS', 'ExcHet')] + [hwe_all]
              + [f'%HWE_{g}' for g in Hg + XH] + [f'%ExcHet_{g}' for g in XF]
              + [f'%{t}_{g}' for g in XM for t in ('AC_Het', 'NS')])
    cols = (['CHROM', 'POS', 'REF', 'ALT'] + [f'{t}_{g}' for g in B for t in ('AC', 'AN', 'NS', 'ExcHet')] + ['HWE_all']
            + [f'HWE_{g}' for g in Hg + XH] + [f'ExcHet_{g}' for g in XF] + [f'{t}_{g}' for g in XM for t in ('AC_Het', 'NS')])
    fmt = '%CHROM\\t%POS\\t%REF\\t%ALT\\t' + '\\t'.join(fields) + '\\n'
    tags = 'AN,AC,NS,HWE,ExcHet' + (',AC_Het' if R.with_x else '')

    def command(region, inp):
        c = [R.view_cmd(region, inp, extra='-M2 -e \'ALT="*"\'')] + M.cmds()
        if recompute:
            c.append(f'bcftools +fill-tags --no-version -Ou -- -S {q(G.path)} -t {tags}')
        c.append(f"bcftools query -f '{fmt}'")
        return ' | '.join(c)

    chunks = R.chunks(a, inputs)
    a.jobs = max(1, min(a.jobs, len(chunks)))
    log(f'{len(B)} batches, {len(R.chroms)} chromosomes, {len(chunks)} chunks, {a.jobs} jobs')
    cfg = dict(batches=B, hwe_groups=Hg, x_hwe=XH, x_f=XF, x_m=XM, cols=cols, tmp=tmp, chunksize=a.chunksize,
               size=dict(G.size), x=R.x if R.with_x else None, par=R.par, max_miss=a.max_miss,
               max_dmiss=a.max_dmiss, hwe=a.hwe, exchet=a.exchet, max_male_het=a.max_male_het, p_batch=a.p_batch)
    parts = {}
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(process_chunk, (i, command(r, f), cfg)) for i, (r, f) in enumerate(chunks)]
        for k, fu in enumerate(cf.as_completed(futs), 1):
            i, path, first, last, n = fu.result()
            parts[i] = (path, first, last, n)
            if k % max(1, len(chunks) // 20) == 0 or k == len(chunks):
                log(f'{k}/{len(chunks)} chunks done')

    # --- ゲノム順に連結し、隣接 chunk の重複行を除去しながら集計 ---
    corder = {c: j for j, c in enumerate(R.lengths)}
    order = sorted((i for i in parts if parts[i][1] is not None),
                   key=lambda i: (corder.get(parts[i][1][0], 1e9), parts[i][1][1]))
    k = len(B)
    med = chi2.ppf(.5, k - 1)
    cnt = collections.Counter()
    chis = collections.defaultdict(list)
    dafs = collections.defaultdict(list)
    pooled_only = dropped = 0
    prev = None
    first_out = True
    with gzip.open(f'{a.prefix}.tsv.gz', 'wt') as fo:
        for i in order:
            for d in pd.read_csv(parts[i][0], sep='\t', chunksize=a.chunksize, dtype={'CHROM': str, 'REF': str, 'ALT': str}):
                if prev is not None:
                    dup = (d.CHROM == prev[0]) & (d.POS <= prev[1])
                    dropped += int(dup.sum())
                    d = d[~dup]
                if not len(d):
                    continue
                d.to_csv(fo, sep='\t', index=False, header=first_out, float_format='%.4g')
                first_out = False
                rg = np.where(d.region.values == 'X', 'chrX', 'autosome')
                b = np.digitize(d.MAF_pooled.fillna(0).values, [.001, .01, .05])
                x = chi2.isf(np.clip(d.P_batch.values, 1e-300, 1), k - 1)
                okp = np.isfinite(d.P_batch.values)
                qcpass = ~d.qc_fail.values.astype(bool)
                fl = {f: d[f].values.astype(bool) for f in FLAGS}
                for r in ('autosome', 'chrX'):
                    for t in ('SNV', 'INDEL'):
                        st = (rg == r) & (d.type.values == t)
                        if not st.any():
                            continue
                        for bi, bn in enumerate(BINS):
                            s = st & (b == bi)
                            if not s.any():
                                continue
                            key = (r, t, bn)
                            cnt[key + ('n',)] += int(s.sum())
                            cnt[key + ('n_qcpass',)] += int((s & qcpass).sum())
                            for f in FLAGS:
                                cnt[key + (f,)] += int(fl[f][s].sum())
                            cnt[key + ('batch_af_qcpass',)] += int((fl['flag_batch_af'] & s & qcpass).sum())
                            if bi >= 1:
                                chis[key + ('all',)].append(x[s & okp].astype(np.float32))
                                chis[key + ('qcpass',)].append(x[s & okp & qcpass].astype(np.float32))
                                dafs[key].append(d.dAF_range.values[s & qcpass].astype(np.float32))
                auto = d.region.values != 'X'
                pooled_only += int((auto & (d.HWE_all.values < a.hwe) & ~fl['flag_hwe']).sum())
            last = parts[i][2]
            if prev is None or last[0] != prev[0] or last[1] > prev[1]:
                prev = last
    if a.input_list:
        log(f'chunks: {len(parts)} (empty {len(parts) - len(order)}), rows dropped at chunk boundaries: {dropped}')

    # --- サマリ ---
    lam = lambda v: np.median(v) / med if len(v) else np.nan
    lam1000 = (lambda l: 1 + (l - 1) * 500 * (1 / G.size[B[0]] + 1 / G.size[B[1]])) if k == 2 else (lambda l: np.nan)
    L = [f'input={a.input or a.input_list}  chroms={",".join(R.chroms)}',
         'batches: ' + '  '.join(f'{g}={G.size[g]}' for g in B) + f'   hwe groups: {",".join(Hg)}',
         f'genotype mask: {M.desc}',
         f'thresholds: miss>={a.max_miss}  dMISS(max-min)>={a.max_dmiss}  HWE<{a.hwe}  ExcHet<{a.exchet}'
         + (f'  maleHet>{a.max_male_het}' if R.with_x else '') + f'  (diagnostic) P_batch<{a.p_batch}',
         f'batch test: {"2x2 chi-square / Fisher" if k == 2 else f"2x{k} chi-square (df={k - 1})"}', '']
    if R.with_x:
        L.insert(3, f'chrX: PAR {";".join(f"{s}-{e}" for s, e in R.par)}; non-PAR HWE groups (female) {",".join(XH)}')
    rows = []
    for r in ('autosome', 'chrX'):
        for t in ('SNV', 'INDEL'):
            for bn in BINS:
                key = (r, t, bn)
                n = cnt[key + ('n',)]
                if not n:
                    continue
                cat = lambda s: np.concatenate(chis[key + (s,)]) if chis[key + (s,)] else np.array([])
                la, lq = (lam(cat('all')), lam(cat('qcpass'))) if bn != '<0.1%' else (np.nan, np.nan)
                dq = np.concatenate(dafs[key]) if dafs[key] else np.array([])
                nq = cnt[key + ('n_qcpass',)]
                row = {'region': r, 'type': t, 'MAF': bn, 'sites': n, 'qc_pass': nq / n}
                row.update({f.replace('flag_', ''): (cnt[key + (f,)] / n if (r == 'chrX' or f != 'flag_malehet') else np.nan)
                            for f in FLAGS[:6] if R.with_x or f != 'flag_malehet'})
                row.update({'batch_af(all)': cnt[key + ('flag_batch_af',)] / n,
                            'batch_af(qc_pass)': cnt[key + ('batch_af_qcpass',)] / max(nq, 1),
                            'lambda(all)': la, 'lambda(qc_pass)': lq})
                if k == 2:
                    row['lambda1000(qc_pass)'] = lam1000(lq)
                row['|dAF|q99(qc_pass)' if k == 2 else 'dAF_range_q99(qc_pass)'] = np.quantile(dq, .99) if len(dq) else np.nan
                rows.append(row)
    if rows:
        tab = pd.DataFrame(rows).set_index(['region', 'type', 'MAF'])
        with pd.option_context('display.width', 300, 'display.max_columns', 40):
            L.append('fraction of sites by region x type x pooled MAF (flags = filter_gwas.py の基準, batch_af = 診断用):')
            L.append(tab.to_string(float_format=lambda v: f'{v:.4g}', na_rep='-'))
        L.append('')
        L.append(f'lambda は頻度帯ごとに計算 (自由度 {k - 1}; MAF<0.1% は検定統計量が離散的で意味をなさないため省略)。'
                 + (' lambda1000 = 1,000 人 vs 1,000 人に換算した値 (各群が 1,000 人未満だと大きく出る)' if k == 2 else ''))
    else:
        L.append('対象サイトがありません')
    L.append(f'pooled HWE fail but per-group HWE pass (autosomes; batch-difference signature): {pooled_only}')
    txt = '\n'.join(L)
    print(txt, file=sys.stderr)
    open(f'{a.prefix}.summary.txt', 'w').write(txt + '\n')
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
