#!/usr/bin/env python3
"""サンプル単位の QC とバッチ間比較

bcftools stats -s のサンプル別カウント (PSC) を chunk ごとに並列で集計し、サンプルごとの指標を
バッチ間で比べる。サイト単位の診断 (batch_af_test_vcf.py) では見えにくい、レアバリアントのバッチ差や
個別の異常サンプルの検出に使う。filter_gwas.py の出力 (pass.bcf) に対して実行するのを想定している。

指標 (常染色体):
  missing_rate    欠測率 (遺伝型が呼ばれなかったサイトの割合)
  n_nonref_snv    非参照の SNV 遺伝型の数 (ヘテロ + 非参照ホモ)
  het_hom         ヘテロ / 非参照ホモ 比
  titv            Ti/Tv (サンプルが持つ SNV)
  n_indel         インデルの非参照遺伝型の数
  n_singleton     コホート全体で 1 人だけが持つ変異の数
  mean_dp         平均深さ (FORMAT/DP がある場合)
  n_rare_nonref   MAF < --rare-maf のサイトでの非参照 SNV 遺伝型の数
  rare_titv       同 Ti/Tv
chrX (--sex 指定時): x_het_ratio = PAR 以外のヘテロ / (ヘテロ + 非参照ホモ)。申告された性別と合わないサンプルを報告。

外れ値: 頑健な z (中央値と MAD。MAD=0 なら平均絶対偏差) をバッチ内で計算し、|z| > --outlier-z (既定 4) のサンプルを報告。
      対象: missing_rate, n_nonref_snv, het_hom, titv, n_singleton, n_rare_nonref
バッチ比較: 指標ごとにバッチ別の中央値・四分位・全体比と Kruskal-Wallis 検定の p。

usage:
  sample_qc.py -i out/gwasqc.pass.bcf -S groups.txt --batches BatchA,BatchB -o out/sampleqc [--sex sex.txt]
  入力の渡し方・遺伝型マスク (既定で有効, マスク済みなら --no-gt-filter) は filter_gwas.py と同じ。
出力:
  <prefix>.samples.tsv   サンプルごとの指標・z・外れ値フラグ
  <prefix>.batches.tsv   指標 x バッチの中央値・四分位・全体比・Kruskal-Wallis p
  <prefix>.summary.txt   サマリ (標準エラーにも出力)
"""
import argparse
import collections
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.stats import kruskal

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from gwasqc_common import Groups, Mask, Regions, add_chunk_args, add_common_args, list_inputs, make_log, q, run  # noqa: E402

PSC = ['nRefHom', 'nNonRefHom', 'nHets', 'nTransitions', 'nTransversions', 'nIndels', 'avg_dp', 'nSingletons',
       'nHapRef', 'nHapAlt', 'nMissing']
METRICS = ['missing_rate', 'n_nonref_snv', 'het_hom', 'titv', 'n_indel', 'n_singleton', 'mean_dp',
           'n_rare_nonref', 'rare_titv']
# rare_titv は 1 人あたりのレア変異が少ないと不安定なので外れ値判定には使わない (バッチ比較のみ)
OUTLIER_METRICS = ['missing_rate', 'n_nonref_snv', 'het_hom', 'titv', 'n_singleton', 'n_rare_nonref']


def stats(path, targets, include=None):
    """bcftools stats -s - → {sample: counts(np.array)}, レコード数"""
    cmd = ['bcftools', 'stats', '-s', '-', '-t', targets] + (['-i', include] if include else []) + [path]
    out = {}
    nrec = 0
    for line in run(cmd).splitlines():
        if line.startswith('SN\t') and 'number of records:' in line:
            nrec = int(line.split('\t')[-1])
        elif line.startswith('PSC\t'):
            f = line.split('\t')
            v = np.array([float(x) for x in f[3:14]])
            v[6] *= v[0] + v[1] + v[2]            # 平均深さ → 深さの合計 (呼ばれた SNV 遺伝型数で重み付け)
            out[f[2]] = v
    return out, nrec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_chunk_args(ap)
    ap.add_argument('-o', '--prefix', required=True)
    add_common_args(ap, 'バッチのグループ名 (カンマ区切り)。省略時は全サンプルを 1 バッチとして扱う')
    ap.add_argument('--rare-maf', type=float, default=0.01, help='レアバリアントの MAF 上限 (既定 0.01)')
    ap.add_argument('--outlier-z', type=float, default=4.0, help='外れ値とする頑健な z の絶対値 (既定 4)')
    ap.add_argument('--tmpdir', help='一時ディレクトリ (既定: <prefix> と同じ場所)')
    a = ap.parse_args()
    log = make_log('sample_qc')
    if not shutil.which('bcftools'):
        sys.exit('bcftools not found in PATH')

    inputs = list_inputs(a)
    hdr = run(['bcftools', 'view', '-h', inputs[0]])
    if '##FORMAT=<ID=GT,' not in hdr:
        sys.exit('FORMAT/GT がありません')
    samples = run(['bcftools', 'query', '-l', inputs[0]]).split()
    os.makedirs(os.path.dirname(os.path.abspath(a.prefix)), exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='sample_qc.', dir=a.tmpdir or os.path.dirname(os.path.abspath(a.prefix)))
    G = Groups(a, samples, tmp, log)
    R = Regions(a, hdr, log)
    M = Mask(a, hdr, inputs[0], R, log)
    chunks = R.chunks(a, inputs)
    a.jobs = max(1, min(a.jobs, len(chunks)))
    log(f'{len(G.batches)} batches, {len(R.chroms)} chromosomes, {len(chunks)} chunks, {a.jobs} jobs')

    # ---------- 1. chunk ごとにマスクして一時 BCF に ----------
    def prepare(job):
        i, (region, inp) = job
        out = os.path.join(tmp, f'c{i:05d}.bcf')
        cmd = ' | '.join([R.view_cmd(region, inp, extra='-e \'ALT="*"\'')] + M.cmds()
                         + ['bcftools view --no-version -Ob -o ' + q(out)])
        r = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(f'{region or inp}: {r.stderr.strip()[-1500:]}')
        pos = run(['bash', '-c', f"bcftools query -f '%CHROM\\t%POS\\n' {q(out)} | sed -n '1p;$p'"]).split()
        first = (pos[0], int(pos[1])) if pos else None
        last = (pos[-2], int(pos[-1])) if pos else None
        return i, out, first, last

    with cf.ThreadPoolExecutor(a.jobs) as ex:
        parts = {i: (o, f, l) for i, o, f, l in ex.map(prepare, enumerate(chunks))}

    # ---------- 2. ゲノム順に並べ、隣接 chunk の境界の重複を除去 ----------
    corder = {c: j for j, c in enumerate(R.lengths)}
    order = sorted((i for i in parts if parts[i][1]), key=lambda i: (corder.get(parts[i][1][0], 1e9), parts[i][1][1]))
    prev = None
    trimmed = 0
    files = []
    for i in order:
        out, first, last = parts[i]
        if prev and first[0] == prev[0] and first[1] <= prev[1]:
            fixed = out.replace('.bcf', '.trim.bcf')
            run(['bcftools', 'view', '--no-version', '-t', f'^{prev[0]}:1-{prev[1]}', '--targets-overlap', '0',
                 '-Ob', '-o', fixed, out])
            out = fixed
            trimmed += 1
        files.append(out)
        if prev is None or last[0] != prev[0] or last[1] > prev[1]:
            prev = last
    if trimmed:
        log(f'chunks overlapping their neighbour: {trimmed} (boundary records removed)')

    # ---------- 3. bcftools stats -s を並列で ----------
    auto = [c for c in R.chroms if c != R.x]
    t_auto = ','.join(auto)
    t_x = None
    if R.with_x:
        par = sorted(R.par)
        L = R.lengths[R.x]
        segs = [(1, par[0][0] - 1)] + [(e1 + 1, s2 - 1) for (_, e1), (s2, _) in zip(par, par[1:])] + [(par[-1][1] + 1, L)]
        t_x = ','.join(f'{R.x}:{s}-{e}' for s, e in segs if s <= e)
    rare = f'MAF<{a.rare_maf}'

    def do_stats(path):
        res = {'all': stats(path, t_auto) if auto else ({}, 0), 'rare': stats(path, t_auto, rare) if auto else ({}, 0)}
        if t_x:
            res['x'] = stats(path, t_x)
        return res

    tot = {k: collections.defaultdict(lambda: np.zeros(len(PSC))) for k in ('all', 'rare', 'x')}
    nrec = collections.Counter()
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for res in ex.map(do_stats, files):
            for k, (per, n) in res.items():
                nrec[k] += n
                for s, v in per.items():
                    tot[k][s] += v
    if nrec['all'] == 0:
        log('WARNING: 常染色体のサイトがありません')

    # ---------- 4. サンプルごとの指標 ----------
    idx = {k: i for i, k in enumerate(PSC)}
    batch_of = {}
    for line in open(G.path):
        s, gs = line.split('\t')
        for g in gs.strip().split(','):
            if g in G.batches:
                batch_of[s] = g
    rows = []
    with np.errstate(divide='ignore', invalid='ignore'):
        for s in samples:
            v, r = tot['all'][s], tot['rare'][s]
            called = v[idx['nRefHom']] + v[idx['nNonRefHom']] + v[idx['nHets']]
            row = {'sample': s, 'batch': batch_of.get(s, ''), 'sex': (G.sex or {}).get(s, ''),
                   'missing_rate': v[idx['nMissing']] / nrec['all'] if nrec['all'] else np.nan,
                   'n_nonref_snv': v[idx['nNonRefHom']] + v[idx['nHets']],
                   'het_hom': v[idx['nHets']] / v[idx['nNonRefHom']],
                   'titv': v[idx['nTransitions']] / v[idx['nTransversions']],
                   'n_indel': v[idx['nIndels']],
                   'n_singleton': v[idx['nSingletons']],
                   'mean_dp': v[idx['avg_dp']] / called if called else np.nan,
                   'n_rare_nonref': r[idx['nNonRefHom']] + r[idx['nHets']],
                   'rare_titv': r[idx['nTransitions']] / r[idx['nTransversions']]}
            if t_x:
                xv = tot['x'][s]
                row['x_het_ratio'] = xv[idx['nHets']] / (xv[idx['nHets']] + xv[idx['nNonRefHom']])
            rows.append(row)
    df = pd.DataFrame(rows)
    df = df[df.batch != '']                       # どのバッチにも属さないサンプルは除外
    if not ('FORMAT=<ID=DP,' in hdr) or df.mean_dp.fillna(0).eq(0).all():
        df['mean_dp'] = np.nan

    # ---------- 5. 外れ値 (頑健な z) ----------
    def robust_z(x):
        x = x.astype(float)
        med = np.nanmedian(x)
        dev = np.abs(x - med)
        scale = 1.4826 * np.nanmedian(dev)
        if not scale > 0:                       # 大半が同じ値 (例: 欠測率 0) → 平均絶対偏差で代用
            scale = 1.2533 * np.nanmean(dev)
        return (x - med) / scale if scale > 0 else pd.Series(np.zeros(len(x)), index=x.index)

    flags = []
    for m in OUTLIER_METRICS:
        if df[m].isna().all():
            continue
        df[f'z_{m}'] = df.groupby('batch')[m].transform(robust_z)
        df[f'zall_{m}'] = robust_z(df[m])
        flags.append(m)
    df['outlier'] = df.apply(lambda r: ','.join(m for m in flags if abs(r[f'z_{m}']) > a.outlier_z), axis=1)

    # chrX: 申告された性別との不一致
    if t_x and G.sex:
        mf = df.loc[df.sex == 'F', 'x_het_ratio'].median()
        mm = df.loc[df.sex == 'M', 'x_het_ratio'].median()
        if np.isfinite(mf) and np.isfinite(mm) and mf > mm:
            cut = (mf + mm) / 2
            df['x_sex'] = np.where(df.x_het_ratio.isna(), '', np.where(df.x_het_ratio > cut, 'F', 'M'))
            df['sex_mismatch'] = (df.sex != '') & (df.x_sex != '') & (df.sex != df.x_sex)
        else:
            log('WARNING: chrX のヘテロ比で女性と男性が分かれないため、性別の一致は判定しません')
    df.to_csv(f'{a.prefix}.samples.tsv', sep='\t', index=False, float_format='%.5g')

    # ---------- 6. バッチ比較 ----------
    brow = []
    for m in METRICS + (['x_het_ratio'] if t_x else []):
        if df[m].isna().all():
            continue
        allmed = df[m].median()
        grp = [g[m].dropna().values for _, g in df.groupby('batch')]
        p = kruskal(*grp).pvalue if len(grp) > 1 and all(len(x) > 0 for x in grp) and np.ptp(np.concatenate(grp)) > 0 else np.nan
        for b in G.batches:
            x = df.loc[df.batch == b, m].dropna()
            brow.append({'metric': m, 'batch': b, 'n': len(x), 'median': x.median(), 'q1': x.quantile(.25),
                         'q3': x.quantile(.75), 'ratio_to_all': x.median() / allmed if allmed else np.nan,
                         'kruskal_p': p})
    bt = pd.DataFrame(brow)
    bt.to_csv(f'{a.prefix}.batches.tsv', sep='\t', index=False, float_format='%.5g')

    # ---------- 7. サマリ ----------
    L = [f'input={a.input or a.input_list}  autosomes={len(auto)}' + (f'  chrX non-PAR={t_x}' if t_x else ''),
         'batches: ' + '  '.join(f'{b}={(df.batch == b).sum()}' for b in G.batches),
         f'genotype mask: {M.desc}',
         f'sites: autosomes={nrec["all"]:,}  rare (MAF<{a.rare_maf})={nrec["rare"]:,}' + (f'  chrX non-PAR={nrec["x"]:,}' if t_x else ''),
         '']
    if len(bt):
        piv = bt.pivot(index='metric', columns='batch', values='median').reindex([m for m in bt.metric.unique()])[G.batches]
        piv['kruskal_p'] = bt.groupby('metric').kruskal_p.first()
        with pd.option_context('display.width', 250, 'display.max_columns', 30):
            L.append('median by batch (kruskal_p: バッチ間で分布が違うか):')
            L.append(piv.to_string(float_format=lambda v: f'{v:.4g}'))
    L.append('')
    n_out = (df.outlier != '').sum()
    L.append(f'outliers (|robust z within batch| > {a.outlier_z}): {n_out} samples')
    for b in G.batches:
        sub = df[(df.batch == b) & (df.outlier != '')]
        if len(sub):
            L.append(f'  {b}: {len(sub)}  e.g. ' + ', '.join(f'{r.sample}({r.outlier})' for r in sub.head(5).itertuples()))
    if 'sex_mismatch' in df:
        mm = df[df.sex_mismatch]
        L.append(f'sex mismatch (reported vs chrX het ratio): {len(mm)}' +
                 ('  e.g. ' + ', '.join(f'{r.sample}({r.sex}->{r.x_sex})' for r in mm.head(5).itertuples()) if len(mm) else ''))
    if a.input_list:
        L.append('note: 境界で重複した record は除去済み')
    txt = '\n'.join(L)
    print(txt, file=sys.stderr)
    open(f'{a.prefix}.summary.txt', 'w').write(txt + '\n')
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
