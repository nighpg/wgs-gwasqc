#!/usr/bin/env python3
"""GWAS 用サイトフィルタ (バッチ間のアレル頻度差は使わない)

各基準を FILTER 列のソフトフィルタラベルとして付与し、全基準を通過したサイト (PASS) を出力する。

前処理 (既定で有効): 低品質の遺伝型を欠測 (./.) にする。欠測率のフィルタはこれを前提にしている。
  --gt-filter  'FMT/GQ<20 | FMT/DP<10 | FMT/DP>65'  に当てはまる遺伝型を欠測に
               (DP>65 は平均約 33x のデータで決めた上限。深さの違うデータでは調整する)
  --het-binom  1e-5   ヘテロで AD のアレル比が二項検定 p<1e-5 のものを欠測に (+setGT -t "b:AD<P")
  --no-gt-filter      事前にマスク済みの入力向けに両方を無効化 (同じ条件で二重にかけても結果は変わらない)

  ラベル        基準
  ------------  ----------------------------------------------------------------
  (既存)        入力の FILTER (VQSR など)。PASS 以外はそのまま除外される
  HardFilter    --hard-filter 指定時: GATK 推奨ハードフィルタ (SNP / INDEL 別)
  Monomorphic   AC=0 または AC=AN (遺伝型フィルタ後に単型)
  LowMAC        --min-mac 指定時: min(AC, AN-AC) < MAC
  BatchMiss     いずれかのバッチで欠測率 >= --max-miss
  BatchDiffMiss いずれか 2 バッチ間で |欠測率差| >= --max-dmiss (バッチが 2 つ以上のとき)
  HWE           --hwe-groups の各グループ内で HWE p < --hwe   (chrX の PAR 以外は女性のみで判定)
  ExcHet        各バッチ内で ExcHet p < --exchet               (chrX の PAR 以外は女性のみで判定)
  MaleHet       chrX の PAR 以外で、いずれかのバッチの男性のヘテロ割合 > --max-male-het
  Mask          --mask の BED と重なる (MHC, セントロメア, segdup, 難領域など)

バッチ:
  --batches 省略時は全サンプルを 1 バッチ (COHORT) として扱う (-S も不要)。BatchDiffMiss は判定しない。
chrX:
  --sex sex.txt (sample<TAB>1/M/male か 2/F/female) を指定すると chrX も対象にする。
  PAR は ##contig の chrX の長さから GRCh38/37 を判定 (--par で手動指定も可)。PAR は常染色体と同じ扱い。
  欠測率は「遺伝型が呼ばれたサンプルの割合」で計算するので、男性が 1 倍体でも 2 倍体でも同じ基準になる。
  chrY・ミトコンドリアは対象外。

usage:
  # 1 ファイル (index 付き) を --chunk-mb ごとに分割して並列処理
  filter_gwas.py -i joint.bcf -S groups.txt -o out/gwasqc \\
      --batches BatchA,BatchB [--hwe-groups BatchA_ctrl,BatchB_ctrl] [--sex sex.txt] \\
      [--mask mhc.bed --mask difficult.bed.gz] [--hard-filter] [--split --ref ref.fa] [--jobs 16]
  # 1 バッチのコホート
  filter_gwas.py -i joint.bcf -o out/gwasqc
  # 分割済みファイル群 (1 行 1 パス) をそれぞれ 1 chunk として並列処理
  ls chunks/*.bcf > chunks.list
  filter_gwas.py -I chunks.list -S groups.txt -o out/gwasqc --batches BatchA,BatchB ...
    - 一覧の順序は任意 (先頭 record の位置でゲノム順に並べ直す)
    - 隣接 chunk の境界で重複した record は後ろ側から除去する
    - 全ファイルのサンプルと順序が一致している必要がある
  # 分割済みの 1 ファイルをそのまま 1 chunk として処理 (ジョブ配列向け, index 不要)
  filter_gwas.py -i chunk_0001.bcf --no-chunk -S groups.txt -o out/chunk_0001 --batches BatchA,BatchB ...

groups.txt: sample<TAB>group[,group...] (bcftools +fill-tags -S 形式)。サンプル名は 3 文字以上。
  例) S0001  BatchA,BatchA_ctrl
出力:
  <prefix>.pass.bcf(.csi)       全基準を通過したサイト
  <prefix>.annotated.bcf(.csi)  --keep-annotated 指定時: 全サイト + FILTER ラベル
  <prefix>.filter_summary.tsv   設定、ラベル別・組み合わせ別のサイト数
"""
import argparse
import collections
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from gwasqc_common import Groups, Mask, Regions, add_chunk_args, add_common_args, list_inputs, make_log, q, run  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
add_chunk_args(ap)
ap.add_argument('-o', '--prefix', required=True)
add_common_args(ap, 'バッチのグループ名 (カンマ区切り)。省略時は全サンプルを 1 バッチとして扱う')
f = ap.add_argument_group('site filters')
f.add_argument('--max-miss', type=float, default=0.02, help='バッチ内欠測率の上限 (既定 0.02)')
f.add_argument('--max-dmiss', type=float, default=0.01, help='バッチ間欠測率差の上限 (既定 0.01)')
f.add_argument('--hwe', type=float, default=1e-6, help='HWE p の下限 (既定 1e-6)')
f.add_argument('--exchet', type=float, default=1e-6, help='ExcHet p の下限 (既定 1e-6)')
f.add_argument('--max-male-het', type=float, default=0.01, help='chrX の PAR 以外での男性ヘテロ割合の上限 (既定 0.01)')
f.add_argument('--min-mac', type=int, default=0, help='最小マイナーアレルカウント (既定 0 = 無効)')
f.add_argument('--mask', action='append', default=[], help='除外領域 BED (複数指定可, .gz 可)')
f.add_argument('--hard-filter', action='store_true', help='GATK 推奨ハードフィルタを追加')
p = ap.add_argument_group('preprocessing / output')
p.add_argument('--split', action='store_true', help='多アレルサイトを分割 (bcftools norm -m -any)')
p.add_argument('--ref', help='--split 時に左寄せ正規化する参照 FASTA')
p.add_argument('--keep-annotated', action='store_true', help='FILTER ラベル付きの全サイトも出力')
p.add_argument('--tmpdir', help='一時ディレクトリ (既定: <prefix> と同じ場所)')
a = ap.parse_args()

log = make_log('filter_gwas')
if not shutil.which('bcftools'):
    sys.exit('bcftools not found in PATH')

# ---------- 入力・グループ・領域・マスク ----------
inputs = list_inputs(a)
hdr = run(['bcftools', 'view', '-h', inputs[0]])
if '##FORMAT=<ID=GT,' not in hdr:
    sys.exit('FORMAT/GT がありません')
samples = run(['bcftools', 'query', '-l', inputs[0]]).split()
os.makedirs(os.path.dirname(os.path.abspath(a.prefix)), exist_ok=True)
tmp = tempfile.mkdtemp(prefix='filter_gwas.', dir=a.tmpdir or os.path.dirname(os.path.abspath(a.prefix)))
G = Groups(a, samples, tmp, log)
R = Regions(a, hdr, log)
M = Mask(a, hdr, inputs[0], R, log)
batches, hwe_groups, size = G.batches, G.hwe_groups, G.size
chroms = R.chroms
corder = {c: i for i, c in enumerate(R.lengths)}
chunks = R.chunks(a, inputs)
a.jobs = max(1, min(a.jobs, len(chunks)))
mode = 'input files' if a.input_list else ('no split' if a.no_chunk else f'{a.chunk_mb}Mb')
log(f'{len(chroms)} chromosomes, {len(chunks)} chunks ({mode}), {a.jobs} jobs')

# ---------- マスク BED を 1 つにまとめる ----------
mask = None
if a.mask:
    mask = os.path.join(tmp, 'mask.bed')
    with open(mask, 'w') as fo:
        for b in a.mask:
            src = subprocess.Popen(['gzip', '-dcf', b], stdout=subprocess.PIPE, text=True)
            for line in src.stdout:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 3 and not line.startswith(('#', 'track', 'browser')):
                    fo.write(f'{f[0]}\t{f[1]}\t{f[2]}\n')
            src.wait()
    subprocess.run(f'sort -k1,1 -k2,2n {q(mask)} -o {q(mask)}', shell=True, check=True)

# ---------- フィルタ式 ----------
XNP = R.expr_x_nonpar() if R.with_x else None          # chrX の PAR 以外
NOT_XNP = R.expr_not_x_nonpar() if R.with_x else None  # それ以外 (常染色体 + PAR)


def by_region(auto_terms, x_terms):
    """常染色体 (+PAR) と chrX 非 PAR で別の判定式を使う"""
    auto = ' || '.join(auto_terms)
    if not R.with_x:
        return auto
    parts = [f'({NOT_XNP} && ({auto}))'] if auto_terms else []
    if x_terms:
        parts.append(f'({XNP} && ({" || ".join(x_terms)}))')
    return ' || '.join(parts)


filters = []  # (label, expression)
if a.hard_filter:
    has = lambda t: f'ID={t},' in hdr
    snp = [e for t, e in [('QD', 'INFO/QD<2'), (None, 'QUAL<30'), ('SOR', 'INFO/SOR>3'), ('FS', 'INFO/FS>60'),
                          ('MQ', 'INFO/MQ<40'), ('MQRankSum', 'INFO/MQRankSum<-12.5'),
                          ('ReadPosRankSum', 'INFO/ReadPosRankSum<-8')] if t is None or has(t)]
    ind = [e for t, e in [('QD', 'INFO/QD<2'), (None, 'QUAL<30'), ('FS', 'INFO/FS>200'),
                          ('ReadPosRankSum', 'INFO/ReadPosRankSum<-20'), ('SOR', 'INFO/SOR>10')] if t is None or has(t)]
    missing = [t for t in ('QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum') if not has(t)]
    if missing:
        log('WARNING: ハードフィルタ用の INFO が無いため省略:', ','.join(missing))
    filters.append(('HardFilter', f'(TYPE="snp" && ({" || ".join(snp)})) || (TYPE!="snp" && ({" || ".join(ind)}))'))
# 多アレル (未分割) では AC がアレルごとのベクトルで、比較は「いずれかの要素」で真になるため合計で判定
filters.append(('Monomorphic', 'SUM(INFO/AC)=0 || SUM(INFO/AC)=INFO/AN'))
if a.min_mac > 0:
    filters.append(('LowMAC', f'INFO/AC<{a.min_mac} || (INFO/AN-INFO/AC)<{a.min_mac}'))
miss = {g: f'(1-INFO/NS_{g}/{size[g]})' for g in batches}
filters.append(('BatchMiss', ' || '.join(f'{miss[g]}>={a.max_miss}' for g in batches)))
pairs = [(x, y) for i, x in enumerate(batches) for y in batches[i + 1:]]
if pairs:
    filters.append(('BatchDiffMiss', ' || '.join(f'abs({miss[x]}-{miss[y]})>={a.max_dmiss}' for x, y in pairs)))
filters.append(('HWE', by_region([f'INFO/HWE_{g}<{a.hwe}' for g in hwe_groups],
                                 [f'INFO/HWE_{g}<{a.hwe}' for g in G.x_hwe_groups])))
filters.append(('ExcHet', by_region([f'INFO/ExcHet_{g}<{a.exchet}' for g in batches],
                                    [f'INFO/ExcHet_{g}<{a.exchet}' for g in G.x_batches_f])))
if R.with_x and G.x_batches_m:
    # 男性のヘテロ割合 = AC_Het (ヘテロ遺伝型のアレル数 = ヘテロの人数) / NS
    filters.append(('MaleHet', f'{XNP} && (' + ' || '.join(
        f'INFO/AC_Het_{g}>{a.max_male_het}*INFO/NS_{g}' for g in G.x_batches_m) + ')'))
fill_tags = 'AN,AC,AF,NS,HWE,ExcHet' + (',AC_Het' if R.with_x else '')


def pipeline(region, inp, out):
    nv = '--no-version'
    c = [R.view_cmd(region, inp)]
    if a.split:
        c.append(f'bcftools norm {nv} -m -any {"-f " + q(a.ref) if a.ref else ""} -Ou')
    c.append(f'bcftools view {nv} -e \'ALT="*"\' -Ou')
    c += M.cmds()
    c.append(f'bcftools +fill-tags {nv} -Ou -- -S {q(G.path)} -t {fill_tags}')
    for lab, ex in filters:
        c.append(f'bcftools filter {nv} -m+ -s {lab} -e {q(ex)} -Ou')
    if mask:
        c.append(f'bcftools filter {nv} -m+ -s Mask -M {q(mask)} -Ou')
    c[-1] = c[-1].replace(' -Ou', f' -Ob -o {q(out)}')
    return ' | '.join(c)


def stats_of(path, region=None):
    """FILTER x TYPE の件数と、先頭・末尾 record の (CHROM, POS)"""
    cmd = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%FILTER\t%TYPE\n'] + \
          (['-r', region, '--regions-overlap', '0'] if region else []) + [path]
    cnt, first, last = collections.Counter(), None, None
    for line in run(cmd).splitlines():
        ch, pos, filt, typ = line.split('\t')
        cnt[f'{filt}\t{typ}'] += 1
        last = (ch, int(pos))
        if first is None:
            first = last
    return cnt, first, last


def work(i_chunk):
    i, (region, inp) = i_chunk
    out = os.path.join(tmp, f'chunk{i:05d}.bcf')
    r = subprocess.run(['bash', '-o', 'pipefail', '-c', pipeline(region, inp, out)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f'{region or inp}: {r.stderr.strip()[-2000:]}')
    run(['bcftools', 'index', '-f', out])
    cnt, first, last = stats_of(out)
    return i, out, cnt, first, last


log('pipeline example:', pipeline(*chunks[0], 'chunk.bcf'))
results = {}
with cf.ThreadPoolExecutor(a.jobs) as ex:
    futs = [ex.submit(work, x) for x in enumerate(chunks)]
    for n, fu in enumerate(cf.as_completed(futs), 1):
        i, out, cnt, first, last = fu.result()
        results[i] = dict(out=out, cnt=cnt, first=first, last=last)
        if n % max(1, len(chunks) // 20) == 0 or n == len(chunks):
            log(f'{n}/{len(chunks)} chunks done')

# ---------- 分割済み入力: ゲノム順に並べ、隣接 chunk 間の重複 record を除去 ----------
order = [i for i in sorted(results) if results[i]['first'] is not None]
if a.input_list:
    empty = len(results) - len(order)
    order.sort(key=lambda i: (corder[results[i]['first'][0]], results[i]['first'][1]))
    n_ov = n_drop = n_exact = 0
    for prev, cur in zip(order, order[1:]):
        P, C = results[prev], results[cur]
        if P['last'][0] != C['first'][0] or C['first'][1] > P['last'][1]:
            continue
        # cur の先頭が prev の末尾以前 → cur 側の POS<=prev 末尾 の record を落とす
        ch, end = P['last']
        reg = f'{ch}:1-{end}'
        n_ov += 1
        rm_cnt, _, _ = stats_of(C['out'], reg)
        keys = lambda f, r: set(run(['bcftools', 'query', '-r', r, '--regions-overlap', '0',
                                     '-f', '%CHROM:%POS:%REF:%ALT\n', f]).split())
        dup = keys(C['out'], reg) & keys(P['out'], f'{ch}:{C["first"][1]}-{end}')
        n_rm = sum(rm_cnt.values())
        n_drop += n_rm
        n_exact += len(dup)
        fixed = C['out'].replace('.bcf', '.dedup.bcf')
        run(['bcftools', 'view', '--no-version', '-t', f'^{reg}', '--targets-overlap', '0', '-Ob', '-o', fixed, C['out']])
        run(['bcftools', 'index', '-f', fixed])
        C['out'] = fixed
        C['cnt'] = C['cnt'] - rm_cnt
    log(f'input files: {len(results)} (empty after chromosome selection: {empty}); '
        f'overlapping neighbours: {n_ov}, dropped records: {n_drop} (exact CHROM:POS:REF:ALT duplicates: {n_exact})')
    if n_drop > n_exact:
        log(f'WARNING: 重複ではないのに隣の chunk の範囲に入っていた record を {n_drop - n_exact} 件除去しました。'
            ' chunk が領域ごとに分かれているか確認してください')

files = [results[i]['out'] for i in order]
if not files:   # 対象染色体の record が 1 件も無い → ヘッダーのみの出力を作る
    log('WARNING: 対象染色体の record がありません。ヘッダーのみのファイルを出力します')
    files = [results[min(results)]['out']]
lst = os.path.join(tmp, 'chunks.txt')
open(lst, 'w').write('\n'.join(files) + '\n')

def concat(out, extra=''):
    cmd = (f'bcftools concat --no-version -f {q(lst)} -Ou '
           f'| bcftools view --no-version {extra} --threads 4 -Ob -o {q(out)} && bcftools index -f {q(out)}')
    subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], check=True)


if a.keep_annotated:
    concat(f'{a.prefix}.annotated.bcf')
    log('written', f'{a.prefix}.annotated.bcf')
concat(f'{a.prefix}.pass.bcf', '-f PASS')
log('written', f'{a.prefix}.pass.bcf')

# ---------- サマリ ----------
tot = collections.Counter()
for r in results.values():
    tot.update(r['cnt'])
by_label = collections.Counter()
by_type = collections.Counter()
pass_type = collections.Counter()
for k, v in tot.items():
    filt, typ = k.split('\t')
    by_type[typ] += v
    if filt == 'PASS':
        pass_type[typ] += v
    for lab in filt.split(';'):
        by_label[(lab, typ)] += v
types = sorted(by_type)
with open(f'{a.prefix}.filter_summary.tsv', 'w') as fo:
    fo.write('section\tlabel\t' + '\t'.join(types) + '\ttotal\tfraction\n')
    n_all = sum(by_type.values())
    pad = '\t' * len(types)
    for k, v in [('gt_filter', M.gt_filter or 'none'), ('het_binom', M.het_binom or 'none'),
                 ('max_miss', a.max_miss), ('max_dmiss', a.max_dmiss if pairs else 'n/a (1 batch)'),
                 ('hwe', a.hwe), ('exchet', a.exchet),
                 ('batches', ','.join(f'{g}={size[g]}' for g in batches)), ('hwe_groups', ','.join(hwe_groups)),
                 ('chroms', ','.join(chroms)),
                 ('chrX', (f'PAR={";".join(f"{s}-{e}" for s, e in R.par)} max_male_het={a.max_male_het} '
                           f'hwe_groups(nonPAR)={",".join(G.x_hwe_groups)}') if R.with_x else 'not included')]:
        fo.write(f'setting\t{k}\t{pad}{v}\t\n')
    fo.write('input\tall\t' + '\t'.join(str(by_type[t]) for t in types) + f'\t{n_all}\t1\n')
    labs = sorted({l for l, _ in by_label}, key=lambda l: (-sum(by_label[(l, t)] for t in types), l))
    for l in labs:
        s = sum(by_label[(l, t)] for t in types)
        fo.write(f'label\t{l}\t' + '\t'.join(str(by_label[(l, t)]) for t in types) + f'\t{s}\t{s / n_all:.5f}\n')
    combos = collections.Counter()
    for k, v in tot.items():
        combos[k.split('\t')[0]] += v
    for c, v in sorted(combos.items(), key=lambda x: (-x[1], x[0]))[:30]:
        fo.write(f'combination\t{c}\t' + '\t'.join([''] * len(types)) + f'\t{v}\t{v / n_all:.5f}\n')
n_all = sum(by_type.values())
n_pass = sum(pass_type.values())
log(f'sites: input={n_all:,}  PASS={n_pass:,} ({n_pass / max(n_all, 1):.3f})  ' +
    '  '.join(f'{t}:{pass_type[t]:,}/{by_type[t]:,}' for t in types))
log('summary:', f'{a.prefix}.filter_summary.tsv')
shutil.rmtree(tmp, ignore_errors=True)
