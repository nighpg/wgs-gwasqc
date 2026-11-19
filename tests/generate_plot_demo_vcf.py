#!/usr/bin/env python3
"""plot_batch_effect.py の見本画像を作るための、それらしいバッチエフェクトを仕込んだ合成VCFを生成する。

tests/generate_data.py (ユニットテスト用の小さなVCF) とは別に、README/docs 掲載用の画像を作る
ためのもの。対象は chr22 の短腕 (p arm, 0-13.7Mb; セントロメアの手前まで) で、2 バッチ x
--n-per-batch 人の genotype を次の2パターンで仕込む。

  1. 特定領域 (chr22:9.6-10.2Mb, p11.2 band 内) に大きな AF 差と、同時に高い欠測率差 (--hotspot-miss-rate)
     を作る。filter_gwas.py は AF 差ではサイトを除外しないが、この領域は BatchMiss/BatchDiffMiss で
     サイトごと完全に除外されるため、QC前の診断では強い dAF の山 / P_batch シグナルとして見えるが、
     QC後はこの領域自体が消え、山がきれいになくなる。これは docs/tutorial.ja.md Step 4 が説明する
     「AF差のあるサイトの大半は、片方のバッチで欠測が多いサイト」という典型例の再現。
  2. 残りの領域にうっすらノイズ: 一部のサイト (--noisy-site-frac) で BatchB の 0/0・1/1 を
     同数ずつ 0/1 に変える (--noise-frac) ことで、AC/AF を変えずにヘテロ過剰 (低レベルの
     マッピング/ゲノタイピングアーチファクトを想定) を作る。GQ は高いまま (遺伝型マスクでは
     消えない) なので、QC前の診断では ExcHet/HWE p が広く荒れるが、AF 自体は変わらないため
     flag_batch_af は立たない。filter_gwas.py の ExcHet/HWE 基準でサイトごと除外されるため、
     QC後はそのサイトが消えて落ち着く (= 補正前後でノイズが落ち着く例。hotspot の AF 差とは
     独立に、ExcHet/HWE トラックと pass 率の違いで見える)。

依存: numpy (VCF生成のみ)。--vcf-only を付けない場合はさらに bcftools と、
batch_af_test_vcf.py / filter_gwas.py / plot_batch_effect.py が要求する numpy/pandas/scipy/matplotlib。
Singularity コンテナ内での実行を推奨 (これらを全て含む)。

usage:
  generate_plot_demo_vcf.py OUTDIR [--n-per-batch 1000] [--seed 1] [--vcf-only] [--jobs 2]

  singularity exec wgs-gwasqc.sif python3 tests/generate_plot_demo_vcf.py demo_out

生成物 (OUTDIR 以下):
  demo.vcf, demo.bcf(.csi), groups.txt
  --vcf-only を付けない場合はさらに:
    qc/before.tsv.gz, qc/filtered.pass.bcf, qc/after.tsv.gz   (batch_af_test_vcf.py / filter_gwas.py)
    ideogram/chr22.png                                         (plot_batch_effect.py)
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.dirname(HERE)
BATCH = os.path.join(ROOT, 'batch_af_test_vcf.py')
FILTER = os.path.join(ROOT, 'filter_gwas.py')
PLOT = os.path.join(ROOT, 'plot_batch_effect.py')

CHROM = 'chr22'
CHROM_LEN = 50_818_468       # GRCh38 chr22 (data/cytoBandIdeo.hg38.txt.gz)
PARM_START = 50_000
PARM_END = 13_700_000        # p-arm end (cytoBandIdeo: p11.2 ends / acen (centromere) starts here)
HOTSPOT = (9_600_000, 10_200_000)   # inside p11.2 (gvar band), away from both arm/centromere boundaries
GTS = np.array(['0/0', '0/1', '1/1'])


def gt_header(samples):
    return ['##fileformat=VCFv4.2', f'##contig=<ID={CHROM},length={CHROM_LEN}>',
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples)]


def generate(a):
    """(samples, vcf_lines) を作る。バッチA/Bは前半/後半 --n-per-batch 人ずつ"""
    rng = np.random.default_rng(a.seed)
    nA = nB = a.n_per_batch
    samples = [f'S{i:05d}' for i in range(nA + nB)]
    pos = np.arange(PARM_START, PARM_END, a.step)
    n = len(pos)
    in_hotspot = (pos >= HOTSPOT[0]) & (pos <= HOTSPOT[1])

    # --- アレル頻度: 背景はbin全体に広がるよう対数一様、hotspotはBatchBだけ大きく上げる ---
    maf = 10 ** rng.uniform(np.log10(0.0005), np.log10(0.5), n)
    af_b = maf.copy()
    delta = rng.uniform(0.2, 0.35, n)
    af_b[in_hotspot] = np.clip(maf[in_hotspot] + delta[in_hotspot], 0.02, 0.9)

    gA = rng.binomial(2, np.broadcast_to(maf[:, None], (n, nA)))
    gB = rng.binomial(2, np.broadcast_to(af_b[:, None], (n, nB)))

    # --- 非hotspotの一部サイトでBatchBにヘテロ過剰 (ExcHet) を作る ---
    # 0/0 と 1/1 を同数ずつ 0/1 に変える (AC/AF は変えずヘテロ接合率だけ上げる) ので、
    # GQ を落とさず遺伝型マスクでは消えない、HWE/ExcHet フィルタでのみ除去される高品質の偽シグナルになる。
    noisy_site = np.where(~in_hotspot & (rng.random(n) < a.noisy_site_frac))[0]
    for i in noisy_site:
        idx0, idx2 = np.where(gB[i] == 0)[0], np.where(gB[i] == 2)[0]
        rng.shuffle(idx0)
        rng.shuffle(idx2)
        k = min(len(idx0), len(idx2), round(a.noise_frac * nB / 2))
        gB[i, idx0[:k]] = 1
        gB[i, idx2[:k]] = 1

    # --- 深さ・遺伝型クオリティ: 常に高品質 (遺伝型マスクでは何も落ちない) ---
    dpA = np.clip(rng.normal(30, 5, (n, nA)).round().astype(int), 12, 55)
    dpB = np.clip(rng.normal(30, 5, (n, nB)).round().astype(int), 12, 55)
    gqA = np.full((n, nA), 60)
    gqB = np.full((n, nB), 60)

    # --- hotspot全体で BatchB の一部も GQ<20 にし、欠測率差でサイトごと完全に除外されるようにする ---
    # (AF 差だけでは filter_gwas.py は除外しないが、同じ領域の欠測率差は BatchMiss/BatchDiffMiss で
    # 除外される。これにより QC後はこの領域の dAF/P_batch シグナルが完全になくなる)
    drop = in_hotspot[:, None] & (rng.random((n, nB)) < a.hotspot_miss_rate)
    gqB = np.where(drop, rng.integers(2, 15, (n, nB)), gqB)

    lines = gt_header(samples)
    for i in range(n):
        ref, alt = ('A', 'G') if i % 2 == 0 else ('C', 'T')
        fa = [f'{g}:{d}:{q}' for g, d, q in zip(GTS[gA[i]], dpA[i], gqA[i])]
        fb = [f'{g}:{d}:{q}' for g, d, q in zip(GTS[gB[i]], dpB[i], gqB[i])]
        lines.append(f'{CHROM}\t{pos[i]}\t.\t{ref}\t{alt}\t100\tPASS\t.\tGT:DP:GQ\t' + '\t'.join(fa + fb))
    groups = [f'{s}\tBatchA' for s in samples[:nA]] + [f'{s}\tBatchB' for s in samples[nA:]]
    return lines, groups


def run(cmd, log):
    log('$', ' '.join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('outdir')
    ap.add_argument('--n-per-batch', type=int, default=1000, help='1バッチあたりのサンプル数 (既定 1000)')
    ap.add_argument('--step', type=int, default=2500, help='サイトの間隔 (bp, 既定 2500)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--noisy-site-frac', type=float, default=0.15,
                     help='背景領域のうち、BatchBにヘテロ過剰ノイズを仕込むサイトの割合 (既定 0.15)')
    ap.add_argument('--noise-frac', type=float, default=0.1,
                     help='ノイズサイトで強制的にヘテロにする BatchB サンプルの割合 (既定 0.1)')
    ap.add_argument('--hotspot-miss-rate', type=float, default=0.15,
                     help='hotspot全体 (9.6-10.2Mb) で BatchB の GQ を下げる (=欠測にする) サンプルの割合 (既定 0.15)')
    ap.add_argument('--jobs', type=int, default=2, help='batch_af_test_vcf.py / filter_gwas.py の --jobs (既定 2)')
    ap.add_argument('--vcf-only', action='store_true',
                     help='VCF/groups.txt の生成のみ行い、batch_af_test_vcf.py / filter_gwas.py / plot_batch_effect.py は実行しない')
    a = ap.parse_args()
    log = lambda *x: print('[generate_plot_demo_vcf]', *x, file=sys.stderr, flush=True)

    os.makedirs(a.outdir, exist_ok=True)
    lines, groups = generate(a)
    vcf = os.path.join(a.outdir, 'demo.vcf')
    groups_path = os.path.join(a.outdir, 'groups.txt')
    with open(vcf, 'w') as fo:
        fo.write('\n'.join(lines) + '\n')
    with open(groups_path, 'w') as fo:
        fo.write('\n'.join(groups) + '\n')
    log(f'{len(lines) - 6} sites, {2 * a.n_per_batch} samples -> {vcf}')

    if a.vcf_only:
        return
    if not shutil.which('bcftools'):
        sys.exit('bcftools not found in PATH (Singularity コンテナ内で実行してください。VCFのみで良ければ --vcf-only)')

    bcf = os.path.join(a.outdir, 'demo.bcf')
    run(['bcftools', 'view', vcf, '-Ob', '-o', bcf], log)
    run(['bcftools', 'index', '-f', bcf], log)

    qc = os.path.join(a.outdir, 'qc')
    os.makedirs(qc, exist_ok=True)
    common = ['-S', groups_path, '--batches', 'BatchA,BatchB', '--jobs', str(a.jobs)]
    run([sys.executable, BATCH, '-i', bcf, '-o', os.path.join(qc, 'before'), *common], log)
    run([sys.executable, FILTER, '-i', bcf, '-o', os.path.join(qc, 'filtered'), *common], log)
    run([sys.executable, BATCH, '-i', os.path.join(qc, 'filtered.pass.bcf'), '--no-gt-filter',
         '-o', os.path.join(qc, 'after'), *common], log)

    ideogram = os.path.join(a.outdir, 'ideogram')
    run([sys.executable, PLOT, '--before', os.path.join(qc, 'before.tsv.gz'), '--after', os.path.join(qc, 'after.tsv.gz'),
         '-o', ideogram], log)
    log(f'done: {os.path.join(ideogram, CHROM + ".png")}')


if __name__ == '__main__':
    main()
