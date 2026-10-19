#!/usr/bin/env python3
"""Generate small synthetic VCFs (and expected results) for the test suite.

Only the Python standard library is used. All data are simulated with fixed seeds.

usage: generate_data.py OUTDIR
Creates, in OUTDIR:
  basic/   site filters, input modes            (2 batches, chr1/chr2 + a short chrX)
  mask/    genotype masking (GQ / DP / het AD)   (FORMAT GT:AD:DP:GQ)
  chrx/    chrX handling with sex                 (GRCh38 chrX length, PAR1)
  sample/  sample-level QC with planted anomalies
Each directory has *.vcf, groups.txt (and sex.txt), and expected.json.
"""
import json
import os
import random
import sys


def binom(r, n, p):
    return sum(r.random() < p for _ in range(n))


def gt(a):
    return ('0/0', '0/1', '1/1')[a]


def header(samples, contigs, fmt, info=(), filters=()):
    h = ['##fileformat=VCFv4.2']
    h += [f'##contig=<ID={c},length={n}>' for c, n in contigs]
    h += ['##FILTER=<ID=PASS,Description="All filters passed">']
    h += [f'##FILTER=<ID={f},Description="{f}">' for f in filters]
    h += [f'##INFO=<ID={i},Number=1,Type=Float,Description="{i}">' for i in info]
    desc = {'GT': ('1', 'String'), 'GQ': ('1', 'Integer'), 'DP': ('1', 'Integer'), 'AD': ('R', 'Integer')}
    h += [f'##FORMAT=<ID={f},Number={desc[f][0]},Type={desc[f][1]},Description="{f}">' for f in fmt]
    h.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples))
    return h


def write(path, lines):
    with open(path, 'w') as fo:
        fo.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
def basic(out):
    """300 + 100 samples. Each planted site expects one FILTER outcome (with --hwe-groups ctrl, --mask, --hard-filter, --split)."""
    r = random.Random(7)
    nA, nB = 300, 100
    sm = [f'S{i:04d}' for i in range(nA + nB)]
    ctrl = lambda j: j % 2 == 0
    write(f'{out}/groups.txt', [f'{s}\t{"BatchA" if j < nA else "BatchB"}' +
                                (f',{"BatchA" if j < nA else "BatchB"}_ctrl' if ctrl(j) else '') for j, s in enumerate(sm)])
    write(f'{out}/mask.bed', ['chr2\t10000000\t10100000'])
    info = ('QD', 'FS', 'SOR', 'MQ', 'MQRankSum', 'ReadPosRankSum')
    H = header(sm, [('chr1', 50000000), ('chr2', 50000000), ('chrX', 1000000)], ['GT', 'GQ'], info,
               ['VQSRTrancheSNP99.90to100.00'])
    rows, exp = [], {}

    def gts(kind, missA=0, p=0.3):
        g = []
        for j in range(nA + nB):
            isA = j < nA
            if kind == 'mono':
                x = '0/0'
            elif kind == 'hwe_ctrlA' and isA and ctrl(j):
                x = '0/1'
            elif kind == 'exchetB' and not isA:
                x = '0/1'
            elif kind == 'afdiff':
                x = gt(binom(r, 2, 0.1 if isA else 0.6))
            else:
                x = gt(binom(r, 2, p))
            if isA and j < missA:
                x = './.'
            g.append(x + ':60')
        return g

    def add(ch, pos, ref, alt, kind, expect, filt='PASS', qd=15, fs=1, missA=0):
        inf = f'QD={qd};FS={fs};SOR=0.7;MQ=60;MQRankSum=0;ReadPosRankSum=0'
        rows.append((ch, pos, f'{ch}\t{pos}\t.\t{ref}\t{alt}\t500\t{filt}\t{inf}\tGT:GQ\t' + '\t'.join(gts(kind, missA))))
        if isinstance(expect, dict):
            exp.update(expect)
        else:
            exp[f'{ch}:{pos}:{alt}'] = expect

    for ch in ('chr1', 'chr2'):
        for i in range(60):
            pos = 1000 + i * 700000
            if ch == 'chr2' and 10_000_000 <= pos <= 10_100_000:
                pos += 200000
            add(ch, pos, 'C', 'T', 'normal', 'PASS')
    add('chr1', 20_000_000, 'A', 'G', 'normal', 'PASS')          # chunk boundary (20 Mb)
    add('chr1', 20_000_001, 'A', 'G', 'normal', 'PASS')
    add('chr1', 101, 'G', 'A', 'mono', 'Monomorphic')
    add('chr1', 202, 'G', 'A', 'normal', 'BatchMiss;BatchDiffMiss', missA=30)
    add('chr1', 303, 'G', 'A', 'normal', 'BatchDiffMiss', missA=5)
    add('chr1', 404, 'G', 'A', 'hwe_ctrlA', 'HWE;ExcHet')
    add('chr1', 505, 'G', 'A', 'exchetB', 'HWE;ExcHet')
    add('chr1', 606, 'G', 'A', 'normal', 'HardFilter', qd=1)
    add('chr1', 707, 'G', 'A', 'normal', 'VQSRTrancheSNP99.90to100.00', filt='VQSRTrancheSNP99.90to100.00')
    add('chr1', 808, 'G', 'A', 'afdiff', 'PASS')                 # large AF difference between batches: must PASS
    add('chr1', 909, 'GA', 'G', 'normal', 'PASS', fs=150)        # indel: FS limit is 200
    add('chr1', 1010, 'GA', 'G', 'normal', 'HardFilter', fs=250)
    add('chr2', 10_050_000, 'C', 'T', 'normal', 'Mask')
    add('chr1', 1111, 'C', 'T,G', 'normal', {'chr1:1111:T': 'PASS', 'chr1:1111:G': 'Monomorphic'})
    add('chrX', 500, 'C', 'T', 'normal', {})                     # not an autosome: excluded without --sex
    rows.sort(key=lambda x: ({'chr1': 1, 'chr2': 2, 'chrX': 3}[x[0]], x[1]))
    write(f'{out}/basic.vcf', H + [x[2] for x in rows])
    json.dump({'labels': exp, 'excluded': ['chrX:500']}, open(f'{out}/expected.json', 'w'), indent=1)


# ---------------------------------------------------------------------------
def mask(out):
    """100 + 100 samples with GT:AD:DP:GQ; problems that only genotype masking reveals."""
    r = random.Random(3)
    nA = nB = 100
    sm = [f'S{i:04d}' for i in range(nA + nB)]
    write(f'{out}/groups.txt', [f'{s}\t{"BatchA" if i < nA else "BatchB"}' for i, s in enumerate(sm)])
    H = header(sm, [('chr1', 5000000)], ['GT', 'AD', 'DP', 'GQ'])

    def g(x, dp=30, gq=60, ad=None):
        ad = ad or {'0/0': (dp, 0), '0/1': (dp // 2, dp - dp // 2), '1/1': (0, dp)}[x]
        return f'{x}:{ad[0]},{ad[1]}:{dp}:{gq}'

    rows, exp = [], {}
    for k in range(40):
        rows.append((10000 + k * 1000, [g(gt(binom(r, 2, 0.3))) for _ in sm]))
    gts = [gt(binom(r, 2, 0.3)) for _ in sm]                       # 10 samples of batch A with GQ=5
    rows.append((900000, [g(x, gq=5 if i < 10 else 60) for i, x in enumerate(gts)]))
    exp['900000'] = 'BatchMiss;BatchDiffMiss'
    gts = ['0/1' if i % 2 == 0 else '0/0' for i in range(len(sm))]  # batch B hets with AD 40,2
    rows.append((910000, [g(x, dp=42, ad=(40, 2)) if (x == '0/1' and i >= nA) else g(x) for i, x in enumerate(gts)]))
    exp['910000'] = 'BatchMiss;BatchDiffMiss'
    gts = ['0/1' if i in (3, 150) else '0/0' for i in range(len(sm))]  # all ALT carriers low DP/GQ
    rows.append((920000, [g(x, dp=5, gq=12) if x == '0/1' else g(x) for x in gts]))
    exp['920000'] = 'Monomorphic'
    gts = ['0/1' if i % 3 == 0 else '0/0' for i in range(len(sm))]   # 10 samples of batch B with DP=80
    rows.append((930000, [g(x, dp=80) if 100 <= i < 110 else g(x) for i, x in enumerate(gts)]))
    exp['930000'] = 'BatchMiss;BatchDiffMiss'
    write(f'{out}/mask.vcf', H + [f'chr1\t{p}\t.\tA\tG\t100\tPASS\t.\tGT:AD:DP:GQ\t' + '\t'.join(s) for p, s in rows])
    json.dump({'labels': exp}, open(f'{out}/expected.json', 'w'), indent=1)


# ---------------------------------------------------------------------------
def chrx(out):
    """2 x 200 samples, half female. GRCh38 chrX length so that PAR is detected."""
    r = random.Random(11)
    nb = 200
    sm = [f'S{i:04d}' for i in range(2 * nb)]
    batch = lambda i: 'BatchA' if i < nb else 'BatchB'
    female = lambda i: i % 2 == 0
    write(f'{out}/groups.txt', [f'{s}\t{batch(i)}' + (f',{batch(i)}_ctrl' if i % 4 < 2 else '') for i, s in enumerate(sm)])
    write(f'{out}/sex.txt', [f'{s}\t{"F" if female(i) else "M"}' for i, s in enumerate(sm)])
    H = header(sm, [('chr1', 248956422), ('chrX', 156040895)], ['GT'])
    rows, exp = [], {}

    def site(ch, pos, kind, expect):
        g = []
        for i in range(2 * nb):
            if ch == 'chr1' or kind == 'par':
                g.append(gt(binom(r, 2, .4)))
            elif female(i):
                g.append('0/1' if kind == 'femhwe' else gt(binom(r, 2, .5)))
            elif kind == 'malehet' and r.random() < 0.10:
                g.append('0/1')
            elif kind == 'haploid':
                g.append(str(binom(r, 1, .5)))
            elif kind == 'malemiss' and i < nb and r.random() < 0.2:
                g.append('./.')
            else:
                a = binom(r, 1, .5)
                g.append(f'{a}/{a}')
        rows.append(f'{ch}\t{pos}\t.\tA\tG\t100\tPASS\t.\tGT\t' + '\t'.join(g))
        exp[f'{ch}:{pos}'] = expect

    for k in range(20):
        site('chr1', 1_000_000 + k * 1000, 'auto', 'PASS')
    site('chrX', 1_000_000, 'par', 'PASS')                 # PAR1
    site('chrX', 5_000_000, 'normal', 'PASS')              # males hom only: all-sample HWE would fail
    site('chrX', 5_001_000, 'malehet', 'MaleHet')
    site('chrX', 5_002_000, 'femhwe', 'HWE;ExcHet')
    site('chrX', 5_003_000, 'haploid', 'PASS')             # males coded haploid
    site('chrX', 5_004_000, 'malemiss', 'BatchMiss;BatchDiffMiss')
    write(f'{out}/chrx.vcf', H + rows)
    json.dump({'labels': exp}, open(f'{out}/expected.json', 'w'), indent=1)


# ---------------------------------------------------------------------------
def sample(out):
    """2 x 200 samples. Batch B carries extra rare C>A hets; three samples are anomalous."""
    r = random.Random(5)
    nb = 200
    N = 2 * nb
    sm = [f'S{i:04d}' for i in range(N)]
    female = lambda i: i % 2 == 0
    miss_s, het_s, sex_s = 7, 211, 101
    write(f'{out}/groups.txt', [f'{s}\t{"BatchA" if i < nb else "BatchB"}' for i, s in enumerate(sm)])
    write(f'{out}/sex.txt', [f'{s}\t{"F" if female(i) else "M"}' for i, s in enumerate(sm)])
    H = header(sm, [('chr1', 50000000), ('chrX', 156040895)], ['GT', 'DP'])
    lines = []
    pos = 1000
    for k in range(3000):
        pos += r.randint(5000, 20000)
        maf = r.betavariate(0.3, 3) if k % 3 else r.uniform(0.0005, 0.004)
        ref, alt = ('A', 'G') if r.random() < 2 / 3 else ('C', 'A')
        g = []
        for i in range(N):
            a = binom(r, 2, maf)
            if i == het_s and a == 0 and r.random() < 0.05:
                a = 1
            x = gt(a)
            if i == miss_s and r.random() < 0.2:
                x = './.'
            g.append(x + ':30')
        lines.append(('chr1', pos, ref, alt, g))
    for k in range(400):                                   # batch-B-only rare false positives (C>A hets)
        pos += r.randint(5000, 20000)
        g = [('0/1' if (i >= nb and r.random() < 0.004) else '0/0') + ':30' for i in range(N)]
        if not any(x.startswith('0/1') for x in g):
            g[nb + 1] = '0/1:30'
        lines.append(('chr1', pos, 'C', 'A', g))
    xpos = 5_000_000
    for k in range(300):                                   # chrX non-PAR
        xpos += r.randint(10000, 100000)
        maf = r.uniform(0.1, 0.5)
        g = []
        for i in range(N):
            a = binom(r, 2, maf) if (female(i) or i == sex_s) else 2 * binom(r, 1, maf)
            g.append(gt(a) + ':30')
        lines.append(('chrX', xpos, 'A', 'G', g))
    write(f'{out}/sample.vcf', H + [f'{c}\t{p}\t.\t{a}\t{b}\t100\tPASS\t.\tGT:DP\t' + '\t'.join(g) for c, p, a, b, g in lines])
    json.dump({'missing': sm[miss_s], 'het_excess': sm[het_s], 'sex_mismatch': sm[sex_s], 'n_chr1': 3400},
              open(f'{out}/expected.json', 'w'), indent=1)


def main():
    out = sys.argv[1]
    for name, fn in [('basic', basic), ('mask', mask), ('chrx', chrx), ('sample', sample)]:
        os.makedirs(f'{out}/{name}', exist_ok=True)
        fn(f'{out}/{name}')


if __name__ == '__main__':
    main()
