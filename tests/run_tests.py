#!/usr/bin/env python3
"""Test suite (synthetic data). Requires bcftools; diagnostic-tool tests also need numpy, pandas and scipy.

usage: python3 tests/run_tests.py [-v] [-k PATTERN]       (from anywhere)
       KEEP_TEST_DATA=1 python3 tests/run_tests.py         (keep the temporary directory)
"""
import collections
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from check_input_vcf import count_filters

FILTER = os.path.join(ROOT, 'filter_gwas.py')
BATCH = os.path.join(ROOT, 'batch_af_test_vcf.py')
SAMPLE = os.path.join(ROOT, 'sample_qc.py')
CHECK = os.path.join(ROOT, 'check_input_vcf.py')
PLOT = os.path.join(ROOT, 'plot_batch_effect.py')
try:
    import numpy  # noqa: F401
    import pandas as pd
    import scipy  # noqa: F401
    HAVE_SCI = True
except ImportError:
    HAVE_SCI = False
need_sci = unittest.skipUnless(HAVE_SCI, 'numpy / pandas / scipy not installed')
try:
    import matplotlib  # noqa: F401
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False
need_mpl = unittest.skipUnless(HAVE_SCI and HAVE_MPL, 'numpy / pandas / scipy / matplotlib not installed')
D = None   # data directory (set in setUpModule)


def read(path):
    with open(path) as f:
        return f.read()


def write(path, text):
    with open(path, 'w') as f:
        f.write(text)


def load_json(path):
    with open(path) as f:
        return json.load(f)


def sh(*cmd, check=True):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if check and r.returncode:
        raise AssertionError(f'command failed ({r.returncode}): {" ".join(map(str, cmd))}\n{r.stderr[-3000:]}')
    return r


def tool(script, *args, check=True):
    return sh(sys.executable, script, *args, check=check)


def labels(bcf, key='%CHROM:%POS:%ALT'):
    out = sh('bcftools', 'query', '-f', f'{key}\\t%FILTER\\n', bcf).stdout
    return dict(line.split('\t') for line in out.splitlines())


def records(bcf):
    return sh('bcftools', 'view', '-H', bcf).stdout


def to_bcf(vcf, index=True):
    bcf = vcf[:-4] + '.bcf'
    sh('bcftools', 'view', vcf, '-Ob', '-o', bcf)
    if index:
        sh('bcftools', 'index', '-f', bcf)
    return bcf


def basic_chunks():
    """Split basic.bcf into shuffled, partly overlapping chunks (+ one empty after chromosome selection)"""
    b = f'{D}/basic'
    lst = f'{b}/ck/list.txt'
    if os.path.exists(lst):
        return lst
    os.makedirs(f'{b}/ck', exist_ok=True)
    src = f'{b}/basic.bcf'
    sh('bcftools', 'view', src, '-r', 'chr1:1-10000000', '-Ob', '-o', f'{b}/ck/a.bcf')
    sh('bcftools', 'view', src, '-r', 'chr1:9000000-30000000', '-Oz', '-o', f'{b}/ck/b.vcf.gz')   # overlaps a
    sh('bcftools', 'view', src, '-r', 'chr1:30000001-50000000', '-Ob', '-o', f'{b}/ck/c.bcf')
    sh('bcftools', 'view', src, '-r', 'chr2', '-Ob', '-o', f'{b}/ck/d.bcf')
    sh('bcftools', 'view', src, '-r', 'chrX', '-Ob', '-o', f'{b}/ck/x.bcf')                        # empty after selection
    write(lst, '\n'.join(f'{b}/ck/{f}' for f in ('d.bcf', 'x.bcf', 'b.vcf.gz', 'c.bcf', 'a.bcf')) + '\n')
    return lst


def setUpModule():
    global D
    if not shutil.which('bcftools'):
        raise unittest.SkipTest('bcftools not found in PATH')
    D = tempfile.mkdtemp(prefix='gwasqc_test.')
    sh(sys.executable, os.path.join(HERE, 'generate_data.py'), D)
    for name in ('basic', 'mask', 'chrx', 'sample'):
        to_bcf(f'{D}/{name}/{name}.vcf')


def tearDownModule():
    if D and not os.environ.get('KEEP_TEST_DATA'):
        shutil.rmtree(D, ignore_errors=True)
    elif D:
        print(f'\ntest data kept: {D}', file=sys.stderr)


# ===========================================================================
class CheckInputVCF(unittest.TestCase):
    def test_plain_vcf(self):
        path = f'{D}/sample/sample.vcf'
        r = tool(CHECK, '-i', path)
        self.assertIn('Sample count: 400', r.stdout)
        self.assertIn('##contig=<ID=chrX,length=156040895>', r.stdout)
        self.assertIn('##FORMAT=<ID=GT', r.stdout)
        self.assertIn('##FORMAT=<ID=DP', r.stdout)
        self.assertIn('missing FORMAT fields for genotype masking: GQ, AD', r.stderr)
        self.assertIn('uncompressed VCF cannot be indexed', r.stderr)
        self.assertFalse(os.path.exists(path + '.csi'))
        self.assertFalse(os.path.exists(path + '.gz'))
        expected = collections.Counter(sh('bcftools', 'query', '-f', '%FILTER\n', path).stdout.splitlines())
        for label, count in expected.items():
            self.assertIn(f'{count:7d} {label}', r.stdout)

    def test_create_and_reuse_index(self):
        with tempfile.TemporaryDirectory(prefix='check_input.') as tmp:
            for extension, mode in (('vcf.gz', '-Oz'), ('bcf', '-Ob')):
                with self.subTest(extension=extension):
                    path = os.path.join(tmp, 'input.' + extension)
                    sh('bcftools', 'view', f'{D}/basic/basic.vcf', mode, '-o', path)
                    self.assertIn('Index: created', tool(CHECK, '-i', path).stdout)
                    index = path + '.csi'
                    before = os.stat(index).st_mtime_ns
                    self.assertIn('Index: already exists', tool(CHECK, '-i', path).stdout)
                    self.assertEqual(before, os.stat(index).st_mtime_ns)

    def test_missing_input(self):
        r = tool(CHECK, '-i', f'{D}/missing.vcf', check=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('file not found', r.stderr)

    def test_filter_limit(self):
        path = f'{D}/basic/basic.vcf'
        expected = collections.Counter(
            sh('bcftools', 'query', '-f', '%FILTER\n', path).stdout.splitlines()[:10])
        self.assertEqual(count_filters(path, limit=10), expected)

    def test_million_record_limit(self):
        with tempfile.TemporaryDirectory(prefix='check_input.') as tmp:
            path = os.path.join(tmp, 'large.vcf')
            with open(path, 'w') as f:
                f.write('##fileformat=VCFv4.2\n'
                        '##contig=<ID=chr1,length=50000000>\n'
                        '##FILTER=<ID=LowQual,Description="Low quality">\n'
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                block = 'chr1\t1\t.\tA\tG\t100\tPASS\t.\n' * 1000
                for _ in range(1000):
                    f.write(block)
                f.write('chr1\t2\t.\tA\tG\t1\tLowQual\t.\n')
            r = tool(CHECK, '-i', path)
            self.assertIn('inspected 1000000', r.stdout)
            self.assertIn('1000000 PASS', r.stdout)
            self.assertNotIn('LowQual', r.stdout)

    def test_missing_contig_length(self):
        with tempfile.TemporaryDirectory(prefix='check_input.') as tmp:
            path = os.path.join(tmp, 'input.vcf')
            write(path, read(f'{D}/basic/basic.vcf').replace(',length=50000000', ''))
            r = tool(CHECK, '-i', path)
            self.assertIn('contig length missing or invalid', r.stderr)

    def test_invalid_input(self):
        with tempfile.TemporaryDirectory(prefix='check_input.') as tmp:
            path = os.path.join(tmp, 'invalid.vcf')
            write(path, 'not a VCF\n')
            r = tool(CHECK, '-i', path, check=False)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn('bcftools failed', r.stderr)


class FilterBasic(unittest.TestCase):
    """filter_gwas.py: each criterion and input mode"""

    @classmethod
    def setUpClass(cls):
        b = f'{D}/basic'
        cls.args = ['-S', f'{b}/groups.txt', '--batches', 'BatchA,BatchB', '--hwe-groups', 'BatchA_ctrl,BatchB_ctrl',
                    '--mask', f'{b}/mask.bed', '--hard-filter', '--split', '--keep-annotated', '--jobs', '3']
        cls.exp = load_json(f'{b}/expected.json')
        tool(FILTER, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/i', *cls.args)
        cls.ref = f'{b}/out/i.annotated.bcf'

    def test_labels(self):
        got = labels(self.ref)
        for k, v in self.exp['labels'].items():
            self.assertEqual(got.get(k), v, k)
        self.assertEqual(len(got), len(self.exp['labels']))       # chrX excluded, boundary sites once

    def test_af_difference_not_used(self):
        self.assertEqual(labels(self.ref)['chr1:808:A'], 'PASS')

    def test_pass_file_contains_only_pass(self):
        got = labels(f'{D}/basic/out/i.pass.bcf')
        self.assertTrue(got and all(v == 'PASS' for v in got.values()))
        self.assertEqual(len(got), sum(v == 'PASS' for v in self.exp['labels'].values()))

    def test_summary_settings(self):
        txt = read(f'{D}/basic/out/i.filter_summary.tsv')
        self.assertIn('setting\tbatches', txt)
        self.assertIn('label\tPASS', txt)

    def test_input_list_mode_equals_single_file(self):
        b = f'{D}/basic'
        lst = basic_chunks()
        r = tool(FILTER, '-I', lst, '-o', f'{b}/out/I', *self.args)
        self.assertIn('exact CHROM:POS:REF:ALT duplicates: 2', r.stderr)
        self.assertEqual(records(f'{b}/out/I.annotated.bcf'), records(self.ref))

    def test_no_chunk_then_concat_equals_single_file(self):
        b = f'{D}/basic'
        parts = []
        for name, reg in [('p1', 'chr1:1-25000000'), ('p2', 'chr1:24000000-50000000'), ('p3', 'chr2')]:
            sh('bcftools', 'view', f'{b}/basic.bcf', '-r', reg, '-Ob', '-o', f'{b}/{name}.bcf')
            tool(FILTER, '-i', f'{b}/{name}.bcf', '--no-chunk', '-o', f'{b}/out/{name}', *self.args)
            parts.append(f'{b}/out/{name}.pass.bcf')
        sh('bcftools', 'concat', '-a', '-d', 'exact', *parts, '-Ob', '-o', f'{b}/out/nc.bcf')
        self.assertEqual(records(f'{b}/out/nc.bcf'), records(f'{D}/basic/out/i.pass.bcf'))

    def test_single_batch(self):
        b = f'{D}/basic'
        tool(FILTER, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/one', '--keep-annotated')
        got = labels(f'{b}/out/one.annotated.bcf')
        self.assertFalse(any('BatchDiffMiss' in v for v in got.values()))
        self.assertIn('BatchMiss', got['chr1:202:A'])

    def test_errors(self):
        b = f'{D}/basic'
        r = tool(FILTER, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/e', '--chroms', 'chr1,chrX', check=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('--sex', r.stderr)
        r = tool(FILTER, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/e', '--chroms', 'chr1,chrY', check=False)
        self.assertNotEqual(r.returncode, 0)
        r = tool(FILTER, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/e', '--batches', 'BatchA', check=False)
        self.assertNotEqual(r.returncode, 0)          # --batches without -S


class FilterMask(unittest.TestCase):
    """Genotype masking: default, pre-masked + --no-gt-filter, double masking, warning"""

    @classmethod
    def setUpClass(cls):
        m = f'{D}/mask'
        cls.m = m
        cls.args = ['-S', f'{m}/groups.txt', '--batches', 'BatchA,BatchB', '--keep-annotated', '--jobs', '1']
        tool(FILTER, '-i', f'{m}/mask.bcf', '-o', f'{m}/out/def', *cls.args)
        pre = f'{m}/premasked.bcf'
        r = subprocess.run(['bash', '-o', 'pipefail', '-c',
                            f"bcftools +setGT {m}/mask.bcf -Ou -- -t q -n . -i 'FMT/GQ<20 | FMT/DP<10 | FMT/DP>65' "
                            f"| bcftools +setGT -Ob -o {pre} -- -t 'b:AD<1e-5' -n ."], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        sh('bcftools', 'index', '-f', pre)
        cls.pre = pre

    def test_default_mask_labels(self):
        got = labels(f'{self.m}/out/def.annotated.bcf', '%POS')
        for k, v in load_json(f'{self.m}/expected.json')['labels'].items():
            self.assertEqual(got[k], v, k)

    def test_premasked_with_no_gt_filter_is_identical(self):
        tool(FILTER, '-i', self.pre, '-o', f'{self.m}/out/pre', '--no-gt-filter', *self.args)
        self.assertEqual(records(f'{self.m}/out/pre.annotated.bcf'), records(f'{self.m}/out/def.annotated.bcf'))

    def test_double_masking_is_identical(self):
        tool(FILTER, '-i', self.pre, '-o', f'{self.m}/out/dbl', *self.args)
        self.assertEqual(records(f'{self.m}/out/dbl.annotated.bcf'), records(f'{self.m}/out/def.annotated.bcf'))

    def test_warning_on_unmasked_input(self):
        r = tool(FILTER, '-i', f'{self.m}/mask.bcf', '-o', f'{self.m}/out/raw', '--no-gt-filter', *self.args)
        self.assertIn('WARNING: --no-gt-filter', r.stderr)
        got = labels(f'{self.m}/out/raw.annotated.bcf', '%POS')
        self.assertEqual(got['900000'], 'PASS')     # without masking the problem is invisible


class FilterChrX(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        x = f'{D}/chrx'
        cls.x = x
        cls.r = tool(FILTER, '-i', f'{x}/chrx.bcf', '-S', f'{x}/groups.txt', '--batches', 'BatchA,BatchB',
                     '--hwe-groups', 'BatchA_ctrl,BatchB_ctrl', '--sex', f'{x}/sex.txt', '--keep-annotated',
                     '-o', f'{x}/out/x', '--jobs', '4')

    def test_par_detected(self):
        self.assertIn('PAR (GRCh38) = 10001-2781479, 155701383-156030895', self.r.stderr)

    def test_labels(self):
        got = labels(f'{self.x}/out/x.annotated.bcf', '%CHROM:%POS')
        for k, v in load_json(f'{self.x}/expected.json')['labels'].items():
            self.assertEqual(got[k], v, k)

    def test_female_only_hwe_on_nonpar(self):
        out = sh('bcftools', 'query', '-r', 'chrX:5000000', '-f', '%INFO/HWE_BatchA\\t%INFO/HWE_BatchA_ctrl__F\\n',
                 f'{self.x}/out/x.annotated.bcf').stdout.split()
        self.assertLess(float(out[0]), 1e-6)          # all samples: fails (males hom only)
        self.assertGreater(float(out[1]), 1e-3)       # females (controls): passes


# ===========================================================================
@need_sci
class BatchTest(unittest.TestCase):
    """batch_af_test_vcf.py"""

    def test_flags_match_filter_gwas_chrx(self):
        x = f'{D}/chrx'
        tool(FILTER, '-i', f'{x}/chrx.bcf', '-S', f'{x}/groups.txt', '--batches', 'BatchA,BatchB',
             '--hwe-groups', 'BatchA_ctrl,BatchB_ctrl', '--sex', f'{x}/sex.txt', '--keep-annotated', '-o', f'{x}/out/f2')
        tool(BATCH, '-i', f'{x}/chrx.bcf', '-S', f'{x}/groups.txt', '--batches', 'BatchA,BatchB',
             '--hwe-groups', 'BatchA_ctrl,BatchB_ctrl', '--sex', f'{x}/sex.txt', '-o', f'{x}/out/b', '--jobs', '2')
        self.assert_flags_match(f'{x}/out/b.tsv.gz', f'{x}/out/f2.annotated.bcf')

    def test_flags_match_filter_gwas_mask(self):
        m = f'{D}/mask'
        tool(FILTER, '-i', f'{m}/mask.bcf', '-S', f'{m}/groups.txt', '--batches', 'BatchA,BatchB',
             '--keep-annotated', '-o', f'{m}/out/f3')
        tool(BATCH, '-i', f'{m}/mask.bcf', '-S', f'{m}/groups.txt', '--batches', 'BatchA,BatchB', '-o', f'{m}/out/b')
        self.assert_flags_match(f'{m}/out/b.tsv.gz', f'{m}/out/f3.annotated.bcf')

    def assert_flags_match(self, tsv, bcf):
        d = pd.read_csv(tsv, sep='\t')
        lab = labels(bcf, '%CHROM:%POS')
        m = {'flag_mono': 'Monomorphic', 'flag_miss': 'BatchMiss', 'flag_dmiss': 'BatchDiffMiss', 'flag_hwe': 'HWE',
             'flag_exchet': 'ExcHet', 'flag_malehet': 'MaleHet'}
        self.assertEqual(len(d), len(lab))
        for row in d.itertuples():
            ls = lab[f'{row.CHROM}:{row.POS}'].split(';')
            for col, name in m.items():
                if col in d:
                    self.assertEqual(bool(getattr(row, col)), name in ls, f'{row.CHROM}:{row.POS} {col}')
            self.assertEqual(bool(row.qc_fail), ls != ['PASS'], f'{row.CHROM}:{row.POS} qc_fail')

    def test_three_batches_chi_square(self):
        from scipy.stats import chi2_contingency
        x = f'{D}/chrx'
        g3 = f'{x}/groups3.txt'
        with open(g3, 'w') as fo:
            for k, line in enumerate(read(f'{x}/groups.txt').splitlines()):
                s, gs = line.split()
                b = gs.split(',')[0]
                fo.write(f'{s}\t{"BatchC" if (b == "BatchB" and k % 2) else b}\n')
        tool(BATCH, '-i', f'{x}/chrx.bcf', '-S', g3, '--batches', 'BatchA,BatchB,BatchC', '-o', f'{x}/out/k3')
        d = pd.read_csv(f'{x}/out/k3.tsv.gz', sep='\t')
        self.assertIn('dAF_range', d)
        for row in d.itertuples():
            t = [[getattr(row, f'AC_{g}') for g in ('BatchA', 'BatchB', 'BatchC')],
                 [getattr(row, f'AN_{g}') - getattr(row, f'AC_{g}') for g in ('BatchA', 'BatchB', 'BatchC')]]
            if min(sum(t[0]), sum(t[1])) == 0:
                continue
            p = chi2_contingency(t, correction=False)[1]
            self.assertAlmostEqual(row.P_batch, p, delta=max(1e-3 * p, 1e-12))

    def test_input_list_equals_single_file(self):
        b = f'{D}/basic'
        args = ['-S', f'{b}/groups.txt', '--batches', 'BatchA,BatchB', '--jobs', '3']
        tool(BATCH, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/bi', *args)
        tool(BATCH, '-I', basic_chunks(), '-o', f'{b}/out/bI', *args)
        a = sh('bash', '-c', f'zcat {b}/out/bi.tsv.gz').stdout
        c = sh('bash', '-c', f'zcat {b}/out/bI.tsv.gz').stdout
        self.assertEqual(a, c)

    def test_empty_input(self):
        b = f'{D}/basic'
        sh('bcftools', 'view', '-h', f'{b}/basic.bcf', '-Ob', '-o', f'{b}/empty.bcf')
        r = tool(BATCH, '-i', f'{b}/empty.bcf', '--no-chunk', '-S', f'{b}/groups.txt', '--batches', 'BatchA,BatchB',
                 '-o', f'{b}/out/empty')
        self.assertIn('対象サイトがありません', r.stderr)

    def test_requires_two_batches(self):
        b = f'{D}/basic'
        r = tool(BATCH, '-i', f'{b}/basic.bcf', '-S', f'{b}/groups.txt', '--batches', 'BatchA', '-o', f'{b}/out/z', check=False)
        self.assertNotEqual(r.returncode, 0)


# ===========================================================================
@need_mpl
class PlotBatchEffect(unittest.TestCase):
    """plot_batch_effect.py: before/after ideogram PNGs from batch_af_test_vcf.py diagnostics"""

    @classmethod
    def setUpClass(cls):
        b = f'{D}/basic'
        args = ['-S', f'{b}/groups.txt', '--batches', 'BatchA,BatchB']
        tool(BATCH, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/plot_before', *args)
        tool(FILTER, '-i', f'{b}/basic.bcf', '-o', f'{b}/out/plot_filtered', *args)
        tool(BATCH, '-i', f'{b}/out/plot_filtered.pass.bcf', '--no-gt-filter', '-o', f'{b}/out/plot_after', *args)
        cls.before = f'{b}/out/plot_before.tsv.gz'
        cls.after = f'{b}/out/plot_after.tsv.gz'
        cls.mask_bed = f'{b}/mask.bed'   # basic() が生成する chr2:10000000-10100000 の BED

    def test_creates_png_per_chromosome(self):
        outdir = f'{D}/basic/out/ideogram'
        tool(PLOT, '--before', self.before, '--after', self.after, '-o', outdir)
        for chrom in ('chr1', 'chr2'):
            path = os.path.join(outdir, f'{chrom}.png')
            self.assertTrue(os.path.exists(path), path)
            with open(path, 'rb') as f:
                self.assertEqual(f.read(8), b'\x89PNG\r\n\x1a\n')   # PNG magic bytes
            self.assertGreater(os.path.getsize(path), 10_000)

    def test_window_mb_option(self):
        outdir = f'{D}/basic/out/ideogram_10mb'
        tool(PLOT, '--before', self.before, '--after', self.after, '-o', outdir, '--chroms', 'chr1', '--window-mb', '10')
        self.assertTrue(os.path.exists(os.path.join(outdir, 'chr1.png')))

    def test_read_mask_bed_parses_regions_and_skips_header_lines(self):
        """read_mask_bed: filter_gwas.py の --mask と同じ BED 形式を読み、track/browser/# 行はスキップすること"""
        from plot_batch_effect import read_mask_bed
        bed = f'{D}/basic/out/_mask_with_header.bed'
        write(bed, '\n'.join(['track name=mask', '# comment', 'chr1\t1000\t2000', 'chr2\t10000000\t10100000']) + '\n')
        regions = read_mask_bed([bed])
        self.assertEqual(regions, {'chr1': [(1000, 2000)], 'chr2': [(10000000, 10100000)]})

    def test_read_mask_bed_merges_multiple_files(self):
        """--mask を複数指定した場合、同じ染色体の区間が1つにまとめられること"""
        from plot_batch_effect import read_mask_bed
        bed_a = f'{D}/basic/out/_mask_a.bed'
        bed_b = f'{D}/basic/out/_mask_b.bed'
        write(bed_a, 'chr1\t1000\t2000\nchr2\t500\t600\n')
        write(bed_b, 'chr1\t5000\t6000\n')
        regions = read_mask_bed([bed_a, bed_b])
        self.assertEqual(regions, {'chr1': [(1000, 2000), (5000, 6000)], 'chr2': [(500, 600)]})

    def test_mask_option_adds_png_without_error(self):
        """--mask を (複数回) 指定しても染色体ごとのPNGが作られること"""
        outdir = f'{D}/basic/out/ideogram_mask'
        tool(PLOT, '--before', self.before, '--after', self.after, '-o', outdir, '--chroms', 'chr1,chr2',
             '--mask', self.mask_bed, '--mask', self.mask_bed)   # 複数指定できることも確認
        for chrom in ('chr1', 'chr2'):
            self.assertTrue(os.path.exists(os.path.join(outdir, f'{chrom}.png')))

    def test_af_hwe_exchet_columns_are_per_group(self):
        """AF/HWE/ExcHet はグループ (--batches) ごとの列をそのまま保持し、色分けに使えること"""
        from plot_batch_effect import _metric_columns, group_colors
        header = pd.read_csv(self.before, sep='\t', nrows=0).columns.tolist()
        af_cols, hwe_cols, exchet_cols, _, daf_col = _metric_columns(header)
        self.assertEqual(af_cols, ['AF_BatchA', 'AF_BatchB'])
        self.assertEqual(hwe_cols, ['HWE_BatchA', 'HWE_BatchB'])
        self.assertEqual(exchet_cols, ['ExcHet_BatchA', 'ExcHet_BatchB'])
        self.assertEqual(daf_col, 'dAF')   # 2バッチなので符号付きの dAF 列がある
        colors = group_colors(['BatchA', 'BatchB'])
        self.assertNotEqual(colors['BatchA'], colors['BatchB'])

    def test_af_track_uses_raw_value_unchanged(self):
        """AF トラックは二次加工せず生の値をそのまま使うこと (1bp窓=集計なしで生データと一致することで確認)"""
        from plot_batch_effect import load_window_stats
        raw = pd.read_csv(self.before, sep='\t', usecols=['CHROM', 'POS', 'AF_BatchA', 'AF_BatchB'])
        out, _, _, _ = load_window_stats(self.before, 1, lambda *x: None)
        merged = out.merge(raw, left_on=['CHROM', 'window'], right_on=['CHROM', 'POS'], how='inner')
        self.assertGreater(len(merged), 0)
        for g in ('BatchA', 'BatchB'):
            pd.testing.assert_series_equal(merged[f'af__{g}'].reset_index(drop=True),
                                            merged[f'AF_{g}'].reset_index(drop=True), check_names=False,
                                            check_dtype=False)

    def test_daf_track_uses_raw_value_unchanged(self):
        """dAF トラックは二次加工せず生の値をそのまま使うこと (1bp窓=集計なしで生データと一致することで確認)"""
        from plot_batch_effect import load_window_stats
        raw = pd.read_csv(self.before, sep='\t', usecols=['CHROM', 'POS', 'dAF'])
        out, _, _, _ = load_window_stats(self.before, 1, lambda *x: None)
        merged = out.merge(raw, left_on=['CHROM', 'window'], right_on=['CHROM', 'POS'], how='inner')
        self.assertGreater(len(merged), 0)
        pd.testing.assert_series_equal(merged['daf'].reset_index(drop=True),
                                        merged['dAF'].reset_index(drop=True), check_names=False, check_dtype=False)

    def test_af_plot_spec_mirrors_second_group_only_for_two_groups(self):
        """AFトラックの表示向き: ちょうど2グループなら2番目を反転、それ以外は全て正の向き"""
        from plot_batch_effect import af_plot_spec
        self.assertEqual(af_plot_spec(['BatchA', 'BatchB']), [('BatchA', 1), ('BatchB', -1)])
        self.assertEqual(af_plot_spec(['BatchA']), [('BatchA', 1)])
        self.assertEqual(af_plot_spec(['BatchA', 'BatchB', 'BatchC']),
                          [('BatchA', 1), ('BatchB', 1), ('BatchC', 1)])

    def test_chromosome_without_data_is_skipped(self):
        """複数指定した染色体のうち、データが無いものだけスキップして残りは描画する"""
        outdir = f'{D}/basic/out/ideogram_chr5'
        r = tool(PLOT, '--before', self.before, '--after', self.after, '-o', outdir, '--chroms', 'chr1,chr5')
        self.assertIn('chr5: Before/After いずれにもデータがないためスキップ', r.stderr)
        self.assertTrue(os.path.exists(os.path.join(outdir, 'chr1.png')))
        self.assertFalse(os.path.exists(os.path.join(outdir, 'chr5.png')))

    def test_no_matching_chromosome_fails(self):
        """指定した染色体のどれも描画できない場合はエラーで終了する (chr5 のみ, データなし)"""
        outdir = f'{D}/basic/out/ideogram_chr5_only'
        r = tool(PLOT, '--before', self.before, '--after', self.after, '-o', outdir, '--chroms', 'chr5', check=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('描画対象の染色体がありません', r.stderr)

    def test_no_matching_chromosome_name_fails(self):
        outdir = f'{D}/basic/out/ideogram_none'
        r = tool(PLOT, '--before', self.before, '--after', self.after, '-o', outdir, '--chroms', 'chrNotReal', check=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('描画対象の染色体がありません', r.stderr)


# ===========================================================================
@need_sci
class SampleQC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        s = f'{D}/sample'
        cls.s = s
        cls.exp = load_json(f'{s}/expected.json')
        tool(SAMPLE, '-i', f'{s}/sample.bcf', '-S', f'{s}/groups.txt', '--batches', 'BatchA,BatchB',
             '--sex', f'{s}/sex.txt', '-o', f'{s}/out/q', '--jobs', '4')
        cls.d = pd.read_csv(f'{s}/out/q.samples.tsv', sep='\t').set_index('sample')
        cls.b = pd.read_csv(f'{s}/out/q.batches.tsv', sep='\t')

    def test_planted_anomalies(self):
        self.assertIn('missing_rate', str(self.d.at[self.exp['missing'], 'outlier']))
        self.assertIn('n_nonref_snv', str(self.d.at[self.exp['het_excess'], 'outlier']))
        self.assertTrue(self.d.at[self.exp['sex_mismatch'], 'sex_mismatch'])
        self.assertEqual(int(self.d.sex_mismatch.sum()), 1)

    def test_rare_batch_difference(self):
        p = self.b[self.b.metric == 'rare_titv'].kruskal_p.iloc[0]
        self.assertLess(p, 1e-6)
        med = self.b[self.b.metric == 'rare_titv'].set_index('batch')['median']
        self.assertGreater(med['BatchA'], med['BatchB'])

    def test_counts_match_bcftools_stats(self):
        out = sh('bcftools', 'stats', '-s', '-', '-t', 'chr1', f'{self.s}/sample.bcf').stdout
        for line in out.splitlines():
            if line.startswith('PSC\t'):
                f = line.split('\t')
                self.assertEqual(int(f[4]) + int(f[5]), self.d.at[f[2], 'n_nonref_snv'])
                self.assertEqual(int(f[10]), self.d.at[f[2], 'n_singleton'])
                self.assertAlmostEqual(int(f[13]) / self.exp['n_chr1'], self.d.at[f[2], 'missing_rate'], places=4)

    def test_input_list_equals_single_file(self):
        s = self.s
        os.makedirs(f'{s}/ck', exist_ok=True)
        sh('bcftools', 'view', f'{s}/sample.bcf', '-r', 'chr1:1-20000000', '-Ob', '-o', f'{s}/ck/a.bcf')
        sh('bcftools', 'view', f'{s}/sample.bcf', '-r', 'chr1:19000000-50000000', '-Ob', '-o', f'{s}/ck/b.bcf')
        sh('bcftools', 'view', f'{s}/sample.bcf', '-r', 'chrX', '-Oz', '-o', f'{s}/ck/x.vcf.gz')
        write(f'{s}/ck/list.txt', f'{s}/ck/x.vcf.gz\n{s}/ck/b.bcf\n{s}/ck/a.bcf\n')
        tool(SAMPLE, '-I', f'{s}/ck/list.txt', '-S', f'{s}/groups.txt', '--batches', 'BatchA,BatchB',
             '--sex', f'{s}/sex.txt', '-o', f'{s}/out/I', '--jobs', '3')
        self.assertEqual(read(f'{s}/out/I.samples.tsv'), read(f'{s}/out/q.samples.tsv'))

    def test_rare_artifacts_survive_site_filters(self):
        """Batch-specific rare false positives are not removed by filter_gwas.py (documented limitation)"""
        s = self.s
        tool(FILTER, '-i', f'{s}/sample.bcf', '-S', f'{s}/groups.txt', '--batches', 'BatchA,BatchB',
             '--sex', f'{s}/sex.txt', '-o', f'{s}/out/f')
        tool(SAMPLE, '-i', f'{s}/out/f.pass.bcf', '-S', f'{s}/groups.txt', '--batches', 'BatchA,BatchB',
             '--no-gt-filter', '-o', f'{s}/out/after')
        b = pd.read_csv(f'{s}/out/after.batches.tsv', sep='\t')
        self.assertLess(b[b.metric == 'rare_titv'].kruskal_p.iloc[0], 1e-6)


if __name__ == '__main__':
    unittest.main()
