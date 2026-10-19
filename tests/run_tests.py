#!/usr/bin/env python3
"""Test suite (synthetic data). Requires bcftools; diagnostic-tool tests also need numpy, pandas and scipy.

usage: python3 tests/run_tests.py [-v] [-k PATTERN]       (from anywhere)
       KEEP_TEST_DATA=1 python3 tests/run_tests.py         (keep the temporary directory)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.dirname(HERE)
FILTER = os.path.join(ROOT, 'filter_gwas.py')
BATCH = os.path.join(ROOT, 'batch_af_test_vcf.py')
SAMPLE = os.path.join(ROOT, 'sample_qc.py')
try:
    import numpy  # noqa: F401
    import pandas as pd
    import scipy  # noqa: F401
    HAVE_SCI = True
except ImportError:
    HAVE_SCI = False
need_sci = unittest.skipUnless(HAVE_SCI, 'numpy / pandas / scipy not installed')
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
