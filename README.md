# wgs-gwasqc

[![tests](https://github.com/nighpg/wgs-gwasqc/actions/workflows/tests.yml/badge.svg)](https://github.com/nighpg/wgs-gwasqc/actions/workflows/tests.yml) [![license](https://img.shields.io/github/license/nighpg/wgs-gwasqc)](LICENSE)

**English** | [日本語](README.ja.md)

Tools for **selecting GWAS-ready sites from a joint-called WGS VCF**.
They are designed for cohorts that mix sequencers, library-prep lots or sequencing providers, and they reduce
batch effects **without filtering on allele-frequency (AF) differences between batches**.

Removing sites because their AF differs between batches also removes real signals that are associated with the
phenotype or with ancestry. Instead, these tools select sites by genotype quality (missingness after genotype
masking and its difference between batches, within-batch HWE and excess heterozygosity) and by region masks.
Batch differences that remain are handled in the GWAS model (batch covariate, PCs, mixed models).

| Script | Role | Requires |
|---|---|---|
| `filter_gwas.py` | **Site filter.** Adds one FILTER label per criterion and writes the sites that pass all of them | bcftools, Python standard library |
| `batch_af_test_vcf.py` | **Site-level diagnostics.** AF difference, missingness difference and HWE between batches; λGC by frequency band | bcftools, numpy, pandas, scipy |
| `sample_qc.py` | **Sample-level diagnostics.** Per-sample variant counts, Ti/Tv, het/hom ratio, singletons, etc. compared between batches; outliers and sex mismatches | bcftools, numpy, pandas, scipy |
| `gwasqc_common.py` | Shared module used by the three scripts (keep it in the same directory) | — |

## Requirements

| Item | Tested version | Notes |
|---|---|---|
| OS | Ubuntu 24.04 (Linux x86_64) | macOS not tested |
| bcftools / htslib | 1.19 | minimum version not determined |
| Python | 3.12 | |
| numpy / pandas / scipy | 1.26 / 2.1 / 1.11 | only for `batch_af_test_vcf.py` and `sample_qc.py` |

No installation is needed: clone the repository and run the scripts directly (`bcftools` must be on `PATH`).

## Workflow

```
joint-called VCF
 ├─ 1. diagnose      batch_af_test_vcf.py   size of the batch effect (λ, missingness difference, HWE)
 ├─ 2. filter        filter_gwas.py         → <prefix>.pass.bcf
 ├─ 3. re-diagnose   batch_af_test_vcf.py   on the output of 2 with --no-gt-filter: what remains in QC-passing sites
 ├─ 4. sample QC     sample_qc.py           on the output of 2 with --no-gt-filter: rare-variant batch differences, outliers, sex mismatches
 ├─ 5. convert to plink2 → LD pruning → PCA → GWAS (batch as a covariate; REGENIE, SAIGE, ...)
 └─ 6. check significant sites against the tsv from step 1 (P_batch, dMISS, flag_*)
```

```bash
# 1. diagnose
python3 batch_af_test_vcf.py -i joint.bcf -S groups.txt --batches BatchA,BatchB -o qc/before --jobs 8
# 2. filter
python3 filter_gwas.py -i joint.bcf -S groups.txt --batches BatchA,BatchB -o qc/gwasqc \
    --hwe-groups BatchA_ctrl,BatchB_ctrl --sex sex.txt \
    --mask mhc_centromere_segdup.bed.gz --mask GRCh38_alldifficultregions.bed.gz \
    --split --ref GRCh38.fa --hard-filter --jobs 16
# 3. re-diagnose / 4. sample QC (genotypes were masked in step 2, so use --no-gt-filter)
python3 batch_af_test_vcf.py -i qc/gwasqc.pass.bcf -S groups.txt --batches BatchA,BatchB --no-gt-filter -o qc/after
python3 sample_qc.py -i qc/gwasqc.pass.bcf -S groups.txt --batches BatchA,BatchB --sex sex.txt --no-gt-filter -o qc/sampleqc
```

## Input

### VCF / BCF

| Item | Requirement |
|---|---|
| Type | Joint-genotyped multi-sample VCF (e.g. output of GATK GenotypeGVCFs). gVCF is not accepted |
| FORMAT | GT is required. GQ, DP and AD are used for genotype masking when present (missing ones are skipped with a warning) |
| `##contig` | ID and length are required (used to select chromosomes, to split the input and to detect the chrX PAR) |
| Chromosomes | Autosomes by default (`chr1`–`chr22` or `1`–`22`, detected automatically). chrX with `--sex`. chrY and mitochondria are not supported |
| Multi-allelic sites | Splitting them with `filter_gwas.py --split` is recommended. The two diagnostic scripts use biallelic sites only |
| FILTER | `PASS` or `.`. Other values (e.g. sites failing VQSR) are removed by `filter_gwas.py` |
| INFO | For `--hard-filter`: QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum (GATK standard). Missing ones are skipped |

All three scripts accept the input in the same ways.

| Option | Behaviour | Index |
|---|---|---|
| `-i file` | Splits one file into `--chunk-mb` (default 20 Mb) chunks processed with `--jobs` workers | required |
| `-i file --no-chunk` | Processes one file as a single chunk (e.g. one pre-split chunk per job in a job array) | not needed |
| `-I list.txt` | List of pre-split chunk files (one path per line, any order); each file is one unit of parallel work | not needed |

With `-I`, chunks are put in genomic order by their first record, and records that fall in the range of the
previous chunk (e.g. variants spanning a boundary) are removed from the later chunk. All chunks must have the
same samples in the same order. When you merge outputs of separate `--no-chunk` runs, remove boundary duplicates
with `bcftools concat -a -d exact`.

### groups.txt (`-S`)

Same format as `bcftools +fill-tags -S`:

```
S0001	BatchA,BatchA_ctrl
S0002	BatchA
S0003	BatchB,BatchB_ctrl
```

- Column 1: sample name (exactly as in the VCF), **at least 3 characters** (bcftools 1.19 `+fill-tags` fails on
  names of 1–2 characters; the scripts check this beforehand).
- Column 2: comma-separated group names. Letters, digits and `_` only (they become part of INFO tag names).
  The suffixes `__F` / `__M` are reserved.
- Each sample belongs to at most one batch. Extra groups for HWE (e.g. controls) may be added.
- `--batches` names the batches to compare; `--hwe-groups` names the groups in which HWE is tested.
- `filter_gwas.py` and `sample_qc.py` treat all samples as one batch (`COHORT`) when `--batches` is omitted
  (`-S` is then not needed).

### sex.txt (`--sex`, to include chrX)

```
S0001	F
S0002	M
```

Column 2 is `1`/`M`/`male` or `2`/`F`/`female`. Samples with other values or not listed are treated as
"sex unknown" and are not put in the female/male groups used on chrX.

### Region mask BED (`--mask`, optional, repeatable)

Regions to exclude, e.g. MHC, pericentromeric regions, segmental duplications, GIAB difficult regions.
3 or more columns, `.gz` allowed. Chromosome names must match the VCF.

## Genotype masking (all three scripts, on by default)

The missingness filters assume that **low-quality genotypes have been set to missing (`./.`)**.
GATK joint calling keeps low-quality genotypes as `0/0` or `0/1`, so the following masks are applied by default.

| Option | Default | Effect |
|---|---|---|
| `--gt-filter` | `FMT/GQ<20 \| FMT/DP<10 \| FMT/DP>65` | Set matching genotypes to missing (`bcftools +setGT -t q -n .`) |
| `--het-binom` | `1e-5` | Set heterozygous genotypes to missing when the binomial test on AD allele balance gives p below this value (`+setGT -t "b:AD<P"`). 0 disables |
| `--no-gt-filter` | — | No masking (for input that is already masked). A warning is printed if the input does not look masked |

- Masking twice with the same conditions gives the same result, so keeping the default is safe.
- **`DP>65` is the upper limit chosen for data with a mean depth of about 33x** (about twice the mean).
  Adjust `--gt-filter` for data with a different depth.
- If mean depth differs a lot between batches, a DP threshold alone creates missingness differences and many
  sites are removed as `BatchDiffMiss`. Consider masking on GQ only in that case.
- The defaults of `--max-miss` / `--max-dmiss` were evaluated with these default masking conditions.

## filter_gwas.py

Writes the sites that pass all criteria to `<prefix>.pass.bcf`. No criterion uses AF differences.

| Label | Criterion | Option (default) |
|---|---|---|
| (existing FILTER) | Input FILTER is not PASS / `.` | — |
| `HardFilter` | GATK recommended hard filters (SNP: QD<2, QUAL<30, SOR>3, FS>60, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8 / INDEL: QD<2, QUAL<30, FS>200, ReadPosRankSum<-20, SOR>10) | `--hard-filter` (off) |
| `Monomorphic` | Sum of AC over all alleles is 0 or AN (monomorphic after masking) | — |
| `LowMAC` | min(AC, AN-AC) < MAC | `--min-mac` (off) |
| `BatchMiss` | Missing rate ≥ threshold in any batch | `--max-miss` (0.02) |
| `BatchDiffMiss` | Missing-rate difference ≥ threshold between any two batches (2 or more batches) | `--max-dmiss` (0.01) |
| `HWE` | HWE p < threshold in any HWE group (females only on chrX outside the PAR) | `--hwe` (1e-6), `--hwe-groups` |
| `ExcHet` | Excess-heterozygosity p < threshold in any batch (females only on chrX outside the PAR) | `--exchet` (1e-6) |
| `MaleHet` | On chrX outside the PAR, the fraction of heterozygous males in any batch > threshold | `--max-male-het` (0.01) |
| `Mask` | Overlaps a region in the mask BED | `--mask` |

- The missing rate is the fraction of samples in the batch without a called genotype (`1 - NS/size`).
- For binary traits, give the **controls** of each batch to `--hwe-groups` (true associations can break HWE in cases).
- Preprocessing: `--split` (split multi-allelic sites; left-normalise with `--ref`), removal of `ALT="*"`,
  genotype masking, per-batch tags with `+fill-tags`.

Output:

| File | Contents |
|---|---|
| `<prefix>.pass.bcf` (+ `.csi`) | Sites that pass all criteria. INFO has per-batch AN/AC/AF/NS/HWE/ExcHet (per sex as well when chrX is included) |
| `<prefix>.annotated.bcf` | With `--keep-annotated`: all sites with their FILTER labels (to see why sites were removed) |
| `<prefix>.filter_summary.tsv` | Settings used (`setting` rows); site counts per label and per label combination |

## batch_af_test_vcf.py

Site-level diagnostics comparing two or more batches. **`flag_batch_af` (the AF-difference test) is for
diagnosis only**; do not use it to filter sites for GWAS.

- Batch test: chi-square on the 2×2 allele-count table for two batches (Fisher's exact test for sites with small
  expected counts); 2×k chi-square (df = k-1) for three or more batches.
- The flags `flag_mono` / `flag_miss` / `flag_dmiss` / `flag_hwe` / `flag_exchet` / `flag_malehet` and `qc_fail`
  use the same criteria as the labels of `filter_gwas.py`.
- Main options: `--batches` (required, two or more), `--hwe-groups`, `--sex`; thresholds have the same defaults
  as in `filter_gwas.py`; `--p-batch` (1e-8, diagnostic only).

Output:

| File | Contents |
|---|---|
| `<prefix>.tsv.gz` | Per site: AC/AN/NS/AF/missing rate per batch, ΔAF and Δmissingness (max-min for 3+ batches), MAF, P_batch, HWE/ExcHet, flags |
| `<prefix>.summary.txt` | By region (autosomes / chrX) × variant type × frequency band: QC failure rates, fraction of sites with an AF difference (all / QC-passing), λGC, 99th percentile of \|ΔAF\| |

Reading the results:

- **Look at λGC per frequency band.** λ over all sites is not meaningful because rare variants dominate
  (MAF < 0.1% is not reported).
- The difference between `lambda(all)` and `lambda(qc_pass)` is the batch effect that filtering removes. Inflation
  that remains in QC-passing sites may reflect differences in ancestry or non-systematic differences; handle it in
  the GWAS model rather than by removing sites, and check ancestry with a PCA projected onto a reference panel.
- λ1000 (λ rescaled to 1,000 vs 1,000 samples; two batches only) is large when a group has fewer than 1,000
  samples. For small data sets, look at λ itself.

## sample_qc.py

Compares per-sample metrics between batches. Intended to be run on the output of `filter_gwas.py`.

| Metric | Meaning |
|---|---|
| `missing_rate` | Missing rate |
| `n_nonref_snv` | Number of non-reference SNV genotypes (het + non-ref hom) |
| `het_hom` | Het / non-ref hom ratio |
| `titv` | Ti/Tv |
| `n_indel` | Number of non-reference indel genotypes |
| `n_singleton` | Number of variants carried by only this sample in the cohort |
| `mean_dp` | Mean depth (when FORMAT/DP is present) |
| `n_rare_nonref`, `rare_titv` | Non-reference SNV count and Ti/Tv at sites with MAF < `--rare-maf` (default 0.01) |
| `x_het_ratio` | With `--sex`: het / (het + non-ref hom) on chrX outside the PAR; used to detect mismatches with reported sex |

- Outliers: robust z within each batch (median and MAD; mean absolute deviation when MAD = 0); samples with
  |z| > `--outlier-z` (default 4) are reported. `rare_titv` is not used for outliers because it is unstable when a
  sample carries few rare variants.
- Batch comparison: per metric, median, quartiles and ratio to the overall median per batch, and the Kruskal–Wallis p.

Output: `<prefix>.samples.tsv` (per sample), `<prefix>.batches.tsv` (metric × batch), `<prefix>.summary.txt`.

**Batch effects that site filters cannot catch:** rare false positives that appear in only one batch leave no
trace in missingness or HWE and therefore pass `filter_gwas.py`. If `n_rare_nonref` / `rare_titv` / `n_singleton`
differ between batches in `sample_qc.py`, evaluate their impact explicitly in rare-variant analyses in addition to
using a batch covariate.

## chrX

chrX is included when `--sex` is given.

- The PAR is detected from the length of chrX in `##contig` (GRCh38 / GRCh37). Otherwise specify it, e.g.
  `--par 10001-2781479,155701383-156030895`. The PAR is treated like an autosome.
- Outside the PAR, HWE and excess heterozygosity are tested in **females only**, and `MaleHet` is based on the
  fraction of heterozygous males.
- The missing rate is the fraction of samples without a called genotype, so the criterion is the same whether male
  genotypes are coded as haploid or diploid.
- Heterozygous male genotypes stay in the VCF. Give sex information when you load the data into plink2 or similar
  tools so that males are treated as haploid.

## Use cases and limitations

**Works as is**

- Multi-batch WGS cohorts (different instruments, lots, providers), including three or more batches
- Binary traits (`--hwe-groups` with controls) and quantitative traits
- VCFs whose genotypes are already masked (`--no-gt-filter`)
- Joint-call output split by region (`-I` / `--no-chunk`)
- GWAS site QC for single-batch cohorts (`filter_gwas.py`, `sample_qc.py`)

**Needs adjusted settings**

- WGS with a mean depth other than about 30x: adjust the DP upper limit in `--gt-filter`
- Cohorts with diverse ancestry: create "ancestry × batch" groups in groups.txt and pass them to `--hwe-groups`
- Relatives or duplicated individuals: exclude them from the HWE groups
- Joint calls from tools other than GATK (e.g. DRAGEN): do not use `--hard-filter`; use the FILTER of that tool
- WES / panels: adjust the DP upper limit and restrict to target regions (not tested)

**Not supported**

- chrY, mitochondria, structural variants and CNVs
- Low-coverage WGS (below about 10x; use genotype likelihoods and imputation instead)
- Batch correction before joint calling (gVCF or read level)

**Notes**

- If the case/control ratio differs strongly between batches, remaining artifacts can become association signals.
  Consider per-batch GWAS followed by meta-analysis.
- In GWAS, use the AF difference between batches (`flag_batch_af` / `P_batch` of `batch_af_test_vcf.py`) **to check
  significant sites afterwards, not to remove sites**. For uses without a phenotype (e.g. building an AF reference
  panel), removing such sites does less harm.

## Notes on bcftools 1.19

- `+fill-tags -S` fails on sample names of 1–2 characters (the scripts check this beforehand).
- `F_MISSING` is not computed per group; the missing rate is computed from `NS`.
- Filter expressions do not support negation (`!`), so "chrX outside the PAR" is written as the complement of the PAR ranges.

## Tests

Tests with synthetic data (generated with the Python standard library and fixed seeds) are in `tests/`.

```bash
python3 tests/run_tests.py -v                 # needs bcftools; diagnostic-tool tests are skipped without numpy/pandas/scipy
KEEP_TEST_DATA=1 python3 tests/run_tests.py   # keep the generated data and outputs
```

They check: the outcome of each criterion at planted sites, that AF differences are not used, agreement between
input modes (`-i` / `-I` / `--no-chunk`), idempotent genotype masking and the warning for unmasked input, chrX PAR
detection and female-only HWE, agreement of `batch_af_test_vcf.py` flags with `filter_gwas.py` labels, the 3-batch
test (agreement with scipy `chi2_contingency`), agreement of `sample_qc.py` counts with `bcftools stats`, and
detection of planted anomalous samples. GitHub Actions (`.github/workflows/tests.yml`) runs them on Ubuntu 24.04
with bcftools 1.19.

Run time and memory for whole genomes of ~10,000 samples have not been measured.

## License

[MIT License](LICENSE)
