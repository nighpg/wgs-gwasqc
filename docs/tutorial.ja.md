# チュートリアル: NovaSeq 6000 と NovaSeq X Plus の混在データからバッチエフェクトを除いた VCF を作る

[README (日本語)](../README.ja.md) | [README (English)](../README.md)

NovaSeq 6000 と NovaSeq X Plus で読んだサンプルが混ざった WGS コホートで、joint call 後の VCF から
**シーケンサー間のバッチエフェクトを除いた GWAS 用の VCF** を作るまでの手順です。

## 想定するデータ

| 項目 | 想定 |
|---|---|
| データ | ショートリード WGS (平均 30x 前後)、GRCh38 |
| 処理 | bwa-mem → GATK HaplotypeCaller (gVCF) → GenotypeGVCFs で joint call → **VQSR** |
| 入力 VCF | **フィルタリング前**: FILTER 列に VQSR の判定 (`PASS` / `VQSRTranche...`) が入っていて、遺伝型のマスクはしていない |
| バッチ | NovaSeq 6000 (`NS6000`) と NovaSeq X Plus (`NSX`) の 2 群 |
| 表現型 | 2 値形質 (患者・対照) の例で説明。量的形質の場合の違いは各所に記載 |

典型的な症状: QC 前の VCF で PCA をすると、シーケンサーごとに明確なクラスターが分かれる。

## 全体の流れ

| Step | 内容 | 使うもの |
|---|---|---|
| 0 | 準備 (ツール、作業ディレクトリ) | git, bcftools, Python |
| 1 | 入力 VCF の確認 | bcftools |
| 2 | groups.txt / sex.txt / マスク用 BED の作成 | awk, curl |
| 3 | 深さの確認と遺伝型マスク条件の決定 | bcftools |
| 4 | QC 前の診断 | `batch_af_test_vcf.py` |
| 5 | サイトフィルタ | `filter_gwas.py` |
| 6 | QC 後の診断 | `batch_af_test_vcf.py` |
| 7 | サンプル QC (必要なら除外して Step 5 から再実行) | `sample_qc.py` |
| 8 | PCA でバッチの分離が消えたか確認 | plink2 |
| 9 | GWAS と有意サイトの事後確認 | REGENIE / SAIGE など |

以下、ファイル名は次のように置きます。

```bash
VCF=joint.vqsr.vcf.gz        # joint call + VQSR 後の VCF (index 付き: joint.vqsr.vcf.gz.tbi か .csi)
REF=GRCh38.fa                # 参照配列 (.fai 付き)
TOOLS=$HOME/wgs-gwasqc       # このリポジトリ
mkdir -p qc
```

---

## Step 0. 準備

```bash
git clone https://github.com/nighpg/wgs-gwasqc.git $TOOLS
bcftools --version | head -1          # 1.19 で検証済み
python3 -c "import numpy, pandas, scipy"   # 診断ツール (Step 4, 6, 7) に必要
python3 $TOOLS/tests/run_tests.py     # 動作確認 (合成データ, 十数秒)
```

`plink2` (Step 8) と GWAS のツール (Step 9) は別途用意してください。

## Step 1. 入力 VCF の確認

```bash
# サンプル数
bcftools query -l $VCF | wc -l

# ##contig に length があるか (chrX の PAR 判定にも使う)
bcftools view -h $VCF | grep -m3 '^##contig'

# FORMAT に GT, GQ, DP, AD があるか (遺伝型マスクに使う)
bcftools view -h $VCF | grep '^##FORMAT' | cut -d, -f1

# FILTER の内訳 (先頭 100 万行)。VQSR の判定が入っていることを確認
bcftools query -f '%FILTER\n' $VCF | head -1000000 | sort | uniq -c

# index が無ければ作る (Step 4-6 の並列処理に必要)
[ -e $VCF.tbi ] || [ -e $VCF.csi ] || bcftools index $VCF
```

- VQSR で `PASS` 以外になったサイトは `filter_gwas.py` が除外します。VQSR を使っているので、GATK のハードフィルタ (`--hard-filter`) は不要です。
- 多アレルサイトは Step 5 の `--split` で分割します。
- 領域ごとに分割された VCF (チャンク) の場合は、各ステップの `-i $VCF` を `-I chunks.list` (1 行 1 パス) に置き換えてください。

## Step 2. groups.txt / sex.txt / マスク用 BED

### groups.txt

各サンプルのバッチ (シーケンサー) と、HWE を判定するグループ (2 値形質なら対照群) を書きます。
例として、サンプル情報が `samples.tsv` (列: `sample` `sequencer` `status` `sex`) にあるとします。

```
sample	sequencer	status	sex
S0001	NovaSeq 6000	control	F
S0002	NovaSeq X Plus	case	M
```

```bash
awk -F'\t' 'NR>1 {
  b = ($2 ~ /6000/) ? "NS6000" : "NSX"
  g = b; if ($3 == "control") g = g "," b "_ctrl"
  print $1 "\t" g
}' samples.tsv > groups.txt

awk -F'\t' 'NR>1 {print $1 "\t" $4}' samples.tsv > sex.txt

# 確認: バッチ・対照群ごとの人数、VCF にいないサンプル
cut -f2 groups.txt | tr ',' '\n' | sort | uniq -c
comm -13 <(bcftools query -l $VCF | sort) <(cut -f1 groups.txt | sort) | head
```

- グループ名は英数字と `_` だけにしてください (`NS6000`, `NSX`, `NS6000_ctrl` など)。
- サンプル名は 3 文字以上が必要です (bcftools 1.19 の制約)。
- 量的形質の場合は `_ctrl` のグループは不要です (Step 5 で `--hwe-groups` を省略)。
- **患者・対照の割合がシーケンサー間で大きく偏っていないか**をこの時点で確認してください。偏っていると、残ったアーチファクトがそのまま関連シグナルになりえます。

### マスク用 BED

MHC・セントロメア周辺・segmental duplication・GIAB の difficult regions を除外します。

```bash
mkdir -p mask
# GIAB difficult regions (GRCh38, v3.5)
curl -L -o mask/GRCh38_alldifficultregions.bed.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz"
# セントロメア (UCSC hg38) と segmental duplication
curl -sL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz | zcat | cut -f2-4 > mask/centromeres.bed
curl -sL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz | zcat | cut -f2-4 > mask/segdups.bed
# MHC (GRCh38)
printf 'chr6\t28510119\t33480577\n' > mask/mhc.bed
```

GIAB の difficult regions はセントロメアや segdup の多くを含みますが、重複して指定しても問題ありません。

## Step 3. 深さの確認と遺伝型マスク条件の決定

遺伝型マスクの既定値 `FMT/GQ<20 | FMT/DP<10 | FMT/DP>65` の上限 `DP>65` は、平均約 33x のデータを想定しています。
**バッチごとの平均の深さ**を確認し、必要なら調整します。

```bash
# chr20 の一部 (5 Mb) でバッチ別の平均 DP
for b in NS6000 NSX; do
  bcftools query -r chr20:30000000-35000000 -i 'FILTER="PASS"' \
    -S <(awk -v b=$b -F'\t' '{split($2,g,","); if (g[1]==b) print $1}' groups.txt) \
    -f '[%DP\n]' $VCF | awk -v b=$b '$1!="." {s+=$1; n++} END {printf "%s mean DP = %.1f\n", b, s/n}'
done
```

| 結果 | 対応 |
|---|---|
| 両バッチとも 30x 前後 | 既定のまま |
| 両バッチとも深い/浅い (例 45x) | 上限を平均の約 2 倍に: `--gt-filter 'FMT/GQ<20 \| FMT/DP<10 \| FMT/DP>90'` |
| **バッチ間で大きく違う** (例 25x と 40x) | DP の閾値だけで欠測率に差が出て、多くのサイトが `BatchDiffMiss` で落ちる。`--gt-filter 'FMT/GQ<20'` のように GQ だけで切ることを検討 |

以降の Step 4・5・7 では同じ `--gt-filter` を指定してください (既定値を使う場合は指定不要)。
ヘテロ遺伝型のアレル比の検定 (`--het-binom 1e-5`) も既定で行われます。

## Step 4. QC 前の診断

```bash
python3 $TOOLS/batch_af_test_vcf.py -i $VCF -S groups.txt --batches NS6000,NSX \
    --hwe-groups NS6000_ctrl,NSX_ctrl --sex sex.txt -o qc/before --jobs 16
less qc/before.summary.txt
```

入力はフィルタリング前なので、遺伝型マスクはここで既定どおり行われます (`--no-gt-filter` は付けない)。

`qc/before.summary.txt` の表の見方:

| 列 | 意味 |
|---|---|
| `qc_pass` | `filter_gwas.py` の基準をすべて通るサイトの割合 |
| `miss` / `dmiss` | 欠測率 ≥2% / バッチ間の欠測率の差 ≥1% のサイトの割合 |
| `batch_af(all)` | バッチ間の AF 差の検定 (P<1e-8) で引っかかるサイトの割合 (診断用) |
| `batch_af(qc_pass)` | 同じく、QC を通るサイトの中での割合。フィルタ後に残る AF 差 |
| `lambda(all)` / `lambda(qc_pass)` | 頻度帯ごとの λGC。両者の差がフィルタで除けるバッチ差 |

シーケンサー混在データでは、典型的に次のような傾向が見えます。

- `batch_af(all)` が数 % あるのに `batch_af(qc_pass)` はほぼ 0。AF 差のあるサイトの大半は、片方のバッチで欠測が多いサイト
- 片方のバッチだけで欠測が多い (`dmiss` が高い)。GC の多い染色体 (chr16, 17, 19, 22) で目立つ
- `lambda(qc_pass)` も 1 を少し超えたまま (1.1-1.2 など) 残る。祖先構成のわずかな違いの可能性があり、サイトを削らず GWAS のモデルで扱う

この段階の `qc/before.tsv.gz` は、Step 9 の有意サイトの事後確認に使うので残しておいてください。

## Step 5. サイトフィルタ (バッチエフェクトを除いた VCF の作成)

```bash
python3 $TOOLS/filter_gwas.py -i $VCF -S groups.txt --batches NS6000,NSX \
    --hwe-groups NS6000_ctrl,NSX_ctrl --sex sex.txt \
    --mask mask/GRCh38_alldifficultregions.bed.gz --mask mask/centromeres.bed \
    --mask mask/segdups.bed --mask mask/mhc.bed \
    --split --ref $REF --keep-annotated \
    -o qc/gwasqc --jobs 16
```

| オプション | 理由 |
|---|---|
| `--batches NS6000,NSX` | シーケンサーをバッチとして扱う |
| `--hwe-groups NS6000_ctrl,NSX_ctrl` | 2 値形質なので HWE は各バッチの対照群で判定 (量的形質なら省略) |
| `--sex sex.txt` | chrX も含める (PAR 以外は女性のみで HWE を判定し、男性ヘテロ割合で `MaleHet` を判定) |
| `--mask ...` | 難しい領域を除外 |
| `--split --ref` | 多アレルを分割して左寄せ |
| `--keep-annotated` | 全サイト + 除外理由のラベルも出力 (確認用) |
| (`--hard-filter` なし) | VQSR の判定 (FILTER) をそのまま使う |
| (`--no-gt-filter` なし) | 入力はマスク前なので、既定の遺伝型マスクを行う |

出力:

- **`qc/gwasqc.pass.bcf`**: バッチエフェクトを除いた GWAS 用の VCF (全基準を通過したサイト)
- `qc/gwasqc.annotated.bcf`: 全サイトと除外理由
- `qc/gwasqc.filter_summary.tsv`: 設定と、除外理由ごとのサイト数

```bash
# 除外理由の内訳
grep -E '^(input|label)' qc/gwasqc.filter_summary.tsv | column -t
# 特定のサイトがなぜ除外されたか
bcftools query -r chr1:1000000-1100000 -f '%CHROM:%POS %REF>%ALT %FILTER\n' qc/gwasqc.annotated.bcf | head
```

ラベルの意味 (詳しくは README):
`VQSRTranche...` (VQSR で除外) / `Monomorphic` / `BatchMiss` / `BatchDiffMiss` / `HWE` / `ExcHet` / `MaleHet` / `Mask`。
**バッチ間の AF 差による除外は一切していません。**

## Step 6. QC 後の診断

```bash
python3 $TOOLS/batch_af_test_vcf.py -i qc/gwasqc.pass.bcf -S groups.txt --batches NS6000,NSX \
    --hwe-groups NS6000_ctrl,NSX_ctrl --sex sex.txt --no-gt-filter -o qc/after --jobs 16
less qc/after.summary.txt
```

Step 5 で遺伝型はマスク済みなので `--no-gt-filter` を付けます。確認するポイント:

- `batch_af(all)` (全サイト = QC 後のサイト) が Step 4 の `batch_af(all)` より大きく下がっているか
- 残った AF 差のあるサイトが特定の領域に固まっていないか (`qc/after.tsv.gz` の `flag_batch_af` を染色体・位置で集計)。固まっていれば、その領域をマスクに追加して Step 5 をやり直す
- `lambda` が 1 に近づいたか。1.1-1.2 程度の上昇が頻度帯によらず均一に残る場合は、サイトを削らず Step 8・9 で扱う

## Step 7. サンプル QC

```bash
python3 $TOOLS/sample_qc.py -i qc/gwasqc.pass.bcf -S groups.txt --batches NS6000,NSX \
    --sex sex.txt --no-gt-filter -o qc/sampleqc --jobs 16
less qc/sampleqc.summary.txt
```

確認するポイント:

| 項目 | 見方 |
|---|---|
| `n_rare_nonref` / `rare_titv` / `n_singleton` のバッチ差 (`kruskal_p`) | 片方のシーケンサーにだけ出るレアな偽陽性は、欠測や HWE に兆候が出ないため Step 5 では除けない。差が大きければ、レアバリアント解析 (burden 検定など) ではバッチを共変量に入れたうえで影響を個別に評価する |
| 外れ値 (`outlier` 列) | 欠測率が高い、ヘテロが多すぎる (汚染の疑い)、変異数が極端など |
| 性別の不一致 (`sex_mismatch`) | 申告された性別と chrX のヘテロ比が合わない。取り違えの可能性 |

除外するサンプルを決めたら、**VCF から外して Step 5 からやり直します** (サイトの統計が変わるため)。

```bash
# 除外候補の一覧 (外れ値か性別の不一致があるサンプル)。内容を確認して除外するか決める
python3 - <<'EOF'
import pandas as pd
d = pd.read_csv('qc/sampleqc.samples.tsv', sep='\t')
flag = d.outlier.notna() | (d.sex_mismatch if 'sex_mismatch' in d else False)
cols = ['sample', 'batch', 'sex', 'outlier'] + (['x_sex'] if 'x_sex' in d else [])
print(d.loc[flag, cols].to_string(index=False))
EOF

# 除外するサンプルを exclude.txt に 1 行 1 サンプルで書いたら:
bcftools view -S ^exclude.txt --force-samples -Ob -o joint.vqsr.sqc.bcf $VCF && bcftools index joint.vqsr.sqc.bcf
VCF=joint.vqsr.sqc.bcf   # Step 5 (と Step 6, 7) を再実行
```

## Step 8. PCA でバッチの分離が消えたか確認

```bash
# plink2 形式に変換 (性別を与え、chrX の PAR を分ける)
awk 'BEGIN{print "#IID\tSEX"} {print}' sex.txt > qc/sex.plink.txt
plink2 --bcf qc/gwasqc.pass.bcf --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 truncate \
    --vcf-half-call m --max-alleles 2 --split-par hg38 --update-sex qc/sex.plink.txt \
    --memory 16000 --threads 16 --make-pgen --out qc/gwasqc

# common SNV で LD pruning → PCA
plink2 --pfile qc/gwasqc --autosome --maf 0.05 --indep-pairwise 1000 100 0.1 --memory 16000 --out qc/prune
plink2 --pfile qc/gwasqc --extract qc/prune.prune.in --pca 10 --memory 16000 --out qc/pca
```

`--memory` は必ず指定してください (plink2 は既定で物理メモリの半分を確保しようとします)。

各 PC がシーケンサーをどれだけ分けるかを数値で確認します (AUC 0.5 = 分離なし、1.0 = 完全に分離)。

```bash
cat > qc/pc_batch.py <<'EOF'
import sys, pandas as pd
from scipy.stats import mannwhitneyu
pcs = pd.read_csv(sys.argv[1], sep='\t').rename(columns={'#IID': 'IID'})
grp = pd.read_csv(sys.argv[2], sep='\t', header=None, names=['IID', 'groups'])
grp['batch'] = grp.groups.str.split(',').str[0]
d = pcs.merge(grp[['IID', 'batch']], on='IID')
a, b = sorted(d.batch.unique())[:2]
for pc in [c for c in d.columns if c.startswith('PC')]:
    x, y = d.loc[d.batch == a, pc], d.loc[d.batch == b, pc]
    u = mannwhitneyu(x, y)
    auc = u.statistic / (len(x) * len(y))
    print(f'{pc}\tAUC={max(auc, 1 - auc):.3f}\tp={u.pvalue:.2e}')
EOF
python3 qc/pc_batch.py qc/pca.eigenvec groups.txt
```

QC 前と比べたい場合は、`qc/gwasqc.pass.bcf` の代わりに VQSR の PASS だけを取り出した VCF
(`bcftools view -f PASS $VCF`) で同じ処理をして、AUC を並べてください。

| 結果 | 解釈と対応 |
|---|---|
| 上位の PC で AUC が 0.5 に近い | シーケンサーによる分離は解消 |
| AUC が高い PC が残る | (1) その PC の寄与 (ローディング) が大きいサイトの位置と性質を調べ、特定の領域に集中していればマスクに追加して Step 5 から再実行。(2) 2 群の祖先構成が実際に違う可能性があるので、1000 Genomes などの参照集団に投影した PCA で確認 |

## Step 9. GWAS と有意サイトの事後確認

- `qc/gwasqc.pgen` を使い、**シーケンサー (バッチ) を共変量に入れ**、Step 8 の PC と混合モデル (REGENIE, SAIGE など) で解析します。
- 患者と対照が両方のシーケンサーにいる場合は、シーケンサーごとに GWAS をしてメタ解析し、バッチ間の異質性 (Cochran の Q, I²) を見る方法も有力です。

有意になったサイトは、Step 4 の診断結果と照合してバッチ差の兆候が無いか確認します。

```bash
# hits.tsv: 有意サイトの一覧 (列 CHROM, POS)
python3 - <<'EOF'
import pandas as pd
hits = pd.read_csv('hits.tsv', sep='\t')
diag = pd.read_csv('qc/before.tsv.gz', sep='\t',
                   usecols=['CHROM', 'POS', 'REF', 'ALT', 'AF_NS6000', 'AF_NSX', 'MISS_NS6000', 'MISS_NSX', 'P_batch', 'qc_fail'])
print(hits.merge(diag, on=['CHROM', 'POS'], how='left').to_string(index=False))
EOF
```

`P_batch` が極端に小さい、片方のバッチで欠測が多い、などのサイトは、IGV でリードを確認するか、シーケンサーごとの効果の大きさが一致するかを見てください。
**AF 差は事前の除外ではなく、この事後確認にだけ使います。**

---

## 大規模データでの実行

- 並列数は `--jobs` で指定します。`filter_gwas.py` は 1 チャンクあたり bcftools のプロセスが 10 個前後つながって動くので、CPU コア数の 1/8〜1/4 程度が目安です。
- 一時ファイルは出力先と同じ場所に作られます (入力のおおむね 2 倍の容量)。別の場所にするには `--tmpdir` を指定します。
- 領域ごとに分割された VCF は `-I chunks.list` でそのまま渡せます。ジョブスケジューラで 1 チャンクずつ流す場合は `-i chunk.bcf --no-chunk` で処理し、最後に `bcftools concat -a -d exact` で結合してください。

## よくあるエラーと警告

| メッセージ | 対応 |
|---|---|
| `2 文字以下のサンプル名を扱えません` | bcftools 1.19 の制約。サンプル名を変更する (`bcftools reheader -s`) |
| `index not found` | `bcftools index $VCF` を実行するか、`--no-chunk` を付ける |
| `chrX を対象にするには --sex が必要です` | `--sex sex.txt` を指定するか、`--chroms` から chrX を外す |
| `chrX の長さから参照配列を判定できません` | `--par 10001-2781479,155701383-156030895` (GRCh38) のように PAR を指定 |
| `FORMAT/AD が無いため --het-binom を省略します` | AD が無い VCF。ヘテロのアレル比のマスクは行われない |
| `WARNING: --no-gt-filter ですが、...マスクされていない可能性があります` | マスク前の VCF に `--no-gt-filter` を付けている。外して実行する |
| `グループ ... に該当するサンプルが VCF にありません` | groups.txt のサンプル名やグループ名の綴りを確認 |
