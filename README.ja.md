# wgs-gwasqc

[![tests](https://github.com/nighpg/wgs-gwasqc/actions/workflows/tests.yml/badge.svg)](https://github.com/nighpg/wgs-gwasqc/actions/workflows/tests.yml) [![license](https://img.shields.io/github/license/nighpg/wgs-gwasqc)](LICENSE)

[English](README.md) | **日本語**

joint call 後の VCF から、**GWAS に使うサイトを選ぶ**ためのツール群です。
シーケンサー・ライブラリ調製ロット・委託先などが混在するコホートで生じるバッチエフェクトを、
**バッチ間のアレル頻度 (AF) 差を使わずに**軽減することを目的にしています。

AF 差でサイトを除外すると、表現型や祖先と関連する本物のシグナルまで消えてしまいます。
そこで、遺伝型の品質 (マスク後の欠測率とその差、バッチ内の HWE / ヘテロ過剰) と領域マスクでサイトを選び、
残るバッチ差は GWAS のモデル (バッチ共変量・PC・混合モデル) で扱う、という方針をとります。

| スクリプト | 役割 | 依存 |
|---|---|---|
| `filter_gwas.py` | **サイトフィルタ**。各基準を FILTER ラベルとして付け、全基準を通ったサイトを出力する | bcftools, Python 標準ライブラリ |
| `batch_af_test_vcf.py` | **サイト単位の診断**。バッチ間の AF 差・欠測差・HWE と、頻度帯ごとの λGC | bcftools, numpy, pandas, scipy |
| `sample_qc.py` | **サンプル単位の診断**。1 人あたりの変異数・Ti/Tv・ヘテロ/ホモ比・singleton などをバッチ間で比較、外れ値と性別不一致を検出 | bcftools, numpy, pandas, scipy |
| `gwasqc_common.py` | 上の 3 つが共通で使うモジュール (同じディレクトリに置く) | — |

## 動作環境

| 項目 | 検証した版 | 備考 |
|---|---|---|
| OS | Ubuntu 24.04 (Linux x86_64) | macOS は未検証 |
| bcftools / htslib | 1.19 | 下限の版は未確認 |
| Python | 3.12 | |
| numpy / pandas / scipy | 1.26 / 2.1 / 1.11 | `batch_af_test_vcf.py`, `sample_qc.py` のみ |

インストールは不要です。リポジトリを clone し、スクリプトを直接実行してください (`bcftools` に PATH を通しておく)。

## 基本の流れ

```
joint call 後の VCF
 ├─ ① 診断     batch_af_test_vcf.py   バッチ差の大きさ (λ, 欠測差, HWE) を把握
 ├─ ② フィルタ  filter_gwas.py         → <prefix>.pass.bcf
 ├─ ③ 再診断   batch_af_test_vcf.py   ②の出力に --no-gt-filter。QC 合格サイトに残る差を確認
 ├─ ④ サンプル QC  sample_qc.py       ②の出力に --no-gt-filter。レアバリアントのバッチ差・外れ値・性別不一致
 ├─ ⑤ plink2 に変換 → LD pruning → PCA → GWAS (バッチを共変量に。REGENIE / SAIGE など)
 └─ ⑥ 有意サイトの事後確認: ①の tsv (P_batch, dMISS, flag_*) と照合
```

```bash
# ① 診断
python3 batch_af_test_vcf.py -i joint.bcf -S groups.txt --batches BatchA,BatchB -o qc/before --jobs 8
# ② フィルタ
python3 filter_gwas.py -i joint.bcf -S groups.txt --batches BatchA,BatchB -o qc/gwasqc \
    --hwe-groups BatchA_ctrl,BatchB_ctrl --sex sex.txt \
    --mask mhc_centromere_segdup.bed.gz --mask GRCh38_alldifficultregions.bed.gz \
    --split --ref GRCh38.fa --hard-filter --jobs 16
# ③ 再診断 / ④ サンプル QC (②で遺伝型はマスク済みなので --no-gt-filter)
python3 batch_af_test_vcf.py -i qc/gwasqc.pass.bcf -S groups.txt --batches BatchA,BatchB --no-gt-filter -o qc/after
python3 sample_qc.py -i qc/gwasqc.pass.bcf -S groups.txt --batches BatchA,BatchB --sex sex.txt --no-gt-filter -o qc/sampleqc
```

## 入力

### VCF / BCF

| 項目 | 要件 |
|---|---|
| 種類 | joint genotyping 済みの複数サンプル VCF (GATK GenotypeGVCFs の出力など)。gVCF は不可 |
| FORMAT | GT 必須。GQ・DP・AD があれば遺伝型マスクに使う (無い条件は警告を出して省略) |
| `##contig` | ID と length が必要 (対象染色体の判定、分割、chrX の PAR 判定に使う) |
| 染色体 | 既定は常染色体 (`chr1`-`chr22` / `1`-`22` を自動判定)。`--sex` を指定すると chrX も。chrY・ミトコンドリアは対象外 |
| 多アレル | `filter_gwas.py --split` で分割するのを推奨。診断の 2 スクリプトは 2 アレルのサイトのみ扱う |
| FILTER | `PASS` または `.`。それ以外 (VQSR で落ちたサイトなど) は `filter_gwas.py` で除外される |
| INFO | `--hard-filter` を使うときは QD・FS・SOR・MQ・MQRankSum・ReadPosRankSum (GATK 標準)。無い項目は省略 |

入力の渡し方は 3 スクリプト共通です。

| 指定 | 内容 | index |
|---|---|---|
| `-i file` | 1 ファイルを `--chunk-mb` (既定 20 Mb) ごとに分割して `--jobs` 並列で処理 | 必要 |
| `-i file --no-chunk` | 1 ファイルを分割せずそのまま処理 (分割済みチャンクをジョブ配列で 1 つずつ流す場合など) | 不要 |
| `-I list.txt` | 分割済みチャンクのパス一覧 (1 行 1 ファイル、順序は任意)。各ファイルを 1 単位として並列処理 | 不要 |

`-I` ではチャンクを先頭の変異の位置でゲノム順に並べ直し、隣のチャンクと範囲が重なる record (境界をまたぐ変異など) を後ろ側から除去します。全チャンクのサンプル名と並び順は同じである必要があります。
`--no-chunk` で別々に処理した出力を結合するときは、`bcftools concat -a -d exact` で境界の重複を除いてください。

### groups.txt (`-S`)

`bcftools +fill-tags -S` と同じ形式です。

```
S0001	BatchA,BatchA_ctrl
S0002	BatchA
S0003	BatchB,BatchB_ctrl
```

- 1 列目: サンプル名 (VCF と完全一致)。**3 文字以上** (bcftools 1.19 の `+fill-tags` は 2 文字以下の名前でエラーになるため、事前に確認して止める)
- 2 列目: グループ名 (カンマ区切りで複数可)。英数字と `_` のみ (INFO タグ名の一部になる)。末尾 `__F` / `__M` は予約
- 1 人のサンプルが属するバッチは 1 つだけ。HWE 判定用のグループ (例: 対照群) は追加で付けてよい
- `--batches` に比べるバッチ、`--hwe-groups` に HWE を判定するグループを指定する
- `filter_gwas.py` と `sample_qc.py` は `--batches` を省略すると全サンプルを 1 バッチ (`COHORT`) として扱う (このとき `-S` も不要)

### sex.txt (`--sex`, chrX を含めるとき)

```
S0001	F
S0002	M
```

2 列目は `1`/`M`/`male` か `2`/`F`/`female`。それ以外の値や記載の無いサンプルは「性別不明」として chrX の女性・男性グループに含めません。

### マスク用 BED (`--mask`, 任意・複数可)

MHC、セントロメア周辺、segdup、GIAB difficult regions などの除外領域。3 列以上、`.gz` 可。染色体名を VCF と揃えてください。

## 遺伝型マスク (3 スクリプト共通、既定で有効)

欠測率によるフィルタは、**低品質の遺伝型が欠測 (`./.`) になっていること**を前提にしています。
GATK の joint call では低品質の遺伝型も `0/0` や `0/1` のまま残るため、既定で次のマスクを行います。

| オプション | 既定 | 内容 |
|---|---|---|
| `--gt-filter` | `FMT/GQ<20 \| FMT/DP<10 \| FMT/DP>65` | 当てはまる遺伝型を欠測に (`bcftools +setGT -t q -n .`) |
| `--het-binom` | `1e-5` | ヘテロ遺伝型で、AD のアレル比の二項検定 p がこれ未満なら欠測に (`+setGT -t "b:AD<P"`)。0 で無効 |
| `--no-gt-filter` | — | マスクしない (事前にマスク済みの入力向け)。マスクされていないように見える入力には警告を出す |

- 同じ条件で二重にマスクしても結果は変わりません。迷ったら既定のままで構いません。
- **`DP>65` は平均約 33x のデータで決めた上限**です (平均のおよそ 2 倍)。深さの違うデータでは `--gt-filter` で調整してください。
- バッチ間で平均の深さが大きく違うと、DP の閾値だけで欠測率に差が生じ、`BatchDiffMiss` で多くのサイトが落ちます。その場合は GQ だけで切ることも検討してください。
- `--max-miss` / `--max-dmiss` の既定値は、この既定のマスク条件で試算した値です。

## filter_gwas.py

全基準を通過したサイトを `<prefix>.pass.bcf` に出力します。どの基準にも AF 差は使いません。

| ラベル | 基準 | 関連オプション (既定) |
|---|---|---|
| (既存の FILTER) | 入力の FILTER が PASS / `.` 以外 | — |
| `HardFilter` | GATK 推奨ハードフィルタ (SNP: QD<2, QUAL<30, SOR>3, FS>60, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8 / INDEL: QD<2, QUAL<30, FS>200, ReadPosRankSum<-20, SOR>10) | `--hard-filter` (無効) |
| `Monomorphic` | 全アレルの AC 合計が 0 または AN (マスク後に単型) | — |
| `LowMAC` | min(AC, AN-AC) < MAC | `--min-mac` (無効) |
| `BatchMiss` | いずれかのバッチで欠測率 ≥ 閾値 | `--max-miss` (0.02) |
| `BatchDiffMiss` | いずれか 2 バッチ間で欠測率の差 ≥ 閾値 (バッチが 2 つ以上のとき) | `--max-dmiss` (0.01) |
| `HWE` | HWE 判定グループ内で p < 閾値 (chrX の PAR 以外は女性のみ) | `--hwe` (1e-6), `--hwe-groups` |
| `ExcHet` | いずれかのバッチ内でヘテロ過剰の検定 p < 閾値 (chrX の PAR 以外は女性のみ) | `--exchet` (1e-6) |
| `MaleHet` | chrX の PAR 以外で、いずれかのバッチの男性のヘテロ割合 > 閾値 | `--max-male-het` (0.01) |
| `Mask` | マスク用 BED と重なる | `--mask` |

- 欠測率は「そのバッチで遺伝型が呼ばれた人の割合」(`1 - NS/人数`) です。
- 2 値形質では `--hwe-groups` に各バッチの**対照群**を指定してください (本物の関連があると患者群で HWE が崩れるため)。
- 前処理: `--split` (多アレルの分割、`--ref` で左寄せ)、`ALT="*"` の除去、遺伝型マスク、`+fill-tags` によるバッチ別タグの付与。

出力:

| ファイル | 内容 |
|---|---|
| `<prefix>.pass.bcf` (+ `.csi`) | 全基準を通過したサイト。INFO にバッチ別の AN/AC/AF/NS/HWE/ExcHet (chrX を含めると男女別も) |
| `<prefix>.annotated.bcf` | `--keep-annotated` 指定時。全サイト + FILTER ラベル (除外理由の確認用) |
| `<prefix>.filter_summary.tsv` | 使った設定 (`setting` 行)、ラベル別・組み合わせ別のサイト数 |

## batch_af_test_vcf.py

2 つ以上のバッチを比べる、サイト単位の診断です。**`flag_batch_af` (AF 差の検定) は診断専用**で、GWAS のサイトフィルタには使わないでください。

- バッチ間の検定: 2 バッチは 2×2 アレル数表のカイ二乗 (期待度数が小さいサイトは Fisher)、3 バッチ以上は 2×k カイ二乗 (自由度 k-1)
- 判定フラグ `flag_mono` / `flag_miss` / `flag_dmiss` / `flag_hwe` / `flag_exchet` / `flag_malehet` と `qc_fail` は、`filter_gwas.py` のラベルと同じ基準
- 主なオプション: `--batches` (必須, 2 つ以上)、`--hwe-groups`、`--sex`、閾値は `filter_gwas.py` と同じ既定値、`--p-batch` (1e-8, 診断用)

出力:

| ファイル | 内容 |
|---|---|
| `<prefix>.tsv.gz` | 全サイトの AC/AN/NS/AF/欠測率 (バッチ別)、ΔAF・Δ欠測率 (3 バッチ以上は最大-最小)、MAF、P_batch、HWE/ExcHet、判定フラグ |
| `<prefix>.summary.txt` | 領域 (常染色体 / chrX) × 変異の種類 × 頻度帯ごとの QC 不合格率、AF 差のあるサイトの割合 (全サイト / QC 合格サイト)、λGC、\|ΔAF\| の 99% 点 |

結果の見方:

- **λGC は頻度帯ごとに見る**。全サイトの λ はレアバリアントが大半を占めるため意味をなしません (MAF<0.1% は省略しています)。
- `lambda(all)` と `lambda(qc_pass)` の差が、フィルタで除けるバッチ差です。QC 合格サイトに残る λ の上昇は、祖先構成の違いや系統的でない差の可能性があり、サイトを削るのではなく GWAS のモデルで扱います。参照集団に投影した PCA で祖先構成を確認してください。
- λ1000 (1,000 人 vs 1,000 人に換算した λ, 2 バッチのときのみ) は、各群が 1,000 人未満だと大きく出ます。小規模なデータでは λ をそのまま見てください。

## sample_qc.py

サンプル単位の指標をバッチ間で比べます。`filter_gwas.py` の出力に対して実行するのを想定しています。

| 指標 | 内容 |
|---|---|
| `missing_rate` | 欠測率 |
| `n_nonref_snv` | 非参照の SNV 遺伝型の数 (ヘテロ + 非参照ホモ) |
| `het_hom` | ヘテロ / 非参照ホモ 比 |
| `titv` | Ti/Tv |
| `n_indel` | インデルの非参照遺伝型の数 |
| `n_singleton` | コホート全体で 1 人だけが持つ変異の数 |
| `mean_dp` | 平均深さ (FORMAT/DP がある場合) |
| `n_rare_nonref`, `rare_titv` | MAF < `--rare-maf` (既定 0.01) のサイトでの非参照 SNV 数と Ti/Tv |
| `x_het_ratio` | `--sex` 指定時。chrX の PAR 以外でのヘテロ / (ヘテロ + 非参照ホモ)。申告された性別との不一致を検出 |

- 外れ値: バッチ内で頑健な z (中央値と MAD、MAD=0 なら平均絶対偏差) を計算し、|z| > `--outlier-z` (既定 4) のサンプルを報告します (`rare_titv` は 1 人あたりのレア変異が少ないと不安定なため判定に使いません)。
- バッチ比較: 指標ごとにバッチ別の中央値・四分位・全体比と Kruskal-Wallis 検定の p。

出力: `<prefix>.samples.tsv` (サンプルごと)、`<prefix>.batches.tsv` (指標 × バッチ)、`<prefix>.summary.txt`。

**サイト単位のフィルタで拾えないバッチ差**: 片方のバッチにだけ出るレアな偽陽性は、欠測率や HWE に兆候が出ないため `filter_gwas.py` を通過します。
`sample_qc.py` の `n_rare_nonref` / `rare_titv` / `n_singleton` にバッチ差があれば、レアバリアント解析ではバッチ共変量に加えて、その影響を個別に評価してください。

## chrX

`--sex` を指定すると chrX を対象に含めます。

- PAR は `##contig` の chrX の長さから GRCh38 / GRCh37 を自動判定します (判定できない場合は `--par 10001-2781479,155701383-156030895` のように指定)。PAR は常染色体と同じ扱いです。
- PAR 以外では、HWE とヘテロ過剰を**女性のみ**で判定し、男性のヘテロ割合で `MaleHet` を判定します。
- 欠測率は「遺伝型が呼ばれた人の割合」なので、男性の遺伝型が 1 倍体・2 倍体のどちらで書かれていても同じ基準になります。
- 男性のヘテロ遺伝型そのものは VCF に残ります。plink2 などで読み込むときに、性別情報を与えて 1 倍体として扱ってください。

## ユースケースと制約

**そのまま使えるケース**

- 複数バッチ (装置・ロット・委託先の違い) の WGS コホート。3 バッチ以上も可
- 2 値形質 (`--hwe-groups` に対照群) / 量的形質
- 事前に遺伝型をマスク済みの VCF (`--no-gt-filter`)
- 領域ごとに分割された joint call の出力 (`-I` / `--no-chunk`)
- 1 バッチのコホートの GWAS 用サイト QC (`filter_gwas.py`, `sample_qc.py`)

**設定の調整が必要なケース**

- 平均の深さが約 30x と違う WGS: `--gt-filter` の DP 上限を調整
- 祖先が多様なコホート: groups.txt に「祖先 × バッチ」のグループを作り `--hwe-groups` に指定
- 血縁者・同一個体の重複: HWE 判定グループから外す
- GATK 以外の joint call (DRAGEN など): `--hard-filter` は使わず、各ツールの FILTER を使う
- WES・パネル: DP の上限や対象領域の扱いを調整 (未検証)

**対象外**

- chrY・ミトコンドリア、構造多型 (SV)・CNV
- 低深度 WGS (おおむね 10x 未満。遺伝型尤度とインピュテーションを使うべき)
- joint call 前 (gVCF・リード) のバッチ補正

**注意**

- 患者・対照の割合がバッチ間で大きく偏っている場合、残ったアーチファクトが関連シグナルになりえます。バッチごとの GWAS とメタ解析を検討してください。
- バッチ間の AF 差 (`batch_af_test_vcf.py` の `flag_batch_af` / `P_batch`) は、GWAS では**サイトの除外に使わず**、有意サイトの事後確認に使ってください。表現型の無い用途 (AF 参照パネルの作成など) では除外に使っても害は小さくなります。

## bcftools 1.19 についての注意

- `+fill-tags -S` は 2 文字以下のサンプル名でエラーになります (スクリプトが事前に確認します)。
- `F_MISSING` はグループ別に計算されません。欠測率は `NS` から計算しています。
- フィルタ式で否定 (`!`) が使えないため、chrX の PAR 以外の判定は範囲の補集合で書いています。

## テスト

合成データ (Python 標準ライブラリで生成、シード固定) によるテストを `tests/` に置いています。

```bash
python3 tests/run_tests.py -v                 # bcftools が必要。診断ツールのテストは numpy/pandas/scipy が無ければ省略
KEEP_TEST_DATA=1 python3 tests/run_tests.py   # 生成したデータと出力を残す
```

確認している内容: 各基準の判定 (仕込んだサイトごと)、AF 差を判定に使わないこと、入力モード間 (`-i` / `-I` / `--no-chunk`) の一致、
遺伝型マスクの冪等性と未マスク入力の警告、chrX の PAR 判定と女性のみの HWE、`batch_af_test_vcf.py` のフラグと `filter_gwas.py` のラベルの一致、
3 バッチの検定 (scipy `chi2_contingency` との一致)、`sample_qc.py` の集計値と `bcftools stats` の一致、仕込んだ異常サンプルの検出。
GitHub Actions (`.github/workflows/tests.yml`) で、Ubuntu 24.04 の bcftools 1.19 を使って自動実行します。

1 万人規模の全ゲノムでの処理時間とメモリ使用量は測定していません。

## ライセンス

[MIT License](LICENSE)
