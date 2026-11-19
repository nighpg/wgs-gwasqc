#!/usr/bin/env python3
"""batch_af_test_vcf.py の診断結果 (補正前/補正後) を染色体イデオグラム上に可視化する

入力: batch_af_test_vcf.py が出力する <prefix>.tsv.gz を2つ (補正前と補正後, 通常は
      補正後を --no-gt-filter で再診断したもの)。染色体ごとに1枚のPNGを出力する。

各PNGは次の構成 (上から):
  1. イデオグラム (UCSC cytoBandIdeo の G バンド + セントロメア + テロメア)
  2. (--mask 指定時のみ) filter_gwas.py の --mask と同じ BED で指定した除外領域を塗り領域として表示
  3. 指標ごとに Before/After の2行1組 (AF(+dAFを右軸に重ねる) / P_batch(診断専用, AF差) / dMISS /
     HWE p / ExcHet p): AF/HWE/ExcHet はグループ (--batches / --hwe-groups) ごとに色分けして
     同じ軸に重ねる (同じグループ名はどのトラックでも同じ色)

各トラックは --window-mb (既定 1bp = 集計なし、生のサイト単位の値をそのままプロット) ごとに集計した
生の値をそのまま使う (二次加工なし):
  AFはグループ (--batches で指定したシーケンサー・コホートなど) ごとの列 (AF_<group>) を窓内中央値に
  したものを、グループごとに別系列で点に打つ。ちょうど2グループ (通常のバッチエフェクト診断のケース) の
  ときは、2番目のグループだけ表示を反転し (y = -AF) ゼロ線をはさんで下向きに伸びるようにする (保存
  される値は変えず、表示だけ反転する)。3グループ以上のときは全グループを正の向きのまま重ねて表示する。
  AFの2点 (反転表示) だけだと、2点間のギャップは実は両グループのAFの合計を表してしまい、差の大小は
  中点のゼロからのズレにしか現れず分かりにくいため、dAF (差そのもの) を同じ行に、専用の目盛り
  (「dAF」ラベル) を持つ右側のY軸を使った折れ線として重ねて表示する (範囲はAFの軸と共有し、差の
  実際の大きさがそのまま比較できるようにする)。dAF/dMISS/P_batch は窓内の中央値 (符号維持。
  2バッチ超では符号付きのdAF列が無く
  dAF_range (|差|の範囲) になる)。HWE/ExcHet もグループ (--hwe-groups / --batches) ごとに色分けし、
  各グループの列 (HWE_<group>, ExcHet_<group>) を窓内最小値 (= そのグループで最も悪いサイト) にした
  ものを点に打つ。p 値 (HWE/ExcHet/P_batch) は Y軸だけ対数スケールにする (値自体は変換しない)。
  dAF/dMISS/P_batch はグループ別の生データが元の tsv に無いため、
  従来どおり全体をまとめた1系列のまま。

usage:
  plot_batch_effect.py --before qc/before.tsv.gz --after qc/after.tsv.gz -o qc/ideogram
"""
import argparse
import gzip
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')  # ヘッドレス環境 (Singularityコンテナ含む) で動かすため、GUIバックエンドを使わない
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from gwasqc_common import make_log  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
DEFAULT_CYTOBAND = os.path.join(HERE, 'data', 'cytoBandIdeo.hg38.txt.gz')
CHROM_ORDER = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
IDEOGRAM_COLORS = {
    'gneg': '#ffffff', 'gpos25': '#bdbdbd', 'gpos50': '#808080', 'gpos75': '#4d4d4d',
    'gpos100': '#000000', 'gvar': '#dcdcdc', 'stalk': '#eeeeee',
}
CENTROMERE_COLOR = '#b35959'
MASK_COLOR = 'crimson'
DAF_COLOR = 'magenta'
# (列キー, 表示名, スケール, 色, グループ別に色分けするか)。dMISS/p_batch は窓内の中央値
# (いずれも二次加工はせず生の値をそのまま使う, 符号維持)。
# af/hwe/exchet はグループ別の列 (AF_<group> / HWE_<group> / ExcHet_<group>) をそれぞれグループごとに
# 別色でプロットする (色は None にして、実行時にグループごとの色を割り当てる; af は窓内中央値,
# hwe/exchet は窓内最小値)。af はちょうど2グループのとき2番目を反転表示する (af_plot_spec 参照)。
# dAF は、AF の2点 (反転表示) だけでは差の大小が分かりにくい (ギャップは実は AF の合計を表し、差は
# 中点のズレにしか現れない) ため、AF と同じ行に、専用の右側Y軸 (範囲はAFと共有) を使った折れ線として
# 重ねて表示する。
TRACKS = [
    ('af', 'AF', 'linear', None, True),
    ('p_batch', 'P_batch p', 'log', 'gray', False),
    ('dmiss', 'dMISS', 'linear', 'steelblue', False),
    ('hwe', 'HWE p', 'log', None, True),
    ('exchet', 'ExcHet p', 'log', None, True),
]
GROUP_COLORMAP = 'tab10'


def _open_maybe_gzip(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_cytoband(path):
    """UCSC cytoBandIdeo (chrom start end label gieStain) を読み、染色体ごとの帯・セントロメア・長さを返す"""
    chroms = {}
    with _open_maybe_gzip(path) as f:
        for line in f:
            chrom, start, end, _label, stain = line.rstrip('\n').split('\t')
            if '_' in chrom:
                continue   # alt/random/unplaced contig は除外
            start, end = int(start), int(end)
            c = chroms.setdefault(chrom, {'bands': [], 'centromere': [], 'length': 0})
            c['length'] = max(c['length'], end)
            if stain == 'acen':
                c['centromere'].append((start, end - start))
            else:
                c['bands'].append((start, end - start, IDEOGRAM_COLORS.get(stain, '#ffffff')))
    for c in chroms.values():
        c['centromere'].sort()
    if not chroms:
        sys.exit(f'{path}: cytoBandIdeo のデータが読めません')
    return chroms


def read_mask_bed(paths):
    """--mask の BED (filter_gwas.py の --mask と同じ形式: chrom\\tstart\\tend, 0-based half-open, .gz可,
    #/track/browser で始まる行はスキップ) を読み、染色体ごとの (start, end) 区間一覧を返す。
    複数ファイル指定時はまとめて返す (区間が重なっていても描画で塗り重なるだけなので問題ない)。"""
    regions = {}
    for path in paths:
        with _open_maybe_gzip(path) as f:
            for line in f:
                if line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 3:
                    continue
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                regions.setdefault(chrom, []).append((start, end))
    return regions


def draw_mask_track(ax, regions):
    """--mask の BED 領域を染色体座標上に塗りバーで示す (イデオグラムの直下に1行追加する)"""
    xranges = [(s / 1e6, (e - s) / 1e6) for s, e in regions]
    if xranges:
        ax.broken_barh(xranges, (0.1, 0.8), facecolors=MASK_COLOR, edgecolors='none')
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel('Mask', rotation=0, ha='right', va='center', fontsize=8)


def _metric_columns(columns):
    af_cols = sorted(c for c in columns if re.match(r'^AF_.+$', c))
    hwe_cols = sorted(c for c in columns if re.match(r'^HWE_(?!all$).+$', c))
    exchet_cols = sorted(c for c in columns if re.match(r'^ExcHet_.+$', c))
    dmiss_col = 'dMISS' if 'dMISS' in columns else 'dMISS_range'
    daf_col = 'dAF' if 'dAF' in columns else 'dAF_range'
    return af_cols, hwe_cols, exchet_cols, dmiss_col, daf_col


def load_window_stats(path, window_bp, log):
    """batch_af_test_vcf.py の tsv.gz を読み、染色体 x 窓ごとに集計する (窓内は生の値をそのまま使う)。
    戻り値: (集計済み DataFrame, AF のグループ名一覧, HWE のグループ名一覧, ExcHet のグループ名一覧)。
    AF/HWE/ExcHet はグループ (AF_<group> / HWE_<group> / ExcHet_<group>) ごとに列
    'af__<group>' / 'hwe__<group>' / 'exchet__<group>' として残す
    (グループ別に色分けして描画するため、ここではグループ間でまとめない)。"""
    header = pd.read_csv(path, sep='\t', nrows=0).columns.tolist()
    for required in ('CHROM', 'POS', 'P_batch'):
        if required not in header:
            sys.exit(f'{path}: 列 {required} がありません (batch_af_test_vcf.py の出力ですか?)')
    af_cols, hwe_cols, exchet_cols, dmiss_col, daf_col = _metric_columns(header)
    if dmiss_col not in header:
        sys.exit(f'{path}: dMISS / dMISS_range 列がありません')
    if daf_col not in header:
        sys.exit(f'{path}: dAF / dAF_range 列がありません')
    af_groups = [c[len('AF_'):] for c in af_cols]
    hwe_groups = [c[len('HWE_'):] for c in hwe_cols]
    exchet_groups = [c[len('ExcHet_'):] for c in exchet_cols]
    usecols = ['CHROM', 'POS', 'P_batch', dmiss_col, daf_col] + af_cols + hwe_cols + exchet_cols
    log(f'{os.path.basename(path)}: 読み込み中 ({len(af_cols)} AF列, {len(hwe_cols)} HWE列, {len(exchet_cols)} ExcHet列)')
    df = pd.read_csv(path, sep='\t', usecols=usecols, dtype={'CHROM': str}, na_values='.')
    df['window'] = df['POS'].values // window_bp
    agg = {'daf': (daf_col, 'median'), 'dmiss': (dmiss_col, 'median'),
           'p_batch': ('P_batch', 'median'), 'n': ('POS', 'size')}
    agg.update({f'af__{g}': (c, 'median') for g, c in zip(af_groups, af_cols)})
    agg.update({f'hwe__{g}': (c, 'min') for g, c in zip(hwe_groups, hwe_cols)})
    agg.update({f'exchet__{g}': (c, 'min') for g, c in zip(exchet_groups, exchet_cols)})
    out = df.groupby(['CHROM', 'window'], sort=False).agg(**agg).reset_index()
    out['start'] = out['window'] * window_bp
    log(f'{os.path.basename(path)}: {len(df)} サイト, {len(out)} 窓')
    return out, af_groups, hwe_groups, exchet_groups


def group_colors(groups):
    """グループ名ごとに色を割り当てる (アルファベット順で決定論的に)"""
    cmap = plt.get_cmap(GROUP_COLORMAP)
    return {g: cmap(i % cmap.N) for i, g in enumerate(groups)}


def af_plot_spec(af_groups):
    """AF トラックで各グループをどちらの向きに表示するか (グループ名, 符号) を返す。
    ちょうど2グループ (バッチエフェクト診断の通常のケース) なら2番目のグループを負の向きにして、
    ゼロ線をはさんで下向きに伸びるようにする (dAF を別パネルにせず、AF の高さの違いとして直接比較
    できるようにするため)。保存される値自体は変えず、表示だけ反転する。2グループ以外ではすべて
    正の向きのまま表示する。"""
    if len(af_groups) == 2:
        return [(af_groups[0], 1), (af_groups[1], -1)]
    return [(g, 1) for g in af_groups]


def _calc_ylim(vals, scale):
    """値の列から (下限, 上限) を決める。線形は原点対称、対数は実際の最小/最大の周辺にする"""
    vals = vals.replace([np.inf, -np.inf], np.nan).dropna()
    if scale == 'log':
        vals = vals[vals > 0]
    if len(vals) == 0:
        return (0.0, 1.0)
    if scale == 'linear':
        m = max(abs(vals.min()), abs(vals.max()), 1e-9)
        return (-m * 1.1, m * 1.1)
    if scale == 'log':
        return (vals.min() * 0.5, vals.max() * 2)
    return (0.0, 1.0)


def shared_ylims(before, after, af_groups, hwe_groups, exchet_groups):
    """トラックごとに Before/After 共通のY軸範囲を決める (生の値の最小/最大をそのまま使う)。
    グループ別トラック (af/hwe/exchet) は全グループの値をまとめて範囲を決める。af は表示時の符号
    (af_plot_spec) を反映した値で範囲を決める。dAF は AF と重ねて描くため、同じ範囲 (lims['af']) を
    共有して使う (軸を分けて独立にスケールすると、dAF の実際の大きさがAFに対して誇張されて見えるため)。"""
    group_cols = {'hwe': [f'hwe__{g}' for g in hwe_groups], 'exchet': [f'exchet__{g}' for g in exchet_groups]}
    af_spec = af_plot_spec(af_groups)
    lims = {}
    for key, _, scale, _, grouped in TRACKS:
        if key == 'af':
            vals = (pd.concat([sign * before[f'af__{g}'] for g, sign in af_spec]
                               + [sign * after[f'af__{g}'] for g, sign in af_spec])
                    if af_spec else pd.Series(dtype=float))
        elif grouped:
            cols = group_cols[key]
            vals = pd.concat([before[c] for c in cols] + [after[c] for c in cols]) if cols else pd.Series(dtype=float)
        else:
            vals = pd.concat([before[key], after[key]])
        lims[key] = _calc_ylim(vals, scale)
    return lims


def draw_ideogram(ax, cyto, pad_mb):
    y0, height = 0.1, 0.8
    xranges = [(s / 1e6, w / 1e6) for s, w, _ in cyto['bands']]
    colors = [c for _, _, c in cyto['bands']]
    # matplotlib.collections.BrokenBarHCollection は 3.9 で削除されたため Axes.broken_barh を使う (同等の実装)
    ax.broken_barh(xranges, (y0, height), facecolors=colors, edgecolors='black', linewidths=0.3)
    cent = cyto['centromere']
    if len(cent) == 2:
        (s0, w0), (s1, w1) = cent
        x0, x1, x2 = s0 / 1e6, (s0 + w0) / 1e6, (s1 + w1) / 1e6
        mid = y0 + height / 2
        for tri in ([[x0, y0], [x0, y0 + height], [x1, mid]],
                    [[x1, mid], [x2, y0 + height], [x2, y0]]):
            ax.add_patch(mpatches.Polygon(tri, closed=True, facecolor=CENTROMERE_COLOR,
                                           edgecolor='black', linewidth=0.3))
    length_mb = cyto['length'] / 1e6
    # テロメア: 両端を半円でキャップする (座標範囲の外側にはみ出す分だけ軸をpad_mb広げる)
    for x, theta in ((0.0, (90, 270)), (length_mb, (270, 450))):
        ax.add_patch(mpatches.Wedge((x, y0 + height / 2), pad_mb, *theta, facecolor='white',
                                     edgecolor='black', linewidth=0.3, zorder=3))
    ax.set_ylim(0, 1)
    ax.set_yticks([])


def plot_chromosome(chrom, cyto, before, after, thresholds, window_mb, af_groups, hwe_groups, exchet_groups,
                     mask_regions, out_path, dpi):
    length_mb = cyto['length'] / 1e6
    pad_mb = min(0.6, length_mb * 0.004)
    has_mask = mask_regions is not None
    height_ratios = [1.4] + ([0.5] if has_mask else []) + [1] * (len(TRACKS) * 2)
    fig, axes = plt.subplots(len(height_ratios), 1, figsize=(12, 1.1 * sum(height_ratios)),
                              sharex=True, gridspec_kw={'height_ratios': height_ratios})
    draw_ideogram(axes[0], cyto, pad_mb)
    axes[0].set_ylabel(chrom, rotation=0, ha='right', va='center', fontsize=9)
    row0 = 1
    if has_mask:
        draw_mask_track(axes[1], mask_regions)
        row0 = 2

    ylims = shared_ylims(before, after, af_groups, hwe_groups, exchet_groups)
    group_lists = {'hwe': hwe_groups, 'exchet': exchet_groups}
    af_spec = af_plot_spec(af_groups)
    # 同じグループ名は AF/HWE/ExcHet のどのトラックでも同じ色にする (全グループ名の和集合で1つの色マップを作る)
    shared_colors = group_colors(sorted(set(af_groups) | set(hwe_groups) | set(exchet_groups)))
    # 染色体, 指標1 Before, 指標1 After, 指標2 Before, 指標2 After, ... の順に並べる
    for row, (key, title, scale, color, grouped) in enumerate(TRACKS):
        for block, (label, df) in enumerate((('Before', before), ('After', after))):
            ax = axes[row0 + row * 2 + block]
            x_mb = df['start'] / 1e6 + window_mb / 2   # 窓の中心に点を打つ
            if key == 'af':
                for g, sign in af_spec:
                    ax.plot(x_mb, sign * df[f'af__{g}'], color=shared_colors[g], marker='o', markersize=2.5,
                             markeredgewidth=0, linestyle='none', alpha=0.75, label=g)
                if len(af_spec) == 2:
                    # 反転表示している2番目のグループの目盛りが負の値のままだと分かりにくいので、
                    # 「<グループ名> AF <絶対値>」の形式にする (値自体は変えず、ラベルの見せ方だけ変更)
                    g1, g2 = af_spec[0][0], af_spec[1][0]

                    def _af_tick(v, _pos, g1=g1, g2=g2):
                        if abs(v) < 1e-9:
                            return '0'
                        return f'{g1} AF {v:.2f}' if v > 0 else f'{g2} AF {-v:.2f}'
                    ax.yaxis.set_major_formatter(FuncFormatter(_af_tick))
                # dAF を AF と同じ範囲 (ylims['af']) を共有する右側の専用Y軸に折れ線で重ねる。
                # 軸を分けて独立にスケールすると、dAF の実際の大きさがAFに対して誇張されて見えるため、
                # 目盛りの範囲はAFと統一し、差分としての実際の大小がそのまま見えるようにする
                ax2 = ax.twinx()
                ax2.plot(x_mb, df['daf'], color=DAF_COLOR, linewidth=1.0, linestyle='-', label='dAF')
                ax2.set_ylim(*ylims['af'])
                ax2.tick_params(labelsize=7, labelcolor=DAF_COLOR)
                ax2.set_ylabel('dAF', fontsize=6.5, color=DAF_COLOR)
                # AF (グループ) の凡例は軸の右上外側、dAF の凡例は軸の左上外側に分けて
                # データ点と重ならないようにする (どちらもプロット枠の外に配置)
                h1, l1 = ax.get_legend_handles_labels()
                leg1 = ax.legend(h1, l1, loc='lower right', bbox_to_anchor=(1.0, 1.0), fontsize=5,
                                  framealpha=0.6, markerscale=1.5, handletextpad=0.3, borderaxespad=0.1)
                ax.add_artist(leg1)
                h2, l2 = ax2.get_legend_handles_labels()
                ax2.legend(h2, l2, loc='lower left', bbox_to_anchor=(0.0, 1.0), fontsize=5, framealpha=0.6,
                           handletextpad=0.3, borderaxespad=0.1)
            elif grouped:
                for g in group_lists[key]:
                    ax.plot(x_mb, df[f'{key}__{g}'], color=shared_colors[g], marker='o', markersize=2.5,
                             markeredgewidth=0, linestyle='none', alpha=0.75, label=g)
            else:
                ax.plot(x_mb, df[key], color=color, marker='o', markersize=2.5, markeredgewidth=0, linestyle='none')
            if scale == 'log':
                ax.set_yscale('log')
                # 通常の目盛りも "1e-06" 形式にし、点線の脇の数値と表記を揃える (既定は "10^-6" 形式で混在する)
                ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f'{v:g}' if v > 0 else '0'))
            else:
                ax.axhline(0, color='black', linewidth=0.5, linestyle=':')
            ax.set_ylim(*ylims[key])
            thr = thresholds.get(key)
            if thr is not None:
                # 点線 = QCで使う基準線 (この値より悪い側がQC不合格になる境界)。
                # 凡例に加えて、線のすぐ左側にも値を直接表示する (凡例は軸の外側で離れているため)
                ax.axhline(thr, color='black', linewidth=0.6, linestyle='--', label=f'QC threshold ({thr:g})')
                ax.annotate(f'{thr:g}', xy=(0, thr), xycoords=('axes fraction', 'data'), xytext=(-4, 0),
                            textcoords='offset points', ha='right', va='center', fontsize=6, clip_on=False)
            if key != 'af' and ax.get_legend_handles_labels()[0]:
                # 枠内に重ならないよう、軸の右上外側に配置する
                ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.0), fontsize=5, framealpha=0.6,
                          markerscale=1.5, handletextpad=0.3, borderaxespad=0.1)
            ax.set_ylabel(label, fontsize=7)   # Before/After のみ (指標名はペアごとに大きくまとめて表示)
            ax.tick_params(labelsize=7)
    axes[-1].set_xlabel('Position (Mb)')
    axes[-1].set_xlim(-pad_mb, length_mb + pad_mb)
    window_bp = round(window_mb * 1_000_000)
    window_label = f'{window_bp} bp' if window_bp < 1_000_000 else f'{window_mb:g} Mb'
    fig.suptitle(f'{chrom} (window={window_label})', fontsize=11)
    fig.tight_layout(rect=(0.07, 0, 1, 0.98))   # 指標名ラベル用に左側へ余白を確保

    # Before/After 1組ごとに、指標名をまたがる大きなラベルを付ける (比較が1セットだと分かるように)
    for row, (key, title, scale, color, grouped) in enumerate(TRACKS):
        p_before = axes[row0 + row * 2].get_position()
        p_after = axes[row0 + row * 2 + 1].get_position()
        y_center = (p_before.y1 + p_after.y0) / 2
        fig.text(0.02, y_center, title, rotation=90, ha='center', va='center', fontsize=11, fontweight='bold')

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--before', required=True, help='補正前の batch_af_test_vcf.py 出力 (<prefix>.tsv.gz)')
    ap.add_argument('--after', required=True, help='補正後の batch_af_test_vcf.py 出力 (<prefix>.tsv.gz, 通常 --no-gt-filter)')
    ap.add_argument('-o', '--outdir', required=True, help='染色体ごとのPNGを書き出すディレクトリ')
    ap.add_argument('--window-mb', type=float, default=1e-6, help='集計する窓のサイズ (Mb, 既定 1e-6 = 1bp、つまり集計なし)')
    ap.add_argument('--cytoband', default=DEFAULT_CYTOBAND,
                     help=f'UCSC cytoBandIdeo ファイル (既定の配置先: {DEFAULT_CYTOBAND})')
    ap.add_argument('--chroms', help='対象染色体 (カンマ区切り)。既定: 両方の入力にデータがある染色体すべて')
    ap.add_argument('--mask', action='append', default=[],
                     help='イデオグラム直下に塗り領域として重ねる除外領域 BED (filter_gwas.py の --mask と'
                          '同じ形式, 複数指定可, .gz 可)。指定しなければマスク行自体を表示しない')
    ap.add_argument('--hwe', type=float, default=1e-6, help='HWEの基準線 (batch_af_test_vcf.py の既定と同じ 1e-6)')
    ap.add_argument('--exchet', type=float, default=1e-6, help='ExcHetの基準線 (既定 1e-6)')
    ap.add_argument('--p-batch', type=float, default=1e-8, help='P_batchの基準線 (既定 1e-8, 診断専用)')
    ap.add_argument('--dpi', type=int, default=150, help='PNGの解像度 (既定 150)')
    a = ap.parse_args()
    log = make_log('plot_batch_effect')

    window_bp = int(round(a.window_mb * 1_000_000))
    if window_bp <= 0:
        sys.exit('--window-mb は正の値にしてください')
    cyto = read_cytoband(a.cytoband)
    before, before_af_groups, before_hwe_groups, before_exchet_groups = load_window_stats(a.before, window_bp, log)
    after, after_af_groups, after_hwe_groups, after_exchet_groups = load_window_stats(a.after, window_bp, log)
    af_groups = sorted(set(before_af_groups) | set(after_af_groups))
    hwe_groups = sorted(set(before_hwe_groups) | set(after_hwe_groups))
    exchet_groups = sorted(set(before_exchet_groups) | set(after_exchet_groups))

    if a.chroms:
        chroms = [c.strip() for c in a.chroms.split(',') if c.strip()]
    else:
        present = set(before['CHROM']) | set(after['CHROM'])
        chroms = [c for c in CHROM_ORDER if c in present and c in cyto]
    if not chroms:
        sys.exit('描画対象の染色体がありません (--before / --after / --cytoband / --chroms を確認してください)')

    mask_by_chrom = read_mask_bed(a.mask) if a.mask else None

    thresholds = {'hwe': a.hwe, 'exchet': a.exchet, 'p_batch': a.p_batch}
    os.makedirs(a.outdir, exist_ok=True)
    empty_cols = {k: pd.Series(dtype=float) for k, _, _, _, grouped in TRACKS if not grouped}
    empty_cols['daf'] = pd.Series(dtype=float)   # af パネルに重ねる折れ線用 (TRACKS には無い列)
    for g in af_groups:
        empty_cols[f'af__{g}'] = pd.Series(dtype=float)
    for g in hwe_groups:
        empty_cols[f'hwe__{g}'] = pd.Series(dtype=float)
    for g in exchet_groups:
        empty_cols[f'exchet__{g}'] = pd.Series(dtype=float)
    empty_cols['start'] = pd.Series(dtype=float)
    rendered = 0
    for chrom in chroms:
        if chrom not in cyto:
            log(f'{chrom}: cytoBandIdeo にないためスキップ')
            continue
        b = before[before.CHROM == chrom].sort_values('window')
        af = after[after.CHROM == chrom].sort_values('window')
        if not len(b) and not len(af):
            log(f'{chrom}: Before/After いずれにもデータがないためスキップ')
            continue
        if not len(b):
            b = pd.DataFrame(empty_cols)
        if not len(af):
            af = pd.DataFrame(empty_cols)
        out_path = os.path.join(a.outdir, f'{chrom}.png')
        mask_regions = (mask_by_chrom.get(chrom, []) if mask_by_chrom is not None else None)
        plot_chromosome(chrom, cyto[chrom], b, af, thresholds, a.window_mb, af_groups, hwe_groups, exchet_groups,
                         mask_regions, out_path, a.dpi)
        rendered += 1
        log(f'{chrom}: {out_path}')
    if rendered == 0:
        sys.exit('描画対象の染色体がありません (--before / --after / --cytoband / --chroms を確認してください)')


if __name__ == '__main__':
    main()
