"""filter_gwas.py / batch_af_test_vcf.py / sample_qc.py の共通処理

- グループ (バッチ, HWE 判定グループ) と性別の読み込み・検証、+fill-tags -S 用グループファイルの生成
- 対象染色体と chrX の PAR の判定
- 遺伝型マスク (GQ/DP, ヘテロの AD 二項検定) の設定と、未マスク入力の検出
- 並列処理の単位 (chunk) の計画
"""
import collections
import os
import re
import shlex
import subprocess
import sys

q = shlex.quote
DEF_GT_FILTER = 'FMT/GQ<20 | FMT/DP<10 | FMT/DP>65'
SINGLE_BATCH = 'COHORT'           # --batches 省略時に全サンプルを入れるグループ名
FEMALE, MALE = '__F', '__M'       # chrX 用に自動で作るグループの接尾辞
# chrX の長さ → PAR (1-based, 両端含む)
PAR_BY_LENGTH = {156040895: ('GRCh38', [(10001, 2781479), (155701383, 156030895)]),
                 155270560: ('GRCh37', [(60001, 2699520), (154931044, 155260560)])}
TAG_RE = re.compile(r'[A-Za-z_][0-9A-Za-z_]*')


def run(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def make_log(tag):
    return lambda *x: print(f'[{tag}]', *x, file=sys.stderr, flush=True)


def add_common_args(ap, batches_help, need_two=False):
    """3 スクリプト共通の引数"""
    g = ap.add_argument_group('groups / sex')
    g.add_argument('-S', '--groups', help='sample<TAB>group[,group...] (bcftools +fill-tags -S 形式)。--batches / --hwe-groups の定義に使う')
    g.add_argument('--batches', required=need_two, help=batches_help)
    g.add_argument('--hwe-groups', help='HWE を判定するグループ (既定: バッチ)。2値形質ならバッチ別の対照群を指定')
    g.add_argument('--sex', help='sample<TAB>sex (1/M/male または 2/F/female)。指定すると chrX も対象にする')
    g.add_argument('--par', help='chrX の PAR を手動指定 (例 10001-2781479,155701383-156030895)。既定は contig 長から GRCh38/37 を判定')
    c = ap.add_argument_group('regions')
    c.add_argument('--chroms', help='対象染色体 (カンマ区切り)。既定: 常染色体 1-22 (--sex 指定時は X も)')
    m = ap.add_argument_group('genotype mask (既定で有効)')
    m.add_argument('--gt-filter', default=DEF_GT_FILTER, help=f"遺伝型を欠測にする条件 (既定 '{DEF_GT_FILTER}')")
    m.add_argument('--het-binom', type=float, default=1e-5,
                   help='ヘテロ遺伝型の AD アレル比の二項検定 p がこれ未満なら欠測に (既定 1e-5, 0 で無効)')
    m.add_argument('--no-gt-filter', action='store_true', help='遺伝型のマスクを行わない (事前にマスク済みの入力向け)')


def add_chunk_args(ap):
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('-i', '--input', help='VCF/BCF 1 ファイル (index があれば染色体/領域ごとに分割して並列処理)')
    src.add_argument('-I', '--input-list', help='分割済み VCF/BCF のパス一覧 (1 行 1 ファイル; 各ファイルを 1 chunk として処理)')
    ap.add_argument('--no-chunk', action='store_true', help='-i のファイルを分割せず 1 chunk として処理 (index 不要)')
    ap.add_argument('--chunk-mb', type=int, default=20, help='-i を分割するサイズ Mb (既定 20)')
    ap.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 2) // 4), help='並列ジョブ数')


# ---------------------------------------------------------------------------
def list_inputs(a):
    """-i / -I から入力ファイル一覧を作り、存在とサンプルの一致を確認"""
    if a.input_list:
        inputs = [l.strip() for l in open(a.input_list) if l.strip() and not l.startswith('#')]
        if not inputs:
            sys.exit(f'{a.input_list} にファイルがありません')
    else:
        inputs = [a.input]
    for f in inputs:
        if not os.path.exists(f):
            sys.exit(f'file not found: {f}')
    if len(inputs) > 1:
        import concurrent.futures as cf
        ref = run(['bcftools', 'query', '-l', inputs[0]])
        with cf.ThreadPoolExecutor(8) as ex:
            diff = [f for f, s in zip(inputs, ex.map(lambda f: run(['bcftools', 'query', '-l', f]), inputs)) if s != ref]
        if diff:
            sys.exit(f'サンプル (または順序) が 1 ファイル目と異なる chunk があります: {diff[:5]}')
    return inputs


def has_index(path):
    return os.path.exists(path + '.csi') or os.path.exists(path + '.tbi')


def read_sex(path, samples):
    sex = {}
    code = {'1': 'M', 'm': 'M', 'male': 'M', '2': 'F', 'f': 'F', 'female': 'F'}
    for line in open(path):
        f = line.split()
        if len(f) < 2 or f[0].startswith('#'):
            continue
        s = code.get(f[1].lower())
        if s and f[0] in samples:
            sex[f[0]] = s
    return sex


class Groups:
    """バッチ・HWE グループ・性別を解釈し、+fill-tags -S 用のファイルを書く

    属性: batches, hwe_groups, size{group: n}, sex{sample: 'M'/'F'} (or None), path (fill-tags 用ファイル)
          x_hwe_groups / x_batches_f / x_batches_m : chrX の非 PAR で使うグループ (人数 0 のものは除く)
    """

    def __init__(self, a, samples, tmpdir, log, min_batches=1, max_batches=None):
        samples_set = set(samples)
        member = collections.defaultdict(list)       # sample -> [group]
        if a.groups:
            for line in open(a.groups):
                f = line.split()
                if len(f) < 2 or f[0].startswith('#'):
                    continue
                if f[0] in samples_set:
                    member[f[0]] += f[1].split(',')
        if a.batches:
            self.batches = a.batches.split(',')
            if not a.groups:
                sys.exit('--batches を使うときは -S (groups.txt) が必要です')
        else:
            self.batches = [SINGLE_BATCH]
            for s in samples:
                member[s].append(SINGLE_BATCH)
        self.hwe_groups = a.hwe_groups.split(',') if a.hwe_groups else list(self.batches)
        if a.hwe_groups and not a.groups:
            sys.exit('--hwe-groups を使うときは -S (groups.txt) が必要です')
        if len(self.batches) < min_batches:
            sys.exit(f'--batches には {min_batches} つ以上のグループを指定してください')
        if max_batches and len(self.batches) > max_batches:
            sys.exit(f'--batches は {max_batches} つまでです')
        for g in set(self.batches + self.hwe_groups):
            if not TAG_RE.fullmatch(g) or g.endswith((FEMALE, MALE)):
                sys.exit(f'グループ名は英数字と _ のみ (末尾 {FEMALE}/{MALE} は予約) です: {g}')
        for s in samples:
            if len(s) <= 2 and member.get(s):
                sys.exit(f'bcftools 1.19 の +fill-tags は 2 文字以下のサンプル名を扱えません: {s}')
        for s, gs in member.items():
            if sum(g in self.batches for g in gs) > 1:
                sys.exit(f'サンプル {s} が複数のバッチに属しています')
        self.sex = read_sex(a.sex, samples_set) if a.sex else None

        size = collections.Counter()
        self.path = os.path.join(tmpdir, 'fill_tags_groups.txt')
        with open(self.path, 'w') as fo:
            for s in samples:
                gs = list(dict.fromkeys(member.get(s, [])))
                if self.sex and s in self.sex:
                    sx = FEMALE if self.sex[s] == 'F' else MALE
                    gs += [g + sx for g in gs if g in self.batches]
                    if sx == FEMALE:
                        gs += [g + FEMALE for g in gs if g in self.hwe_groups and g not in self.batches]
                gs = list(dict.fromkeys(gs))
                if gs:
                    fo.write(f'{s}\t{",".join(gs)}\n')
                    size.update(gs)
        self.size = size
        for g in set(self.batches + self.hwe_groups):
            if size[g] == 0:
                sys.exit(f'グループ {g} に該当するサンプルが VCF にありません')
        in_batch = sum(1 for s in samples if any(g in self.batches for g in member.get(s, [])))
        msg = ', '.join(f'{g}={size[g]}' for g in dict.fromkeys(self.batches + self.hwe_groups))
        log(f'groups: {msg} | VCF samples={len(samples)}, in a batch={in_batch}')
        if self.sex is not None:
            nf = sum(v == 'F' for v in self.sex.values())
            log(f'sex: female={nf}, male={len(self.sex) - nf}, unknown={len(samples) - len(self.sex)}')
            self.x_hwe_groups = [g + FEMALE for g in self.hwe_groups if size[g + FEMALE] > 0]
            self.x_batches_f = [g + FEMALE for g in self.batches if size[g + FEMALE] > 0]
            self.x_batches_m = [g + MALE for g in self.batches if size[g + MALE] > 0]
            for g in self.hwe_groups:
                if size[g + FEMALE] == 0:
                    log(f'WARNING: {g} に女性がいないため、chrX (PAR 以外) の HWE 判定から外します')
        else:
            self.x_hwe_groups = self.x_batches_f = self.x_batches_m = []


# ---------------------------------------------------------------------------
class Regions:
    """対象染色体・chrX の名前と PAR"""

    def __init__(self, a, hdr, log):
        self.lengths = dict((k, int(v)) for k, v in re.findall(r'##contig=<ID=([^,>]+),length=(\d+)', hdr))
        names = re.findall(r'##contig=<ID=([^,>]+)', hdr)
        pre = 'chr' if any(c.startswith('chr') for c in names) else ''
        self.x = f'{pre}X'
        if a.chroms:
            chroms = a.chroms.split(',')
        else:
            chroms = [f'{pre}{i}' for i in range(1, 23)] + ([self.x] if a.sex else [])
        bad = [c for c in chroms if re.fullmatch(r'(chr)?(Y|M|MT)', c)]
        if bad:
            sys.exit(f'{",".join(bad)} は対象外です (倍数性の扱いが異なるため)')
        if names:
            chroms = [c for c in chroms if c in names]
        if not chroms:
            sys.exit('対象染色体がヘッダーの ##contig に見つかりません (--chroms を確認)')
        self.chroms = chroms
        self.with_x = self.x in chroms
        self.par = []
        if self.with_x:
            if not a.sex:
                sys.exit(f'{self.x} を対象にするには --sex (性別ファイル) が必要です')
            if a.par:
                self.par = [tuple(int(v) for v in r.split('-')) for r in a.par.split(',')]
                build = 'manual'
            elif self.lengths.get(self.x) in PAR_BY_LENGTH:
                build, self.par = PAR_BY_LENGTH[self.lengths[self.x]]
            else:
                sys.exit(f'{self.x} の長さから参照配列を判定できません。--par で PAR を指定してください')
            log(f'{self.x}: PAR ({build}) = ' + ', '.join(f'{s}-{e}' for s, e in self.par))

    def expr_x_nonpar(self):
        """bcftools 式: chrX の PAR 以外の record (1.19 は否定 ! が使えないため PAR の補集合を書く)"""
        par = sorted(self.par)
        if not par:
            return f'(CHROM="{self.x}")'
        parts = [f'POS<{par[0][0]}']
        parts += [f'(POS>{e1} && POS<{s2})' for (_, e1), (s2, _) in zip(par, par[1:])]
        parts.append(f'POS>{par[-1][1]}')
        return f'(CHROM="{self.x}" && ({" || ".join(parts)}))'

    def expr_not_x_nonpar(self):
        """bcftools 式: chrX の PAR 以外 *ではない* record (常染色体 + PAR)"""
        inpar = ''.join(f' || (POS>={s} && POS<={e})' for s, e in sorted(self.par))
        return f'(CHROM!="{self.x}"{inpar})'

    def is_x_nonpar(self, chrom, pos):
        """numpy/pandas 用: chrom, pos は配列"""
        import numpy as np
        inpar = np.zeros(len(pos), bool)
        for s, e in self.par:
            inpar |= (pos >= s) & (pos <= e)
        return (chrom == self.x) & ~inpar

    def chunks(self, a, inputs):
        """(region or None, input) のリスト"""
        if a.input_list or a.no_chunk:
            return [(None, f) for f in inputs]
        if not has_index(a.input):
            sys.exit(f'index not found: {a.input}.csi / .tbi (bcftools index で作成するか --no-chunk を指定)')
        step = a.chunk_mb * 1_000_000
        return [(f'{c}:{s + 1}-{min(s + step, self.lengths[c])}', a.input)
                for c in self.chroms for s in range(0, self.lengths.get(c, 0), step)]

    def view_cmd(self, region, inp, extra=''):
        """chunk の切り出し (POS が領域内の record のみ → 境界で重複しない)"""
        if region:
            return f'bcftools view --no-version -r {q(region)} --regions-overlap 0 {extra} -Ou {q(inp)}'
        return f'bcftools view --no-version -t {q(",".join(self.chroms))} {extra} -Ou {q(inp)}'


# ---------------------------------------------------------------------------
class Mask:
    """遺伝型マスク: gt_filter (式), het_binom (p), desc (ログ用)"""

    def __init__(self, a, hdr, first_input, regions, log):
        has = lambda t: f'##FORMAT=<ID={t},' in hdr
        self.gt_filter = self.het_binom = None
        if a.no_gt_filter:
            if has('GQ') or has('DP'):
                frac = self.unmasked_fraction(first_input, regions)
                if frac > 0.001:
                    log(f'WARNING: --no-gt-filter ですが、先頭 200 サイトの呼ばれた遺伝型の {frac:.1%} が '
                        'GQ<20・DP<10・DP>65 のいずれかです。入力がマスクされていない可能性があります '
                        '(欠測率の判定が効かなくなります)')
            self.desc = 'none (--no-gt-filter)'
        else:
            if a.gt_filter == DEF_GT_FILTER:
                terms = [e for t, e in [('GQ', 'FMT/GQ<20'), ('DP', 'FMT/DP<10'), ('DP', 'FMT/DP>65')] if has(t)]
                if len(terms) < 3:
                    log('WARNING: FORMAT/GQ または FORMAT/DP が無いため、既定の遺伝型マスクの一部を省略します')
                self.gt_filter = ' | '.join(terms) or None
            elif a.gt_filter:
                self.gt_filter = a.gt_filter
            if a.het_binom > 0:
                if has('AD'):
                    self.het_binom = a.het_binom
                else:
                    log('WARNING: FORMAT/AD が無いため --het-binom を省略します')
            self.desc = f'{self.gt_filter or "-"} ; het AD binomial p<{self.het_binom or "-"}'
        log('genotype mask:', self.desc)

    @property
    def active(self):
        return bool(self.gt_filter or self.het_binom)

    def cmds(self):
        c = []
        if self.gt_filter:
            c.append(f'bcftools +setGT --no-version -Ou -- -t q -n . -i {q(self.gt_filter)}')
        if self.het_binom:
            c.append(f'bcftools +setGT --no-version -Ou -- -t {q(f"b:AD<{self.het_binom}")} -n .')
        return c

    @staticmethod
    def unmasked_fraction(path, regions, n_sites=200):
        """先頭 n_sites サイトで、欠測でないのに既定のマスク条件 (GQ<20 / DP<10 / DP>65) に当たる遺伝型の割合"""
        cmd = (f'bcftools view -Ou -t {q(",".join(regions.chroms))} {q(path)} '
               "| bcftools query -f '[%GT:%GQ:%DP\\t]\\n'")
        p = subprocess.Popen(['bash', '-c', cmd], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        n = low = 0
        num = lambda x: None if x in ('', '.') else float(x)
        for k, line in enumerate(p.stdout):
            if k >= n_sites:
                break
            for g in line.rstrip('\t\n').split('\t'):
                gt, gq, dp = (g.split(':') + ['', ''])[:3]
                if gt.startswith('.'):
                    continue
                n += 1
                gq, dp = num(gq), num(dp)
                low += (gq is not None and gq < 20) or (dp is not None and (dp < 10 or dp > 65))
        p.kill()
        p.wait()
        return low / n if n else 0.0
