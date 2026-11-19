#!/usr/bin/env python3
"""Inspect a VCF/BCF before QC and create a missing index when possible."""
import argparse
import collections
import itertools
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile

from gwasqc_common import has_index, make_log, run


def count_filters(path, limit=1_000_000):
    counts = collections.Counter()
    cmd = ['bcftools', 'query', '-f', '%FILTER\n', path]
    # Keep stderr off the pipe so warnings cannot deadlock a large query.
    with tempfile.TemporaryFile(mode='w+t') as err:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True) as proc:
            assert proc.stdout is not None
            try:
                counts.update(line.rstrip('\n') for line in itertools.islice(proc.stdout, limit))
                stopped = sum(counts.values()) == limit and proc.poll() is None
                if stopped:
                    # Close before waiting: bcftools may still be blocked writing to stdout.
                    proc.stdout.close()
                    proc.terminate()
                code = proc.wait()
            finally:
                proc.stdout.close()
                if proc.poll() is None:
                    proc.terminate()
                    proc.wait()
            err.seek(0)
            diagnostics = err.read()
            if code != 0 and not (stopped and code in (-signal.SIGTERM, -signal.SIGPIPE)):
                raise subprocess.CalledProcessError(code, cmd, stderr=diagnostics)
            if diagnostics:
                print(diagnostics, end='', file=sys.stderr)
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-i', '--input', required=True, help='VCF/BCF to inspect')
    a = ap.parse_args()
    log = make_log('check_input_vcf')
    if not shutil.which('bcftools'):
        sys.exit('bcftools not found in PATH')
    if not os.path.isfile(a.input):
        sys.exit(f'file not found: {a.input}')

    try:
        samples = run(['bcftools', 'query', '-l', a.input]).splitlines()
        header = run(['bcftools', 'view', '-h', a.input]).splitlines()
        print(f'Sample count: {len(samples)}')

        contigs = [line for line in header if line.startswith('##contig=')]
        print('\nContig headers (first 3):')
        print('\n'.join(contigs[:3]) or '(none)')
        if not contigs:
            log('WARNING: no contig headers; contig lengths are needed for chrX PAR inference')
        for line in contigs:
            if not re.search(r'[,<]length=[1-9][0-9]*[,>]', line):
                log(f'WARNING: contig length missing or invalid: {line}')

        formats = [line.split(',', 1)[0] for line in header if line.startswith('##FORMAT=')]
        print('\nFORMAT headers:')
        print('\n'.join(formats) or '(none)')
        format_ids = {line.removeprefix('##FORMAT=<ID=') for line in formats}
        missing = [name for name in ('GT', 'GQ', 'DP', 'AD') if name not in format_ids]
        if missing:
            log(f'WARNING: missing FORMAT fields for genotype masking: {", ".join(missing)}')

        counts = count_filters(a.input)
        print(f'\nFILTER counts (first 1,000,000 records; inspected {sum(counts.values())}):')
        for label, count in sorted(counts.items()):
            print(f'{count:7d} {label}')

        if has_index(a.input):
            print('\nIndex: already exists')
        else:
            with open(a.input, 'rb') as f:
                plain_vcf = f.read(2) == b'##'
            if plain_vcf:
                print('\nIndex: not created (uncompressed VCF)')
                log('WARNING: uncompressed VCF cannot be indexed; use bcftools view -Oz '
                    'to create a BGZF VCF, then bcftools index. Input was not modified.')
            else:
                subprocess.run(['bcftools', 'index', a.input], check=True)
                print('\nIndex: created')
    except subprocess.CalledProcessError as exc:
        if exc.stderr:
            print(exc.stderr, end='', file=sys.stderr)
        sys.exit(f'bcftools failed (exit {exc.returncode}): {exc.cmd}')


if __name__ == '__main__':
    main()
