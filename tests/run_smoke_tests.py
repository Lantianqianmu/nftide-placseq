#!/usr/bin/env python3
"""Exercise the real workflow with small fixtures.

cutadapt, MACS2, and bedGraphToBigWig are real executables. HiC-Pro and
FitHiChIP are explicit test doubles, so this checks routing, input selection,
publishing, and process handoffs, not alignment or FitHiChIP statistics.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile


MOCK_HICPRO = '''#!/usr/bin/env python3
import argparse,gzip,json,os,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('-i');p.add_argument('-o');p.add_argument('-c');a=p.parse_args()
sample=next(Path(a.i).iterdir()).name
config=dict(line.split(' = ',1) for line in Path(a.c).read_text().splitlines() if ' = ' in line)
assert config['RM_DUP']=='1' and config['GET_ALL_INTERACTION_CLASSES']=='1'
out=Path(a.o)
shutil.copytree(Path(os.environ['PLAC_TEST_FIXTURES'])/'hicpro'/sample,out)
counts={}
for f in (Path(a.i)/sample).glob('*.gz'):
    with gzip.open(f,'rt') as h:counts[f.name]=sum(1 for _ in h)//4
(out/'trimmed_read_counts.json').write_text(json.dumps(counts))
'''

MOCK_FITHICHIP = '''#!/usr/bin/env bash
set -euo pipefail
python3 - "$@" <<'PYTHON'
import sys,json,hashlib
from pathlib import Path
config=Path(sys.argv[sys.argv.index('-C')+1])
settings=dict(line.split('=',1) for line in config.read_text().splitlines() if '=' in line)
for key in ['ValidPairs','Interval','Matrix','ChrSizeFile','PeakFile']:
    p=Path(settings[key]);assert p.is_file() and p.stat().st_size>0,(key,p)
assert settings['BINSIZE']=='10000'
assert settings['UseP2PBackgrnd'] in ['0','1']
out=Path(settings['OutDir']);out.mkdir(exist_ok=True)
prefix=settings['PREFIX']
(out/('Summary_results_'+prefix+'.html')).write_text('<html>fixture</html>')
(out/(prefix+'.interactions_FitHiC.bed')).write_text('chr1\\ts1\\te1\\tchr2\\ts2\\te2\\tcc\\tq\\nchr1\\t10000\\t20000\\tchr1\\t50000\\t60000\\t2\\t0.5\\n')
(out/'input_check.json').write_text(json.dumps({'peak_sha256':hashlib.sha256(Path(settings['PeakFile']).read_bytes()).hexdigest(),'background':settings['UseP2PBackgrnd']}))
print('FINISHED calling significant interactions')
PYTHON
'''


def executable(path, content):
    path.write_text(content)
    path.chmod(0o755)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_fixture(root):
    root.mkdir(parents=True, exist_ok=True)
    sizes = root / 'chrom.sizes'
    sizes.write_text('chr1\t1000000\nchr2\t100000\n')
    (root / 'fragments.bed').write_text('chr1\t0\t1000000\tf1\nchr2\t0\t100000\tf2\n')
    index = root / 'index'
    for part in ['1', '2', '3', '4', 'rev.1', 'rev.2']:
        Path(str(index) + '.' + part + '.bt2').write_text('fixture\n')
    (root / 'provided').mkdir()
    provided = root / 'provided' / "owner's narrow peaks.bed.gz"
    provided.write_bytes(gzip.compress(b'chr1\t700000\t701000\n'))
    broad = root / 'provided' / 'external.broadPeak'
    broad.write_text('chr2\t5000\t15000\tref\t100\t.\t10\t20\t30\n')
    samples = [
        ('control-27ac-1', '', 2),
        ('control-27ac-2', str(provided), 1),
        ('control-27ac-3', 'shared/common narrow.bed.gz', 1),
        ('control-27ac-4', './shared/common narrow.bed.gz', 1),
        ('control-27me3-1', '', 1),
        ('control-27me3-2', str(broad), 1),
        ('control-27me3-3', 'shared/common.broadPeak', 1),
        ('control-27me3-4', './shared/common.broadPeak', 1),
        ('control-9ac-1', 'shared/solo.bed', 1),
    ]
    rows = []
    for number, (sample, peaks, lanes) in enumerate(samples):
        for lane in range(lanes):
            directory = root / 'reads' / sample / ('lane %d' % (lane + 1))
            directory.mkdir(parents=True)
            for mate in [1, 2]:
                # Repeated basenames across lanes exercise collision-safe staging.
                path = directory / ('reads_R%d.fq.gz' % mate)
                text = '@read%d/%d\n%s\n+\n%s\n' % (lane, mate, 'ACGT' * 20, 'I' * 80)
                path.write_bytes(gzip.compress(text.encode()))
            rows.append([sample, str(directory / 'reads_R1.fq.gz'), str(directory / 'reads_R2.fq.gz'), peaks])
        hic = root / 'hicpro' / sample / 'hic_results'
        pairs = hic / 'data' / sample
        pairs.mkdir(parents=True)
        offset = 10000 + 20000 * number
        enriched = ['e%d\tchr1\t%d\t+\tchr1\t%d\t-\n' % (i, offset + i + 1, offset + i + 601) for i in range(400)]
        (pairs / (sample + '.DEPairs')).write_text(''.join(enriched + [enriched[0]]))
        (pairs / (sample + '.SCPairs')).write_text(enriched[0])
        (pairs / (sample + '.REPairs')).write_text('re\tchr1\t400001\t+\tchr1\t400501\t-\n')
        (pairs / (sample + '.allValidPairs')).write_text(
            enriched[0] +
            'short\tchr1\t300001\t+\tchr1\t301000\t-\n' +
            'boundary\tchr1\t400001\t+\tchr1\t401001\t-\n' +
            'long\tchr1\t500001\t+\tchr1\t550001\t-\n' +
            'trans\tchr1\t600001\t+\tchr2\t60001\t-\n')
        matrix = hic / 'matrix' / sample / 'raw' / '10000'
        matrix.mkdir(parents=True)
        (matrix / (sample + '_10000.matrix')).write_text('1\t5\t2\n')
        (matrix / (sample + '_10000_abs.bed')).write_text('chr1\t0\t10000\t1\nchr1\t40000\t50000\t5\n')
    sheet = root / 'samplesheet.csv'
    with sheet.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample', 'fastq_1', 'fastq_2', 'peak_file'])
        writer.writerows(rows)
    tools = root / 'tools'
    tools.mkdir()
    executable(tools / 'HiC-Pro', MOCK_HICPRO)
    fithichip = root / 'fithichip'
    fithichip.mkdir()
    executable(fithichip / 'FitHiChIP_HiCPro.sh', MOCK_FITHICHIP)
    return samples, sheet, provided, broad


def groovy_string(value):
    return "'" + str(value).replace('\\', '\\\\').replace("'", "\\'") + "'"


def validate_outputs(root, samples, provided_hashes):
    out = root / 'output'
    for sample, peak_path, lanes in samples:
        directory = out / sample
        signal = directory / 'PLAC_signal' / 'signal'
        anchors = directory / 'FitHiChIP_peaks' / 'peaks'
        info = json.loads((anchors / 'source.json').read_text())
        expected = 'sample' if not peak_path else ('provided' if 'provided/' in peak_path else 'shared')
        assert info['source'] == expected, (sample, info)
        assert (anchors / 'peaks.bed').stat().st_size > 0
        stats = dict(line.split('\t') for line in (signal / 'read_selection.tsv').read_text().splitlines()[1:])
        assert int(stats['selected_pairs_after_dedup']) == 402, (sample, stats)
        assert int(stats['retained_read_ends']) == 804
        bed = (signal / 'MACS2_input.bed').read_text().splitlines()
        assert len(bed) == 804
        assert not any(line.startswith('chr2\t') for line in bed)
        assert not any(line.startswith('chr1\t500000\t') for line in bed)
        assert stats['macs2_peak_type'] == ('broad' if '27me3' in sample else 'narrow')
        bigwig = directory / 'bw' / (sample + '.plac_1d.cpm.bw')
        assert bigwig.stat().st_size > 100
        assert struct.unpack('<I', bigwig.read_bytes()[:4])[0] == 0x888FFC26
        for line in (signal / 'plac_1d.cpm.bedGraph').read_text().splitlines():
            chrom, start, end, value = line.split()
            assert 0 <= int(start) < int(end) <= (1000000 if chrom == 'chr1' else 100000)
            assert float(value) >= 0
        reads = json.loads((directory / 'hicpro' / 'trimmed_read_counts.json').read_text())
        assert list(reads.values()) == [lanes, lanes]
        for mode in ['L', 'S']:
            checked = json.loads((directory / ('FitHiChIP_' + mode) / 'input_check.json').read_text())
            assert checked['peak_sha256'] == digest(anchors / 'peaks.bed')
        if expected == 'provided':
            original = Path(peak_path)
            text = gzip.decompress(original.read_bytes()).decode() if original.suffix == '.gz' else original.read_text()
            assert (anchors / 'peaks.bed').read_text() == text
        print('PASS:', sample, expected, 'peaks, 1D bigWig, L/S handoffs', flush=True)
    groups = list((out / 'shared_peaks').glob('*/peaks/group.json'))
    assert len(groups) == 3, groups
    for group_file in groups:
        info = json.loads(group_file.read_text())
        pooled = group_file.parent / 'MACS2_input.bed'
        expected_count = sum(len((out / sample / 'PLAC_signal/signal/MACS2_input.bed').read_text().splitlines()) for sample in info['source_samples'])
        assert len(pooled.read_text().splitlines()) == expected_count
        emitted = Path(info['path'])
        assert emitted.exists()
        text = gzip.decompress(emitted.read_bytes()).decode() if emitted.suffix == '.gz' else emitted.read_text()
        assert text == (group_file.parent / 'peaks.bed').read_text()
        for sample in info['source_samples']:
            assert (out / sample / 'FitHiChIP_peaks/peaks/peaks.bed').read_text() == text
        print('PASS: shared peak group', info['source_samples'], 'pooled reads and exported file', flush=True)
    for filename, original_hash in provided_hashes.items():
        assert digest(Path(filename)) == original_hash


def validate_resume(trace, expected_shared):
    with trace.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    upstream = [row for row in rows if row['name'].startswith(
        ('MERGE_FQ', 'CUTADAPT', 'HICPRO', 'PLAC_SIGNAL'))]
    assert len(upstream) == 36 and all(row['status'] == 'CACHED' for row in upstream)
    shared = [row for row in rows if row['name'].startswith('SHARED_PEAKS')]
    assert len(shared) == expected_shared and all(row['status'] == 'CACHED' for row in shared)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path)
    parser.add_argument('--bedgraph-to-bigwig', required=True, type=Path)
    args = parser.parse_args()
    root = (args.work_dir or Path(tempfile.mkdtemp(prefix='placseq-smoke-'))).resolve()
    pipeline = Path(__file__).resolve().parents[1]
    samples, sheet, provided, broad = make_fixture(root)
    original = {str(p): digest(p) for p in [provided, broad]}
    config = root / 'test.config'
    overrides = dict(input_csv=sheet, bowtie2_index=root / 'index', chrom_sizes=root / 'chrom.sizes',
                     genome_fragment=root / 'fragments.bed', hicpro=root / 'tools/HiC-Pro',
                     fithichip_dir=root / 'fithichip', bedgraph_to_bigwig=args.bedgraph_to_bigwig.resolve())
    config.write_text('\n'.join('params.%s = %s' % (key, groovy_string(value)) for key, value in overrides.items()) +
                      "\nparams.hicpro_cpus = 1\nparams.fithichip_cpus = 1\nparams.fithichip_peak_genome = '1000000'\nprocess.cpus = 1\nprocess.maxForks = 2\nworkflow.output.mode = 'copy'\n" +
                      'outputDir = ' + groovy_string(root / 'output') + '\n')
    env = os.environ.copy()
    env['PLAC_TEST_FIXTURES'] = str(root)
    env['NXF_HOME'] = str(root / 'nxf-home')
    env['NXF_OFFLINE'] = 'true'
    env['NXF_DISABLE_CHECK_LATEST'] = 'true'

    # Run the actual planning helpers in Nextflow, without launching analysis tasks.
    helpers = (pipeline / 'main.nf').read_text().split('\nprocess MERGE_FQ {', 1)[0]
    helper_file = root / 'planning-helpers.groovy'
    helper_file.write_text('\n'.join(line for line in helpers.splitlines()
                                    if not line.startswith(('#!', 'nextflow.enable.dsl'))))
    planning = root / 'planning-tests.nf'
    cases = (pipeline / 'tests' / 'sample_planning.groovy').read_text()
    # Evaluate the helpers as plain Groovy. Nextflow's function wrapper aborts
    # the session on intentional validation errors before a caller can catch them.
    planning.write_text('nextflow.enable.dsl = 2\nworkflow {\n' +
                        'def root = new File(' + groovy_string(root) + ')\n' +
                        'def helper = new GroovyShell().parse(new File(root, "planning-helpers.groovy"))\n' +
                        cases + '\n}\n')
    subprocess.run(['nextflow', '-log', str(root / 'planning.log'), 'run', str(planning),
                    '-ansi-log', 'false'], cwd=root, env=env, check=True)

    command = ['nextflow', '-log', str(root / 'nextflow.log'), 'run', str(pipeline / 'main.nf'),
               '-c', str(config), '-w', str(root / 'work'), '-ansi-log', 'false',
               '-with-trace', str(root / 'trace.tsv')]
    subprocess.run(command, cwd=root, env=env, check=True)
    validate_outputs(root, samples, original)

    # Generated shared paths become supplied inputs on the next launch. Alignment
    # and signal must stay cached, and the shared MACS2 tasks must not run again.
    resume_command = command[:-2] + ['-with-trace', str(root / 'resume-trace.tsv'), '-resume']
    subprocess.run(resume_command, cwd=root, env=env, check=True)
    validate_resume(root / 'resume-trace.tsv', expected_shared=0)
    print('PASS: generated peak paths reused with all 36 upstream tasks cached', flush=True)

    # A cached shared call must restore its requested output when it is removed.
    shared_peak = root / 'shared' / 'common narrow.bed.gz'
    original_content = shared_peak.read_bytes()
    shared_peak.unlink()
    republish_command = command[:-2] + ['-with-trace', str(root / 'republish-trace.tsv'), '-resume']
    subprocess.run(republish_command, cwd=root, env=env, check=True)
    validate_resume(root / 'republish-trace.tsv', expected_shared=1)
    assert shared_peak.read_bytes() == original_content
    print('PASS: cached shared peak call restores the requested output file', flush=True)
    print('PASS: smoke workflow. Artifacts:', root, flush=True)


if __name__ == '__main__':
    main()
