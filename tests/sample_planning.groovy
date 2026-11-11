// Embedded by run_smoke_tests.py after extracting the actual main.nf helpers.
def types = ['27ac':'narrow', '27me3':'broad', 'pol2':'narrow']
def csv = new File(root, 'samplesheet.csv').toPath()
def fixture = new File(root, 'reads/control-27ac-1/lane 1')
def row = [sample:'test-27ac-1',
           fastq_1:new File(fixture,'reads_R1.fq.gz').toString(),
           fastq_2:new File(fixture,'reads_R2.fq.gz').toString(), peak_file:'']
def expectError = { String label, String message, Closure action ->
    try { action(); assert false: "Expected failure: ${label}" }
    catch (IllegalArgumentException error) { assert error.message.contains(message): error.message }
    println "PASS: rejected ${label}"
}
expectError('missing peak_file column','CSV requires') {
    helper.prepareSamples([row.findAll { k,v -> k!='peak_file' }], csv, types)
}
expectError('inconsistent lane peak values','All rows for sample') {
    helper.prepareSamples([new LinkedHashMap(row), row + [peak_file:'missing.bed']], csv, types)
}
expectError('mixed broad/narrow group','mixes broad and narrow') {
    helper.prepareSamples([row + [peak_file:'never_created/common.bed'],
                    row + [sample:'test-27me3-1',peak_file:'never_created/common.bed']], csv, types)
}
expectError('existing directory as peaks','not a readable, nonempty file') {
    helper.prepareSamples([row + [peak_file:fixture.toString()]], csv, types)
}
def empty = new File(root, 'empty.bed'); empty.text=''
expectError('empty supplied peaks','not a readable, nonempty file') {
    helper.prepareSamples([row + [peak_file:empty.toString()]], csv, types)
}
expectError('unknown mark','Cannot determine MACS2 peak type') {
    helper.prepareSamples([row + [sample:'unknown-1']], csv, types)
}
expectError('broad call summits','incompatible') { helper.macs2Settings('broad','-q 0.05 --call-summits',0.1) }
expectError('broad options in shared settings','Set fithichip_peak_types') { helper.macs2Settings('narrow','-q 0.05 --broad',0.1) }
expectError('invalid broad cutoff','between 0 and 1') { helper.macs2Settings('broad','-q 0.05',1) }
assert helper.resolvePeakType('test-H3K27ME3-2',types)=='broad'
assert helper.resolvePeakType('test-27ac-1',types + ['test-27ac-1':'broad'])=='broad'
def samples=helper.prepareSamples([row + [peak_file:'never_created/../common_new.bed'],
                           row + [sample:'test-27ac-2',peak_file:'./common_new.bed']],csv,types)
assert samples[0][3].path==samples[1][3].path
assert samples[0][3].mode=='shared'
println 'PASS: relative-path normalization, full-sample overrides and numeric replicate suffixes'
