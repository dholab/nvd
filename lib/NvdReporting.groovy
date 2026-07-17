import nextflow.util.ArrayTuple

class NvdReporting {
    static List processReadyFastqcTuples(meta, readsValue) {
        return expandFastqcUnits(meta, readsValue).collect { planned ->
            new ArrayTuple([planned.unit, planned.read])
        }
    }

    static List<Map> expandFastqcUnits(meta, readsValue) {
        def reads = readsValue
        def readMode = meta.read_mode.toString()
        if (readMode == 'single') {
            return reads.withIndex().collect { read, index ->
                unit(meta, read, 'single', index + 1)
            }
        }
        if (readMode == 'paired') {
            def r1Count = meta.r1_count as Integer
            if (r1Count <= 0 || reads.size() != r1Count * 2) {
                throw new IllegalArgumentException("Paired FastQC planning for ${meta.id} expected ${r1Count * 2} files from r1_count=${r1Count}, received ${reads.size()}")
            }
            def units = []
            reads[0..<r1Count].eachWithIndex { read, index ->
                units << unit(meta, read, 'R1', index + 1)
            }
            reads[r1Count..<reads.size()].eachWithIndex { read, index ->
                units << unit(meta, read, 'R2', index + 1)
            }
            return units
        }
        throw new IllegalArgumentException("Unsupported read_mode for FastQC planning on ${meta.id}: ${readMode}")
    }

    private static Map unit(meta, read, String readEnd, Integer ordinal) {
        def sampleId = meta.id.toString()
        def alias = aliasFor(sampleId, readEnd, ordinal)
        def zipAlias = "${alias}_fastqc.zip".toString()
        def htmlAlias = "${alias}_fastqc.html".toString()
        def packageName = "${alias}.fastqc".toString()
        [
            unit: [
                sample_id: sampleId,
                platform: meta.platform.toString(),
                source: meta.source.toString(),
                read_end: readEnd,
                input_ordinal: ordinal,
                alias: alias,
                zip_alias: zipAlias,
                html_alias: htmlAlias,
                staged_filename: "${alias}${fastqExtension(read)}".toString(),
                package_name: packageName,
            ],
            read: read,
        ]
    }

    private static String aliasFor(String sampleId, String readEnd, Integer ordinal) {
        return "${sampleId}.raw.${readEnd}.${String.format('%03d', ordinal)}".toString()
    }

    private static String fastqExtension(read) {
        def name = read.getName()
        if (name.endsWith('.fastq.gz')) {
            return '.fastq.gz'
        }
        if (name.endsWith('.fq.gz')) {
            return '.fq.gz'
        }
        if (name.endsWith('.fastq')) {
            return '.fastq'
        }
        if (name.endsWith('.fq')) {
            return '.fq'
        }
        throw new IllegalArgumentException("Unsupported FASTQ suffix for raw FastQC planning: ${name}")
    }
}
