process MERGE_PAIRS {

	/* Merge overlapping paired-end reads into single longer reads */

	tag "${sample_id}"
	label "medium"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), path(reads1), path(reads2)

	output:
	tuple val(sample_id), val(platform), val("merged"), path("${sample_id}.merged.fastq.gz")

	script:
	"""
	bbmerge.sh \
	in=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${sample_id}.merged.fastq.gz \
	outu=${sample_id}.unmerged.fastq.gz \
	threads=${task.cpus} \
	-eoom
	"""
}

process INTERLEAVE_PAIRS {

	/* Interleave paired-end reads into a single file */

	tag "${sample_id}"
	label "medium"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), path(reads1), path(reads2)

	output:
	tuple val(sample_id), val(platform), val("interleaved"), path("${sample_id}.interleaved.fastq.gz")

	when:
	params.tools

	script:
	"""
	reformat.sh \
	in=${reads1} \
	in2=${reads2} \
	out=${sample_id}.interleaved.fastq.gz \
	threads=${task.cpus} \
	-eoom
	"""
}

process DEDUP_WITH_CLUMPIFY {

    /* Deduplicate reads using clumpify */

	tag "${sample_id}"
	label "ludicrous"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(reads)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.dedup.fastq.gz")

	script:
	def int_flag = read_structure == "interleaved" ? "int=t" : ""
	"""
	clumpify.sh \\
	in=${reads} \\
	out="${sample_id}.dedup.fastq.gz" \\
	dedupe=2 \\
	reorder=p \\
	subs=2 \\
	${int_flag} \\
	threads=${task.cpus} -eoom
	"""

}

process TRIM_ADAPTERS {

	/* Trim Illumina adapters using bbduk */

	tag "${sample_id}"
	label "ludicrous"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(reads)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.trimmed.fastq.gz")

	script:
	"""
	bbduk.sh \\
		in=${reads} \\
		out=${sample_id}.trimmed.fastq.gz \\
		ref=adapters \\
		ktrim=r \\
		k=23 \\
		mink=11 \\
		hdist=1 \\
		tpe tbo \\
		threads=${task.cpus} \\
		-eoom
	"""
}

process FILTER_READS {

	/* Filter reads by quality and length (pair-aware for interleaved reads) */

	tag "${sample_id}"
	label "ludicrous"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(reads), val(min_quality)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.filtered.fastq.gz")

	script:
	def max_len_arg = params.max_read_length ? "maxlength=${params.max_read_length}" : ""
	def interleaved_arg = read_structure == "interleaved" ? "int=t" : ""
	"""
	reformat.sh \\
		in=${reads} \\
		out=${sample_id}.filtered.fastq.gz \\
		${interleaved_arg} \\
		minavgquality=${min_quality} \\
		minlength=${params.min_read_length} \\
		${max_len_arg} \\
		threads=${task.cpus} \\
		-eoom
	"""
}

process REPAIR_PAIRS {

	/* Repair interleaved paired-end reads, discarding orphans */

	tag "${sample_id}"
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(reads)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.repaired.fastq.gz")

	script:
	"""
	repair.sh \\
		in=${reads} \\
		out=${sample_id}.repaired.fastq.gz \\
		outs=${sample_id}.singletons.fastq.gz \\
		interleaved=t \\
		repair=t \\
		threads=${task.cpus} \\
		-eoom
	"""
}

process MASK_LOW_COMPLEXITY {

	/* Mask low-complexity regions in contigs */

	tag "${sample_id}"
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4
	memory 8.GB

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(contigs)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.masked.fasta")

	script:
	"""
	# Check if input has any sequences (BBMask crashes on empty input)
	if [ ! -s ${contigs} ] || ! grep -q "^>" ${contigs}; then
	    # Empty or no sequences - create empty output
	    touch ${sample_id}.masked.fasta
	else
	    bbmask.sh -Xmx8g \
	        in=${contigs} \
	        out=${sample_id}.masked.fasta \
	        entropy=${params.entropy} \
	        threads=${task.cpus} \
	        -eoom
	fi
	"""
}

process FILTER_SHORT_CONTIGS {

	/* Filter out short contigs */

	tag "${sample_id}"
	label "medium"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4
	memory 8.GB

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(masked_contigs)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.short_filtered.fasta")

	script:
	"""
	reformat.sh -Xmx8g \
    in=${masked_contigs} \
    out=${sample_id}.short_filtered.fasta \
    qtrim=${params.qtrim} \
    minconsecutivebases=${params.min_consecutive_bases} 
	"""
}

process CLUMP_READS {

    /* */

	tag "${sample_id}"
	label "ludicrous"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
    tuple val(sample_id), path("${sample_id}.clumped.fastq.gz")

    when:
	params.tools && (params.tools.contains("clump") || params.tools.contains("all"))

	script:
	"""
	clumpify.sh \
	in=${reads} out=${sample_id}.clumped.fastq.gz \
	reorder \
	threads=${task.cpus} -eoom
	"""

}

process SANITIZE_EXTRACTED_FASTA {

    /* Reformat GOTTCHA2 extracted FASTA for downstream consumption.
     * Taxonomic-level deduplication (strain > species) is handled by
     * each consumer: register_gottcha2_hits.py and labkey_upload_gottcha2_fasta.py. */

	tag "${sample_id}"
	label "medium"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(extracted_reads), path(full_tsv)

	output:
    tuple val(sample_id), path("${sample_id}.no_ambig.fasta"), path(full_tsv)

	script:
	"""
	reformat.sh \
	in=${extracted_reads} out=${sample_id}.no_ambig.fasta \
	nullifybrokenquality=t \
	threads=${task.cpus} -eoom
	"""
}

