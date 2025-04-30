process MERGE_PAIRS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads1), path(reads2)

	output:
	tuple val(sample_id), val("illumina"), path("${sample_id}.merged.fastq.gz")

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

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads1), path(reads2)

	output:
	tuple val(sample_id), val("illumina"), path("${sample_id}.merged.fastq.gz")

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

process MASK_LOW_COMPLEXITY {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4
	memory 8.GB

	input:
	tuple val(sample_id), path(contigs)

	output:
	tuple val(sample_id), path("${sample_id}.masked.fasta")

	script:
	"""
	bbmask.sh -Xmx8g \
    in=${contigs} \
    out=${sample_id}.masked.fasta \
    entropy=${params.entropy} \
	-eoom
	"""
}

process FILTER_SHORT_CONTIGS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4
	memory 8.GB

	input:
	tuple val(sample_id), path(masked_contigs)

	output:
	tuple val(sample_id), path("${sample_id}.short_filtered.fasta")

	script:
	"""
	reformat.sh -Xmx8g \
    in=${masked_contigs} \
    out=${sample_id}.short_filtered.fasta \
    qtrim=${params.qtrim} \
    minconsecutivebases=${params.min_consecutive_bases} 
	"""
}

