process SCRUB_HOST_READS {

	/* Remove host (human) reads using hostile */

	tag "${sample_id}"
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), path(reads)

	output:
	tuple val(sample_id), val(platform), path("${sample_id}.scrubbed.fastq.gz")

	script:
	def aligner = platform == "illumina" ? "bowtie2" : "minimap2"
	def index_arg = params.hostile_index
		? "--index ${params.hostile_index}"
		: "--index ${params.hostile_index_name}"
	def airplane_flag = params.hostile_index ? "--airplane" : ""
	"""
	hostile clean \\
		--fastq1 ${reads} \\
		--aligner ${aligner} \\
		${index_arg} \\
		${airplane_flag} \\
		--threads ${task.cpus} \\
		--force

	mv *.clean.fastq.gz ${sample_id}.scrubbed.fastq.gz
	"""
}
