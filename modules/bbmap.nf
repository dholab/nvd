process MERGE_PAIRS {

	/* Merge overlapping paired-end reads into single longer reads */

	tag "${sample_id}"
	label "medium"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), path(reads)

	output:
	tuple val(sample_id), val(platform), val("single"), val("overlap_merged_pair"), path("${sample_id}.overlap_merged_pair.fastq.gz"), emit: merged
	tuple val(sample_id), val(platform), val("single"), val("single_read"), path("${sample_id}.single_read.fastq.gz"), emit: unmerged

	script:
	"""
	bbmerge.sh \
	in=${reads} \
	out=${sample_id}.overlap_merged_pair.fastq.gz \
	outu=${sample_id}.single_read.fastq.gz \
	interleaved=t \
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
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path(reads)

	output:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path("${sample_id}.${query_class}.dedup.fastq.gz")

	script:
	def int_flag = read_structure == "interleaved" ? "int=t" : ""
	"""
	clumpify.sh \\
	in=${reads} \\
	out="${sample_id}.${query_class}.dedup.fastq.gz" \\
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
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path(reads)

	output:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path("${sample_id}.${query_class}.trimmed.fastq.gz")

	script:
	def pair_trim_args = read_structure == "interleaved" ? "tpe tbo" : ""
	"""
	bbduk.sh \\
		in=${reads} \\
		out=${sample_id}.${query_class}.trimmed.fastq.gz \\
		ref=adapters \\
		ktrim=r \\
		k=23 \\
		mink=11 \\
		hdist=1 \\
		${pair_trim_args} \\
		threads=${task.cpus} \\
		-eoom
	"""
}

process FILTER_READS {

	/* Independently filter reads by quality/length and normalized k-mer entropy */

	tag "${sample_id}"
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path(reads), val(min_quality)

	output:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path("${sample_id}.${query_class}.filtered.fastq.gz")

	script:
	def max_len_arg = params.max_read_length ? "maxlength=${params.max_read_length}" : ""
	def quality_args = params.filter_reads
		? "minavgquality=${min_quality} minlength=${params.min_read_length} ${max_len_arg}"
		: "minavgquality=0 minlength=0"
	def entropy_args = params.filter_low_complexity_reads
		? "entropy=${params.min_read_entropy} entropywindow=50 entropyk=5"
		: ""
	def pair_args = read_structure == "interleaved"
		? "int=t removeifeitherbad=t"
		: "int=f"
	"""
	bbduk.sh \\
		in=${reads} \\
		out=${sample_id}.${query_class}.filtered.fastq.gz \\
		${pair_args} \\
		${quality_args} \\
		${entropy_args} \\
		ordered=t \\
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
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path(reads)

	output:
	tuple val(sample_id), val(platform), val(read_structure), val(query_class), path("${sample_id}.${query_class}.repaired.fastq.gz")

	script:
	"""
	repair.sh \\
		in=${reads} \\
		out=${sample_id}.${query_class}.repaired.fastq.gz \\
		outs=${sample_id}.${query_class}.singletons.fastq.gz \\
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
	tuple val(sample_id), val(platform), val(read_structure), path(contigs), path(query_lookup)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.masked.fasta"), path(query_lookup)

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
	tuple val(sample_id), val(platform), val(read_structure), path(masked_contigs), path(query_lookup)

	output:
	tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.short_filtered.fasta"), path(query_lookup)

	script:
	"""
	reformat.sh -Xmx8g \
    in=${masked_contigs} \
    out=${sample_id}.short_filtered.fasta \
    qtrim=${params.qtrim} \
    minconsecutivebases=${params.min_consecutive_bases} 
	"""
}
