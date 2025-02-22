import os

# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config["r2_fastq"]
clumped_output = config['global']["clumped_out_dir"]

# Define output paths
final_r1 = f"{clumped_output}/{sample}/{sample}_R1.clumped.fastq.gz"
final_r2 = f"{clumped_output}/{sample}/{sample}_R2.clumped.fastq.gz"

rule clumpify:
    input:
        r1 = r1_fastq,
        r2 = r2_fastq
    output:
        final_r1 = final_r1,
        final_r2 = final_r2
    params:
        sample = sample,
        tmp_dir = "tmp",
        tmp_r1 = "tmp/R1.tmp.fastq.gz",
        tmp_r2 = "tmp/R2.tmp.fastq.gz",
        output_dir = lambda wildcards, output: os.path.dirname(output.final_r1)
    shell:
        """
        # Create temporary directory
        mkdir -p {params.tmp_dir}

        # Copy input files to temp directory
        cp {input.r1} {params.tmp_dir}/input_R1.fastq.gz
        cp {input.r2} {params.tmp_dir}/input_R2.fastq.gz

        # Run clumpify with explicit paths
        clumpify.sh \
            in={params.tmp_dir}/input_R1.fastq.gz \
            in2={params.tmp_dir}/input_R2.fastq.gz \
            out={params.tmp_r1} \
            out2={params.tmp_r2} \
            zl=9 pigz \
            fastawrap=100000 \
            threads=16 \
            -Xmx100g

        # Create output directory and move files
        mkdir -p {params.output_dir}
        mv {params.tmp_r1} {output.final_r1}
        mv {params.tmp_r2} {output.final_r2}

        # Cleanup
        rm -rf {params.tmp_dir}
        """