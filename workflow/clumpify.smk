import os

# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config["r2_fastq"]
clumped_output = config['global']["clumped_out_dir"]

# Define final output paths
final_r1 = f"{clumped_output}/{sample}/{sample}_R1.clumped.fastq.gz"
final_r2 = f"{clumped_output}/{sample}/{sample}_R2.clumped.fastq.gz"

rule clumpify:
    input:
        r1 = r1_fastq,
        r2 = r2_fastq
    output:
        # Mark intermediate files as temporary with temp()
        tmp_r1 = temp(f"{sample}_R1.tmp.fastq.gz"),
        tmp_r2 = temp(f"{sample}_R2.tmp.fastq.gz"),
        final_r1 = final_r1,
        final_r2 = final_r2
    params:
        sample = sample,
        output_dir = lambda wildcards, output: os.path.dirname(output.final_r1)
    threads: 16
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p {params.output_dir}
        
        # Run clumpify using the original input files, writing to temporary outputs
        clumpify.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.tmp_r1} \
            out2={output.tmp_r2} \
            zl=9 pigz \
            fastawrap=100000 \
            threads={threads} \
            -Xmx100g
        
        # Copy the results from the temporary files to the final destination
        cp {output.tmp_r1} {output.final_r1}
        cp {output.tmp_r2} {output.final_r2}
        """