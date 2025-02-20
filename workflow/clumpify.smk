# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config["r2_fastq"]
clumped_output = config["clumped_dir"]

# Define temporary and final output file names
tmp_r1 = f"{sample}_R1.tmp.fastq.gz"
tmp_r2 = f"{sample}_R2.tmp.fastq.gz"

final_r1 = f"{clumped_output}/{sample}_R1.clumped.fastq.gz"
final_r2 = f"{clumped_output}/{sample}_R2.clumped.fastq.gz"

rule all:
    input:
        final_r1,
        final_r2

rule clumpify:
    input:
        r1 = r1_fastq,
        r2 = r2_fastq
    output:
        final_r1 = final_r1,
        final_r2 = final_r2
    params:
        sample = sample,
        tmp_r1 = tmp_r1,
        tmp_r2 = tmp_r2,
        clumped_output = clumped_output
    shell:
        """
        # Copy input FASTQ files into working directory
        cp {input.r1} .
        cp {input.r2} .

        # Obtain base names for clumpify
        r1_basename=$(basename {input.r1})
        r2_basename=$(basename {input.r2})

        # Run clumpify
        clumpify.sh in=$r1_basename in2=$r2_basename \
            out={params.tmp_r1} out2={params.tmp_r2} \
            zl=9 pigz fastawrap=100000 threads=16 -Xmx100g

        # Create output directory and move clumped files
        mkdir -p {params.clumped_output}
        mv {params.tmp_r1} {output.final_r1}
        mv {params.tmp_r2} {output.final_r2}
        """
