# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
# From command line
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config["r2_fastq"]

# From parameters file - updated path name
clumped_output = config['global']["clumped_out_dir"]

# Define temporary and final output file names
tmp_r1 = f"{sample}_R1.tmp.fastq.gz"  # Removed tmp/ prefix
tmp_r2 = f"{sample}_R2.tmp.fastq.gz"  # Removed tmp/ prefix

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
        tmp_r1 = tmp_r1,
        tmp_r2 = tmp_r2,
        output_dir = lambda wildcards, output: os.path.dirname(output.final_r1)
    shell:
        """
        # Create and move to temporary directory
        mkdir -p tmp
        cd tmp

        # Copy input FASTQ files to current directory
        cp {input.r1} .
        cp {input.r2} .

        # Get base names for input files
        r1_basename=$(basename {input.r1})
        r2_basename=$(basename {input.r2})

        # Run clumpify
        clumpify.sh \
            in=$r1_basename \
            in2=$r2_basename \
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

        # Return to original directory and cleanup
        cd ..
        rm -rf tmp/
        """