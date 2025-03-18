import os

# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config.get("r2_fastq", None)  # Optional for ONT
is_nanopore = config.get("is_nanopore", False)  # New parameter to identify ONT data
clumped_output = config['global']["clumped_out_dir"]

# Define final output paths based on data type
if is_nanopore:
    final_outputs = [f"{clumped_output}/{sample}.clumped.fastq.gz"]
else:
    final_outputs = [
        f"{clumped_output}/{sample}_R1.clumped.fastq.gz",
        f"{clumped_output}/{sample}_R2.clumped.fastq.gz"
    ]

rule all:
    input:
        final_outputs

rule clumpify:
    input:
        r1 = r1_fastq,
        r2 = r2_fastq if not is_nanopore else []
    output:
        # For ONT data
        tmp_ont = temp(f"{sample}.tmp.fastq.gz") if is_nanopore else [],
        final_ont = f"{clumped_output}/{sample}.clumped.fastq.gz" if is_nanopore else [],
        # For Illumina data
        tmp_r1 = temp(f"{sample}_R1.tmp.fastq.gz") if not is_nanopore else [],
        tmp_r2 = temp(f"{sample}_R2.tmp.fastq.gz") if not is_nanopore else [],
        final_r1 = f"{clumped_output}/{sample}_R1.clumped.fastq.gz" if not is_nanopore else [],
        final_r2 = f"{clumped_output}/{sample}_R2.clumped.fastq.gz" if not is_nanopore else []
    params:
        sample = sample,
        output_dir = clumped_output
    threads: 16
    run:
        # Create output directory if it doesn't exist
        shell(f"mkdir -p {params.output_dir}")
        
        if is_nanopore:
            # Run clumpify on ONT data (single-end)
            shell("""
            clumpify.sh \
                in={input.r1} \
                out={output.tmp_ont} \
                zl=9 pigz \
                fastawrap=100000 \
                threads={threads} \
                -Xmx100g
            
            # Copy the results from the temporary file to the final destination
            cp {output.tmp_ont} {output.final_ont}
            """)
        else:
            # Run clumpify on Illumina data (paired-end)
            shell("""
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
            """)