import os
import sys
import pandas as pd

def contains_non_phage_viruses(group):
    virus_hits = group[group['rank'].str.contains('superkingdom:Viruses')]
    non_phage_viruses = virus_hits[~virus_hits['stitle'].str.contains('phage', case=False)]
    return len(non_phage_viruses) > 0

def main(snakemake):
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    # Check if the input file is empty
    if os.path.getsize(input_file) == 0:
        # Create an empty output file if the input is empty
        open(output_file, 'w').close()
        print("Input file is empty. Created an empty output file.", file=sys.stderr)
    else:
        # Read the input file
        df = pd.read_csv(input_file, sep="\t")

        # Filter the dataframe
        filtered_df = df.groupby('qseqid').filter(contains_non_phage_viruses)

        # Write the filtered data to the output file
        filtered_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    main(snakemake)