# NVD2: Identify Human Viruses (and more) in your Metagenomic Samples

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Further Documentation](#further-documentation)
- [Citation](#citation)

## Overview

NVD2 is a Nextflow pipeline focused primarily though not exclusively on finding
human viruses in metagenomic samples. It leverages the battle-tested NCBI
toolchain, including [stat]() and [BLAST](), to identify human viruses while
minimizing false positives and negatives. And for more general purpose
classification, NVD2 implements a [GOTTCHA2]() subworkflow as well, which uses a
carefully curated database of taxonomically diagnostic reference sequences.

NVD2 was designed from the ground up to handle enormous datasets and performs
particularly well with Illumina deep sequencing datasets from wastewater
sewersheds, the context for which it was originally designed. To perform well in
this context, it must:

1. Handle highly fragmented genomes.
2. Be resilient to wild fluctuations in depth-of-coverage.
3. Resolve ambiguities between closely related organisms with high sequence
   identity.

In addition to wastewater sequenced with Illumina, NVD2 also supports datasets
sequenced on Oxford Nanopore instruments and has been tested on Nanopore
libraries generated from air samples.

## Quick Start

A minimal NVD2 setup requires that [Nextflow](https://nextflow.io/) and
[Docker](https://www.docker.com/) are installed.

## Further Documentation

## Citation

Coming soon!

## Acknowledgments

- [Marc Johnson](https://medicine.missouri.edu/faculty/marc-johnson-phd) and
  [Shelby O'Connor](https://slo.pathology.wisc.edu), our partners in innovative
  pathogen monitoring from environmental samples.
- Kenneth Katz, NCBI, for developing NCBI STAT, maintaining pre-built databases
  for STAT, and helpful discusssions
- [C. Titus Brown](https://www.vetmed.ucdavis.edu/faculty/c-titus-brown), for
  helpful discussions of using kmer classifiers as part of metagenomic workflows
- Development funded by [Inkfish](https://ink.fish)

## License

See `LICENSE` for more information.
