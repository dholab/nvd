#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

import logging
import urllib.request

TAXONOMY_URL = "https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/gettax.sqlite"
PROGRESS_GRANULARITY = 100000

# Set up logging
logger = logging.getLogger("tax_analysis")


def main() -> None:
    logger.info(f"Downloading taxonomy file from {TAXONOMY_URL}")
    urllib.request.urlretrieve(TAXONOMY_URL, "gettax.sqlite")  # noqa: S310
    logger.info("Downloaded taxonomy file to 'gettax.sqlite'")


if __name__ == "__main__":
    main()
