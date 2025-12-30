#!/usr/bin/env python3

import argparse

from labkey.api_wrapper import APIWrapper
from loguru import logger

# Define the expected schema for the BLAST hit lists
blast_fields = [
    "Experiment",
    "Blast Task",
    "Sample Id",
    "Qseqid",
    "Qlen",
    "Sseqid",
    "Stitle",
    "Tax Rank",
    "Length",
    "Pident",
    "Evalue",
    "Bitscore",
    "Blast Db Version",
    "Snakemake Run Id",
    "Sscinames",
    "Mapped Reads",
    "Stat Db Version",
    "Total Reads",
    "Exclude",
    "Adjusted Taxid",
    "Adjustment Method",
    "Adjusted Taxid Name",
    "Adjusted Taxid Rank",
    "Staxids",
]

# Define the expected schema for the FASTA output lists
blast_fasta_fields = [
    "Experiment",
    "Sample Id",
    "Contig Id",
    "Contig Sequence",
    "Notes",
    "Snakemake Run Id",
]

gottcha2_full_fields = [
    "SAMPLE",
    "EXPERIMENT",
    "LEVEL",
    "NAME",
    "TAXID",
    "READ COUNT",
    "TOTAL BP MAPPED",
    "TOTAL BP MISMATCH",
    "LINEAR LEN",
    "LINEAR DOC",
    "ROLLUP DOC",
    "REL ABUNDANCE",
    "LINEAR COV",
    "LINEAR COV MAPPED SIG",
    "BEST LINEAR COV",
    "MAPPED SIG LENGTH",
    "TOL SIG LENGTH",
    "ABUNDANCE",
    "ZSCORE",
    "NOTE",
    "GOTTCHA2 DB VERSION",
]

gottcha2_fasta_fields = [
    "Experiment",
    "Sample Id",
    "Header",
    "Sequence",
    "Notes",
    "Snakemake Run Id",
    "Upload Type",
]


def main() -> None:
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Validate LabKey list schema and permissions",
    )
    parser.add_argument("--server", required=True, help="LabKey server URL")
    parser.add_argument(
        "--container",
        required=True,
        help="LabKey container path (e.g., project folder)",
    )
    parser.add_argument(
        "--list",
        required=True,
        help="Name of the LabKey list to validate",
    )
    parser.add_argument(
        "--api_key",
        required=True,
        help="API key with insert/delete permissions",
    )
    parser.add_argument(
        "--type",
        choices=["blast", "blast_fasta", "gottcha2_full", "gottcha2_fasta"],
        required=True,
        help="Type of list to validate",
    )
    parser.add_argument(
        "--experiment_id",
        type=int,
        required=True,
        help="Experiment ID for validation testing",
    )
    args = parser.parse_args()

    # Initialize LabKey API wrapper
    labkey = APIWrapper(args.server, args.container, api_key=args.api_key)

    logger.info("=" * 80)
    logger.info("🔍 LABKEY LIST VALIDATION")
    logger.info("=" * 80)
    logger.info(f"Server: {args.server}")
    logger.info(f"Container: {args.container}")
    logger.info(f"List: {args.list}")
    logger.info(f"Type: {args.type}")
    logger.info(f"Experiment ID: {args.experiment_id}")
    logger.info("=" * 80)

    # Use a unique experiment ID that's unlikely to conflict with real data for dummy row
    unique_experiment_id = 123456789

    # Construct a dummy row with properly typed fields matching the expected schema
    if args.type == "blast":
        dummy = {
            "Experiment": unique_experiment_id,
            "Blast Task": "__DUMMY__TASK__",
            "Sample Id": "__DUMMY__SAMPLE__",
            "Qseqid": "__DUMMY__QSEQID__",
            "Qlen": 0,
            "Sseqid": "__DUMMY__SSEQID__",
            "Stitle": "__DUMMY__STITLE__",
            "Tax Rank": "__DUMMY__TAX__",
            "Length": 0,
            "Pident": 0.0,
            "Evalue": 0.0,
            "Bitscore": 0,
            "Blast Db Version": "unknown",
            "Snakemake Run Id": "__DUMMY__SMK__",
            "Sscinames": 0.0,
            "Mapped Reads": 0,
            "Stat Db Version": "unknown",
            "Total Reads": 0.0,
            "Exclude": False,
            "Adjusted Taxid": 11775,
            "Adjustment Method": "lca",
            "Adjusted Taxid Name": "Feline sarcoma virus (STRAIN GARDNER-RASHEED)",
            "Adjusted Taxid Rank": "no rank",
            "Staxids": 128951,
        }
    elif args.type == "gottcha2_full":
        dummy = {
            "SAMPLE": "__DUMMY__SAMPLE__",
            "EXPERIMENT": unique_experiment_id,
            "LEVEL": "__DUMMY__LEVEL__",
            "NAME": "__DUMMY__NAME__",
            "TAXID": 0,
            "READ COUNT": 0,
            "TOTAL BP MAPPED": 0,
            "TOTAL BP MISMATCH": 0,
            "LINEAR LEN": 0,
            "LINEAR DOC": 0.0,
            "ROLLUP DOC": 0.0,
            "REL ABUNDANCE": 0.0,
            "LINEAR COV": 0.0,
            "LINEAR COV MAPPED SIG": 0.0,
            "BEST LINEAR COV": 0.0,
            "MAPPED SIG LENGTH": 0,
            "TOL SIG LENGTH": 0,
            "ABUNDANCE": 0.0,
            "ZSCORE": 0.0,
            "NOTE": "__DUMMY__NOTE__",
            "GOTTCHA2 DB VERSION": "v0",
        }
    elif args.type == "gottcha2_fasta":
        dummy = {
            "Experiment": unique_experiment_id,
            "Sample Id": "__DUMMY__SAMPLE__",
            "Header": "__DUMMY__HEADER__",
            "Sequence": "ATGCATGC",
            "Notes": "__DUMMY__NOTES__",
            "Snakemake Run Id": "__DUMMY__SMK__",
            "Upload Type": "fasta",
        }
    elif args.type == "blast_fasta":
        dummy = {
            "Experiment": unique_experiment_id,
            "Sample Id": "__DUMMY__SAMPLE__",
            "Contig Id": "__DUMMY__CONTIG__",
            "Contig Sequence": "ATGCATGC",
            "Notes": 0.0,
            "Snakemake Run Id": "__DUMMY__SMK__",
        }

    try:
        logger.info(f"Testing insert/delete operations on list '{args.list}'...")

        # Attempt to insert the dummy row into the specified LabKey list
        insert_result = labkey.query.insert_rows("lists", args.list, [dummy])
        inserted_row = insert_result["rows"][0]

        # Get the Key field if it exists for logging
        row_key = inserted_row.get("Key", "N/A")
        logger.info(
            f"[OK] Inserted dummy row with Key={row_key}, Experiment={unique_experiment_id}",
        )

        # Delete the dummy row using the complete inserted row data
        # LabKey requires the full row data (including auto-generated fields like Key)
        delete_result = labkey.query.delete_rows("lists", args.list, [inserted_row])
        deleted_key = delete_result["rows"][0].get("Key", "N/A")
        logger.info(f"[OK] Deleted dummy row successfully with Key={deleted_key}")

        print()
        logger.info("=" * 80)
        logger.info("✅ LIST VALIDATION SUCCESSFUL")
        logger.info("=" * 80)
        print()
        logger.info("• LabKey connection and permissions validated")
        logger.info(f"• List '{args.list}' is accessible and writable")
        logger.info(f"• Schema validation passed for type '{args.type}'")
        logger.info("• Insert and delete operations working correctly")
        print()

    except Exception as e:
        logger.info(f"[FAIL] Insert/Delete failed: {e}")
        print()
        logger.info("=" * 80)
        logger.info("⚠️  LIST VALIDATION FAILED - Troubleshooting Guide:")
        logger.info("=" * 80)
        print()
        logger.info("Common issues and solutions:")
        print()
        logger.info("1. API Key Issues:")
        logger.info(
            "   • Check that LABKEY_API_KEY secret is configured: nvd secrets set",
        )
        logger.info("   • Verify the API key has INSERT and DELETE permissions")
        logger.info("   • Ensure the API key hasn't expired")
        logger.info("   • Test the API key directly in LabKey's web interface")
        print()
        logger.info("2. Server/Container Issues:")
        logger.info(f"   • Verify server URL is correct: {args.server}")
        logger.info(f"   • Check container/project path: {args.container}")
        logger.info("   • Ensure you have access to the specified project")
        logger.info("   • Confirm the server URL includes 'https://' if required")
        print()
        logger.info("3. List Issues:")
        logger.info(f"   • Verify list '{args.list}' exists in LabKey")
        logger.info("   • Check that the list has the expected schema/columns")
        logger.info("   • Ensure you have read/write permissions to the list")
        logger.info("   • Confirm list is in the correct folder/container")
        print()
        logger.info("4. Schema Issues:")
        logger.info(f"   • Verify list schema matches expected '{args.type}' format")
        logger.info("   • Check column names and data types")
        logger.info("   • Ensure required fields are not missing")
        print()
        logger.info("5. Network/Connectivity Issues:")
        logger.info("   • Check firewall or proxy settings")
        logger.info("   • Verify network connectivity to the LabKey server")
        logger.info("   • Try accessing the LabKey server through a web browser")
        print()
        logger.info("6. Configuration Issues:")
        logger.info("   • Double-check your nextflow.config parameters:")
        logger.info(f"     - labkey_server = '{args.server}'")
        logger.info(f"     - labkey_project_name = '{args.container}'")
        logger.info(f"     - labkey_{args.type}_list = '{args.list}'")
        print()
        logger.info("For more help:")
        logger.info("   • Check LabKey server logs for detailed error messages")
        logger.info("   • Verify your LabKey user permissions with an administrator")
        logger.info("   • Test API access using LabKey's built-in API documentation")
        print()
        exit(1)


if __name__ == "__main__":
    main()
