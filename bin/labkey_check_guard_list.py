#!/usr/bin/env python3

import argparse
from labkey.api_wrapper import APIWrapper
from loguru import logger
import sys


def check_guard_list(labkey, guard_list_name, experiment_id, unique_id) -> bool:
    """
    Check and manage the experiment ID guard list to prevent concurrent pipeline runs
    with the same experiment ID across different HTCondor clusters.

    Args:
        labkey: APIWrapper instance
        guard_list_name: Name of the LabKey guard list
        experiment_id: Experiment ID to check/register
        unique_id: Custom chosen uuid for each pipeline run and if not specified, fall back to nextflow session ID

    Returns:
        bool: True if experiment ID is safe to use, False if blocked by another cluster
    """
    try:
        logger.info(
            f"Checking guard list for experiment ID: {experiment_id} - {unique_id}"
        )

        # Try different approaches to query the data based on available API methods
        result = None

        # Method 1: Try with QueryFilter import
        try:
            from labkey.query import QueryFilter

            result = labkey.query.select_rows(
                schema_name="lists",
                query_name=guard_list_name,
                filter_array=[QueryFilter("experiment_id", experiment_id, "eq")],
            )
            logger.info("Used QueryFilter method successfully")
        except (ImportError, AttributeError) as e:
            logger.info(f"QueryFilter method not available: {e}")

        # Method 2: Try with filter parameter (some LabKey API versions)
        if result is None:
            try:
                result = labkey.query.select_rows(
                    schema_name="lists",
                    query_name=guard_list_name,
                    filter_array=[f"experiment_id~eq={experiment_id}"],
                )
                logger.info("Used filter string method successfully")
            except Exception as e:
                logger.info(f"Filter string method failed: {e}")

        # Method 3: Get all rows and filter in Python (fallback)
        if result is None:
            logger.info("Falling back to client-side filtering")
            result = labkey.query.select_rows(
                schema_name="lists", query_name=guard_list_name
            )
            # Filter the results in Python
            if result["rows"]:
                original_rows = result["rows"]
                filtered_rows = [
                    row
                    for row in original_rows
                    if str(row.get("experiment_id", "")) == str(experiment_id)
                ]
                result["rows"] = filtered_rows
                logger.info(
                    f"Filtered {len(original_rows)} rows down to {len(filtered_rows)} matching rows"
                )

        if result["rows"]:
            # Experiment ID exists, check if it's from the same cluster
            existing_row = result["rows"][0]
            existing_unique_id = existing_row["uuid"]

            logger.info(
                f"Found existing entry: experiment_id={existing_row['experiment_id']}, unique_id={existing_unique_id}"
            )

            if str(existing_unique_id) == str(unique_id):
                logger.info(
                    f"Experiment ID {experiment_id} already registered for this cluster ({unique_id}) - proceeding"
                )
                return True
            else:
                logger.error(
                    f"Experiment ID {experiment_id} is already in use by cluster {existing_unique_id}"
                )
                return False
        else:
            # Experiment ID doesn't exist, add it to the guard list
            logger.info(
                f"No existing entry found. Registering experiment ID {experiment_id} for cluster {unique_id}"
            )
            guard_row = {"experiment_id": experiment_id, "uuid": unique_id}

            insert_result = labkey.query.insert_rows(
                "lists", guard_list_name, [guard_row]
            )
            inserted_row = insert_result["rows"][0]
            logger.info(
                f"Successfully registered experiment ID {experiment_id} in guard list with Key={unique_id}"
            )
            return True

    except Exception as e:
        # If the guard list doesn't exist, try to handle gracefully
        error_msg = str(e).lower()
        if (
            "does not exist" in error_msg
            or "table" in error_msg
            or "querynotfound" in error_msg
            or "not found" in error_msg
        ):
            logger.warning(
                f"Guard list '{guard_list_name}' doesn't exist - attempting to create entry"
            )
            logger.info(
                "Note: Please ensure the guard list exists with columns: experiment_id (integer), unique_id (text)"
            )
            # Try to proceed anyway - the insert might create the structure
            try:
                guard_row = {"experiment_id": experiment_id, "uuid": unique_id}
                insert_result = labkey.query.insert_rows(
                    "lists", guard_list_name, [guard_row]
                )
                inserted_row = insert_result["rows"][0]
                logger.info(
                    f"Successfully created guard list entry for experiment ID {experiment_id} with Key={inserted_row.get('Key', 'N/A')}"
                )
                return True
            except Exception as create_error:
                logger.error(f"Failed to create guard list entry: {create_error}")
                logger.error(
                    "Please create the guard list manually in LabKey with columns:"
                )
                logger.error("  - experiment_id (integer)")
                logger.error("  - unique_id (text)")
                return False
        else:
            logger.error(f"Error checking guard list: {e}")
            logger.error("This could be due to:")
            logger.error("  - Network connectivity issues")
            logger.error("  - API key permissions")
            logger.error("  - LabKey server problems")
            logger.error("  - Incorrect guard list name or schema")
            logger.error(f"  - Guard list name used: '{guard_list_name}'")
            logger.error(f"  - Schema: 'lists'")
            return False


def main() -> None:
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Check HTCondor cluster guard list for experiment ID conflicts"
    )
    parser.add_argument("--server", required=True, help="LabKey server URL")
    parser.add_argument(
        "--container",
        required=True,
        help="LabKey container path (e.g., project folder)",
    )
    parser.add_argument(
        "--guard_list", required=True, help="Name of the LabKey guard list"
    )
    parser.add_argument(
        "--api_key", required=True, help="API key with insert/delete permissions"
    )
    parser.add_argument("--unique_id", required=True, help="unique ID")
    parser.add_argument(
        "--experiment_id",
        type=int,
        required=True,
        help="Experiment ID to check/register",
    )

    args = parser.parse_args()

    # Initialize LabKey API wrapper
    labkey = APIWrapper(args.server, args.container, api_key=args.api_key)

    logger.info("=" * 80)
    logger.info("🛡️  EXPERIMENT ID GUARD LIST CHECK")
    logger.info("=" * 80)
    logger.info(f"Server: {args.server}")
    logger.info(f"Container: {args.container}")
    logger.info(f"Guard List: {args.guard_list}")
    logger.info(f"Experiment ID: {args.experiment_id}")
    logger.info(f"Unique ID: {args.unique_id}")
    logger.info("=" * 80)

    # Check the guard list
    if not check_guard_list(
        labkey, args.guard_list, args.experiment_id, args.unique_id
    ):
        print()
        logger.info("=" * 80)
        logger.info("❌ EXPERIMENT ID GUARD CHECK FAILED")
        logger.info("=" * 80)
        print()
        logger.info(
            f"Error: Experiment ID '{args.experiment_id}' is already being/been processed by another HTCondor cluster",
        )
        print()
        logger.info("REQUIRED ACTION:")
        logger.info(
            "   This experiment ID is currently in use by another pipeline run."
        )
        logger.info(
            "   Please wait for the other run to complete or use a different experiment_id."
        )
        print()
        logger.info("Troubleshooting:")
        logger.info(
            "   • Check if another pipeline is currently running with this experiment ID"
        )
        logger.info("   • Verify that previous pipeline runs completed successfully")
        logger.info(
            "   • Consider using a new experiment ID if the previous run failed"
        )
        logger.info(
            f"   • Check the guard list '{args.guard_list}' in LabKey for active experiment IDs"
        )
        logger.info(
            "   • You may need to manually remove stale entries from the guard list"
        )
        print()
        logger.info("Configuration:")
        logger.info(
            f"   Current experiment_id: {args.experiment_id} (blocked by another cluster)"
        )
        logger.info(f"   Current unique_id: {args.unique_id}")
        print()
        logger.info("Guard List Management:")
        logger.info(
            f"   • To clear a stale entry, delete the row from the guard list ({args.guard_list}) in LabKey"
        )
        logger.info(
            "   • Guard list should have columns: experiment_id (integer), unique_id (text)"
        )
        logger.info(
            "   • Each experiment_id should only have one active entry at a time"
        )
        print()
        sys.exit(1)
    else:
        print()
        logger.info("=" * 80)
        logger.info("✅ GUARD LIST CHECK SUCCESSFUL")
        logger.info("=" * 80)
        print()
        logger.info(
            f"• Experiment ID {args.experiment_id} is registered for {args.unique_id}"
        )
        logger.info("• No conflicting unique pair found")
        logger.info("• Pipeline can proceed safely")
        logger.info(f"• Guard list '{args.guard_list}' updated successfully")
        print()
        logger.info("Next Steps:")
        logger.info("   • Pipeline will now proceed with sample processing")
        logger.info("   • All samples in this session will use the same experiment ID")
        logger.info("   • Guard list entry will prevent conflicts with other session")
        print()


if __name__ == "__main__":
    main()
