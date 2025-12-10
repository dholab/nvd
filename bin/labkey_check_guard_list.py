#!/usr/bin/env python3

import argparse
from labkey.api_wrapper import APIWrapper
from loguru import logger
import sys


def experiment_exists(labkey, guard_list_name, experiment_id) -> bool:
    """
    Check if an experiment ID exists in the guard list.

    Returns:
        bool: True if experiment exists, False otherwise
    """
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

    return bool(result and result["rows"])


def check_experiment(labkey, guard_list_name, experiment_id) -> bool:
    """
    Check if the experiment ID is available (doesn't exist in guard list).

    Args:
        labkey: APIWrapper instance
        guard_list_name: Name of the LabKey guard list
        experiment_id: Experiment ID to check

    Returns:
        bool: True if experiment ID is available (not in list), False if already used
    """
    try:
        logger.info(f"Checking if experiment ID {experiment_id} exists in guard list")

        if experiment_exists(labkey, guard_list_name, experiment_id):
            logger.error(
                f"Experiment ID {experiment_id} already exists in guard list - data has already been uploaded"
            )
            return False
        else:
            logger.info(
                f"Experiment ID {experiment_id} is available - safe to proceed"
            )
            return True

    except Exception as e:
        error_msg = str(e).lower()
        if any(x in error_msg for x in ["does not exist", "table", "querynotfound", "not found"]):
            logger.warning(
                f"Guard list '{guard_list_name}' doesn't exist - assuming safe to proceed"
            )
            return True
        else:
            logger.error(f"Error checking guard list: {e}")
            return False


def register_experiment(labkey, guard_list_name, experiment_id) -> bool:
    """
    Register an experiment ID in the guard list after successful upload.

    Args:
        labkey: APIWrapper instance
        guard_list_name: Name of the LabKey guard list
        experiment_id: Experiment ID to register

    Returns:
        bool: True if registration succeeded, False otherwise
    """
    try:
        logger.info(f"Registering experiment ID {experiment_id} in guard list")

        # Double-check it doesn't already exist
        if experiment_exists(labkey, guard_list_name, experiment_id):
            logger.warning(
                f"Experiment ID {experiment_id} already exists - skipping registration"
            )
            return True

        # Insert new entry
        guard_row = {"experiment_id": experiment_id}
        insert_result = labkey.query.insert_rows(
            "lists", guard_list_name, [guard_row]
        )
        inserted_row = insert_result["rows"][0]
        logger.info(
            f"Successfully registered experiment ID {experiment_id} (Key={inserted_row.get('Key', 'N/A')})"
        )
        return True

    except Exception as e:
        error_msg = str(e).lower()
        if any(x in error_msg for x in ["does not exist", "table", "querynotfound", "not found"]):
            logger.warning(f"Guard list '{guard_list_name}' doesn't exist")
            logger.info("Attempting to create entry anyway...")
            try:
                guard_row = {"experiment_id": experiment_id}
                insert_result = labkey.query.insert_rows(
                    "lists", guard_list_name, [guard_row]
                )
                logger.info(f"Successfully created guard list entry for experiment ID {experiment_id}")
                return True
            except Exception as create_error:
                logger.error(f"Failed to create guard list entry: {create_error}")
                logger.error("Please create the guard list manually with column: experiment_id (integer)")
                return False
        else:
            logger.error(f"Error registering in guard list: {e}")
            return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check and manage LabKey experiment ID guard list"
    )
    parser.add_argument(
        "--mode",
        choices=["check", "register"],
        required=True,
        help="'check' to verify experiment ID is available, 'register' to add to guard list"
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
        "--api_key", required=True, help="API key with insert permissions"
    )
    parser.add_argument(
        "--experiment_id",
        type=int,
        required=True,
        help="Experiment ID to check/register",
    )

    args = parser.parse_args()

    labkey = APIWrapper(args.server, args.container, api_key=args.api_key)

    mode_label = "CHECK" if args.mode == "check" else "REGISTER"

    logger.info("=" * 80)
    logger.info(f"EXPERIMENT ID GUARD LIST - {mode_label}")
    logger.info("=" * 80)
    logger.info(f"Server: {args.server}")
    logger.info(f"Container: {args.container}")
    logger.info(f"Guard List: {args.guard_list}")
    logger.info(f"Experiment ID: {args.experiment_id}")
    logger.info("=" * 80)

    if args.mode == "check":
        success = check_experiment(labkey, args.guard_list, args.experiment_id)
    else:
        success = register_experiment(labkey, args.guard_list, args.experiment_id)

    if not success:
        logger.info("=" * 80)
        logger.info(f"{mode_label} FAILED")
        logger.info("=" * 80)
        if args.mode == "check":
            logger.error(f"Experiment ID {args.experiment_id} has already been uploaded to LabKey")
            logger.info("Options:")
            logger.info(f"   - Use a different experiment_id")
            logger.info(f"   - Remove the entry from guard list '{args.guard_list}' if re-upload is intended")
        sys.exit(1)
    else:
        logger.info("=" * 80)
        logger.info(f"{mode_label} SUCCESSFUL")
        logger.info("=" * 80)
        sys.exit(0)


if __name__ == "__main__":
    main()
