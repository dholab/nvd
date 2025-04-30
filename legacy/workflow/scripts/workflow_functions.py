"""
Snakemake Workflow Utility Functions

This module contains utility functions used in a Snakemake workflow for
handling various input types, retrieving file paths, and managing remote
file operations.

Functions:
    get_input_type: Determines the input type for a given sample.
    get_input_files: Retrieves the input file paths for a given sample.
    get_sra_accession: Retrieves the SRA accession for a given sample.
    create_directory: Creates a directory on a remote server.
    upload_file: Uploads a file to a remote server.
"""

import time
import requests
from typing import List, Dict, Tuple, Any


def get_input_type(wildcards: Any, config: Dict[str, Any]) -> str:
    try:
        sample = next(s for s in config['samples'] if s['name'] == wildcards.sample)
    except StopIteration:
        available_samples = [s['name'] for s in config['samples']]
        raise ValueError(f"Sample '{wildcards.sample}' not found in config. Available samples: {available_samples}")

    if 'r1_fastq' in sample and 'r2_fastq' in sample:
        return 'paired_end'
    elif 'interleaved' in sample:
        return 'interleaved'
    elif 'ont' in sample:
        return 'ont'
    elif 'sra' in sample:
        return 'sra'
    else:
        raise ValueError(f"No valid input type found for sample {wildcards.sample}")

def get_input_files(wildcards: Any, config: Dict[str, Any]) -> List[str]:
    try:
        sample = next(s for s in config['samples'] if s['name'] == wildcards.sample)
    except StopIteration:
        available_samples = [s['name'] for s in config['samples']]
        raise ValueError(f"Sample '{wildcards.sample}' not found in config. Available samples: {available_samples}")

    if 'r1_fastq' in sample and 'r2_fastq' in sample:
        return [sample['r1_fastq'], sample['r2_fastq']]
    elif 'interleaved' in sample:
        return [sample['interleaved']]
    elif 'ont' in sample:
        return [sample['ont']]
    elif 'sra' in sample:
        return []
    else:
        raise ValueError(f"No valid input files found for sample {wildcards.sample}")

def get_sra_accession(wildcards: Any, config: Dict[str, Any]) -> str:
    sample = next((s for s in config['samples'] if s['name'] == wildcards.sample), None)
    return sample.get('sra', '') if sample else ''


def create_directory(url: str, auth: Tuple[str, str], max_retries: int = 3) -> bool:
    """
    Creates a directory on a remote server.

    Args:
        url (str): The URL of the directory to create.
        auth (Tuple[str, str]): Authentication tuple (username, password).
        max_retries (int): Maximum number of retry attempts.

    Returns:
        bool: True if the directory was created successfully, False otherwise.
    """
    for attempt in range(max_retries):
        try:
            response = requests.request('MKCOL', url, auth=auth)
            if response.status_code in [201, 405]:  # 201: Created, 405: Method Not Allowed (directory might already exist)
                return True
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                print(f"Failed to create directory {url} after {max_retries} attempts: {str(e)}")
                return False
            time.sleep(2 ** attempt)  # Exponential backoff
    return False


def upload_file(local_path: str, remote_path: str, auth: Tuple[str, str], max_retries: int = 3) -> bool:
    """
    Uploads a file to a remote server.

    Args:
        local_path (str): Path to the local file to upload.
        remote_path (str): URL of the remote location to upload the file to.
        auth (Tuple[str, str]): Authentication tuple (username, password).
        max_retries (int): Maximum number of retry attempts.

    Returns:
        bool: True if the file was uploaded successfully, False otherwise.
    """
    for attempt in range(max_retries):
        try:
            with open(local_path, 'rb') as file:
                response = requests.put(remote_path, data=file, auth=auth)
                response.raise_for_status()
            return True
        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                print(f"Failed to upload {local_path} after {max_retries} attempts: {str(e)}")
                return False
            time.sleep(2 ** attempt)  # Exponential backoff
    return False