#!/usr/bin/env python3
"""
Simple functional WebDAV client for LabKey
Usage: python webdav_client.py
"""

import urllib.request
import urllib.parse
import urllib.error
import base64
from functools import partial, wraps
from typing import Optional, Callable, Dict, List
import os
import sys


# Pure functions for authentication
def create_auth_header(username: str, password: str) -> str:
    """Create Basic authentication header value"""
    credentials = f"{username}:{password}"
    encoded = base64.b64encode(credentials.encode()).decode()
    return f"Basic {encoded}"


def create_request(
    url: str,
    method: str,
    auth_header: str,
    data: Optional[bytes] = None,
    headers: Optional[dict] = None,
) -> urllib.request.Request:
    """Create a urllib Request object with authentication"""
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", auth_header)
    if headers:
        for key, value in headers.items():
            req.add_header(key, value)
    return req


# Higher-order function for error handling
def with_error_handling(error_handler: Callable) -> Callable:
    """Decorator to handle HTTP errors functionally"""

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except (urllib.error.HTTPError, urllib.error.URLError) as e:
                return error_handler(e, *args, **kwargs)

        return wrapper

    return decorator


# Error handlers
def handle_exists_error(e: Exception, url: str, auth_header: str) -> bool:
    """Handle errors when checking if directory exists"""
    if isinstance(e, urllib.error.HTTPError):
        if e.code == 404:
            return False  # Directory doesn't exist
        elif e.code == 403:
            print(f"Access forbidden for {url} - assuming it exists")
            return True  # Might exist but no permission
        else:
            print(f"HTTP {e.code} error checking {url}")
            return False
    print(f"URL error checking {url}: {e}")
    return False


def handle_create_error(e: Exception, url: str, auth_header: str) -> bool:
    """Handle errors when creating directory"""
    if isinstance(e, urllib.error.HTTPError):
        if e.code == 405:
            print(f"Directory might already exist: {url}")
            return True
        elif e.code == 409:
            print(f"Parent directory doesn't exist for: {url}")
            return False
        else:
            print(f"Failed to create directory. HTTP {e.code}: {e.reason}")
            return False
    print(f"URL error creating directory at {url}: {e}")
    return False


def handle_upload_error(
    e: Exception, local_path: str, remote_url: str, auth_header: str
) -> bool:
    """Handle errors when uploading file"""
    if isinstance(e, urllib.error.HTTPError):
        print(f"Failed to upload {local_path}. HTTP {e.code}: {e.reason}")
        if e.code == 409:
            print("Parent directory might not exist")
    else:
        print(f"URL error uploading to {remote_url}: {e}")
    return False


# Core WebDAV operations
@with_error_handling(handle_exists_error)
def check_exists(url: str, auth_header: str) -> bool:
    """Check if a WebDAV resource exists"""
    # Use GET method like the working code, not PROPFIND
    req = create_request(url, "GET", auth_header)
    with urllib.request.urlopen(req) as response:
        return response.getcode() == 200


@with_error_handling(handle_create_error)
def create_collection(url: str, auth_header: str) -> bool:
    """Create a WebDAV collection (directory)"""
    req = create_request(url, "MKCOL", auth_header)
    with urllib.request.urlopen(req) as response:
        success = response.getcode() in (201, 204)
        if success:
            print(f"Created directory: {url}")
        return success


@with_error_handling(handle_upload_error)
def put_file(local_path: str, remote_url: str, auth_header: str) -> bool:
    """Upload a file to WebDAV"""
    if not os.path.exists(local_path):
        print(f"Local file not found: {local_path}")
        return False

    with open(local_path, "rb") as f:
        data = f.read()

    # Simplified headers like the working code
    req = create_request(remote_url, "PUT", auth_header, data=data)

    with urllib.request.urlopen(req) as response:
        success = response.getcode() in (200, 201, 204)
        if success:
            print(f"Uploaded: {local_path} -> {remote_url}")
        return success


# Utility functions
def ensure_trailing_slash(url: str) -> str:
    """Ensure URL ends with slash for directories"""
    return url if url.endswith("/") else url + "/"


def get_parent_url(url: str) -> str:
    """Get parent directory URL"""
    return url.rstrip("/").rsplit("/", 1)[0] + "/"


def is_root_webdav(url: str) -> bool:
    """Check if URL is the WebDAV root"""
    return "/_webdav/" in url and url.rstrip("/").endswith("/_webdav")


# Recursive directory creation
def create_directory_recursive(url: str, auth_header: str) -> bool:
    """Recursively create directory and its parents"""
    url = ensure_trailing_slash(url)

    # Check if already exists
    if check_exists(url, auth_header):
        print(f"Directory already exists: {url}")
        return True

    # Create parent if needed (recursive)
    parent = get_parent_url(url)
    if parent != url and not is_root_webdav(parent):
        if not create_directory_recursive(parent, auth_header):
            return False

    # Create the directory
    return create_collection(url, auth_header)


# High-level operations
def upload_files_to_directory(
    directory_url: str, file_paths: List[str], auth_header: str, create_dir: bool = True
) -> Dict[str, bool]:
    """Upload multiple files to a directory"""
    results = {}
    directory_url = ensure_trailing_slash(directory_url)

    # Create directory if requested
    if create_dir and not create_directory_recursive(directory_url, auth_header):
        print(f"Failed to create directory: {directory_url}")
        return results

    # Upload each file
    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            results[file_path] = False
            continue

        filename = os.path.basename(file_path)
        remote_url = urllib.parse.urljoin(directory_url, filename)
        results[file_path] = put_file(file_path, remote_url, auth_header)

    return results


# Simple WebDAV client factory
def create_client(base_url: str, username: str, password: str):
    """Create a simple WebDAV client"""
    auth_header = create_auth_header(username, password)
    base = ensure_trailing_slash(base_url)

    return {
        "check": partial(check_exists, auth_header=auth_header),
        "mkdir": partial(create_directory_recursive, auth_header=auth_header),
        "upload": partial(put_file, auth_header=auth_header),
        "upload_to_dir": partial(upload_files_to_directory, auth_header=auth_header),
        "base_url": base,
    }


# Main example
def main():
    """Simple usage examples"""

    # Configuration
    BASE_URL = "https://placeholder"
    USERNAME = "apikey"
    API_KEY = "<APIKEY>"

    # Create client
    client = create_client(BASE_URL, USERNAME, API_KEY)

    print("Simple WebDAV Client Examples\n" + "=" * 40)

    # Example 1: Create a directory
    dir_url = urllib.parse.urljoin(client["base_url"], "my_folder/subfolder/")
    print(f"\n1. Creating directory: {dir_url}")
    success = client["mkdir"](dir_url)
    print(f"   Result: {'Success' if success else 'Failed'}")

    # Example 2: Upload a single file
    test_file = "test.txt"
    if not os.path.exists(test_file):
        with open(test_file, "w") as f:
            f.write("Test content\n")

    file_url = urllib.parse.urljoin(client["base_url"], "my_folder/test.txt")
    print(f"\n2. Uploading file: {test_file} -> {file_url}")
    success = client["upload"](test_file, file_url)
    print(f"   Result: {'Success' if success else 'Failed'}")

    # Example 3: Upload multiple files to a directory
    files = ["file1.txt", "file2.txt", "file3.txt"]
    for fname in files:
        if not os.path.exists(fname):
            with open(fname, "w") as f:
                f.write(f"Content of {fname}\n")

    target_dir = urllib.parse.urljoin(client["base_url"], "batch_upload/")
    print(f"\n3. Uploading multiple files to: {target_dir}")
    results = client["upload_to_dir"](target_dir, files)
    for file_path, success in results.items():
        print(f"   {'✓' if success else '✗'} {file_path}")

    # Example 4: Check if paths exist
    print(f"\n4. Checking paths:")
    paths = [
        urllib.parse.urljoin(client["base_url"], "my_folder/"),
        urllib.parse.urljoin(client["base_url"], "my_folder/test.txt"),
        urllib.parse.urljoin(client["base_url"], "nonexistent/"),
    ]
    for path in paths:
        exists = client["check"](path)
        print(f"   {path}: {'exists' if exists else 'not found'}")


# Command-line interface
def cli():
    """
    CLI commands

    ## mkdir

    Create a directory or many directories
    >python3 bin/webdav_CLIent.py --server "https://Placeholder" --password XXXX mkdir foo

    --Creates foo dir
    --It will not overwrite an existing directory

    ## upload

    Upload a single file
    >python3 bin/webdav_CLIent.py --server "https://Placeholder" --password XXXX upload sample_krona.html my_folder/sample_krona.html

    --sample_krona.html will be copied from your local dir and placed at that path
    --If the path does not exist it will automatically create the parent dir

    --|my_folder
    --|--sample_krona.html

    ## upload-dir

    Create a dir and upload files
    >python3 bin/webdav_CLIent.py --server "https://Placeholder" --password XXXX upload-dir test_folder test_seq.fasta

    -- test_folder will be the directory and test_seq.fasta will be uploaded with the same name in test_folder
    -- |-testfolder
       |--test_seq.fasta

    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Simple WebDAV client for LabKey",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create directory
  python %(prog)s mkdir my_folder/subfolder/
  
  # Upload single file
  python %(prog)s upload local_file.txt remote_path/file.txt
  
  # Upload multiple files to directory
  python %(prog)s upload-dir target_folder/ file1.txt file2.csv file3.pdf
  
  # Check if path exists
  python %(prog)s check some_folder/file.txt
        """,
    )

    # Authentication
    parser.add_argument(
        "--server", default="https://Placeholder", help="WebDAV server URL"
    )
    parser.add_argument("--user", default="apikey", help="Username")
    parser.add_argument("--password", required=True, help="Password or API key")

    # Commands
    parser.add_argument(
        "command",
        choices=["mkdir", "upload", "upload-dir", "check"],
        help="Command to execute",
    )
    parser.add_argument("args", nargs="+", help="Command arguments")

    args = parser.parse_args()

    # Create client
    client = create_client(args.server, args.user, args.password)

    # Execute command
    if args.command == "mkdir":
        for path in args.args:
            url = urllib.parse.urljoin(client["base_url"], path)
            success = client["mkdir"](url)
            if not success:
                sys.exit(1)

    elif args.command == "upload":
        if len(args.args) != 2:
            print("Error: upload requires exactly 2 arguments: local_file remote_path")
            sys.exit(1)
        local_file, remote_path = args.args
        remote_url = urllib.parse.urljoin(client["base_url"], remote_path)

        # Debug output
        print(f"DEBUG - Local file: {local_file}")
        print(f"DEBUG - Remote path: {remote_path}")
        print(f"DEBUG - Base URL: {client['base_url']}")
        print(f"DEBUG - Full remote URL: {remote_url}")
        print(f"DEBUG - File exists: {os.path.exists(local_file)}")

        success = client["upload"](local_file, remote_url)
        sys.exit(0 if success else 1)

    elif args.command == "upload-dir":
        if len(args.args) < 2:
            print(
                "Error: upload-dir requires at least 2 arguments: target_dir file1 [file2 ...]"
            )
            sys.exit(1)
        target_dir = args.args[0]
        files = args.args[1:]
        dir_url = urllib.parse.urljoin(client["base_url"], target_dir)
        results = client["upload_to_dir"](dir_url, files)
        success = all(results.values())
        sys.exit(0 if success else 1)

    elif args.command == "check":
        for path in args.args:
            url = urllib.parse.urljoin(client["base_url"], path)
            exists = client["check"](url)
            print(f"{'✓ Exists' if exists else '✗ Not found'}: {url}")


def debug_test():
    """Debug test to compare main() vs CLI auth"""
    import sys

    BASE_URL = "https://PLACEHOLDER"
    USERNAME = "apikey"
    API_KEY = "<APIKEY>"

    print("Debug Test - Comparing auth methods\n" + "=" * 40)

    # Method 1: Like main() - hardcoded
    client1 = create_client(BASE_URL, USERNAME, API_KEY)
    print(f"Client1 base_url: {client1['base_url']}")

    # Method 2: Like CLI - from args
    if len(sys.argv) > 1 and sys.argv[1] == "--password":
        password_from_cli = sys.argv[2]
        client2 = create_client(BASE_URL, USERNAME, password_from_cli)
        print(f"Client2 base_url: {client2['base_url']}")

        # Compare auth headers
        auth1 = create_auth_header(USERNAME, API_KEY)
        auth2 = create_auth_header(USERNAME, password_from_cli)
        print(f"\nAuth headers match: {auth1 == auth2}")
        print(
            f"Password lengths - Hardcoded: {len(API_KEY)}, CLI: {len(password_from_cli)}"
        )

        # Try upload with both
        test_file = "test_debug.txt"
        if not os.path.exists(test_file):
            with open(test_file, "w") as f:
                f.write("Debug test\n")

        print("\nTesting upload with hardcoded client:")
        result1 = client1["upload"](
            test_file, urllib.parse.urljoin(BASE_URL, "debug1.txt")
        )
        print(f"Result: {'Success' if result1 else 'Failed'}")

        print("\nTesting upload with CLI client:")
        result2 = client2["upload"](
            test_file, urllib.parse.urljoin(BASE_URL, "debug2.txt")
        )
        print(f"Result: {'Success' if result2 else 'Failed'}")
    else:
        print("Run with: python3 webdav_CLIent.py --password YOUR_KEY")


if __name__ == "__main__":
    # Debug mode - uncomment to test
    # debug_test()
    # sys.exit()

    # Normal CLI mode
    cli()
