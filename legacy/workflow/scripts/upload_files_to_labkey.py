import os
import urllib.request
import urllib.parse
import urllib.error
import base64
import yaml

def directory_exists(url, auth):
    req = urllib.request.Request(url, method='GET')
    req.add_header('Authorization', 'Basic ' + base64.b64encode(f"{auth[0]}:{auth[1]}".encode()).decode())
    
    try:
        with urllib.request.urlopen(req) as response:
            return response.getcode() == 200
    except urllib.error.HTTPError as e:
        if e.code == 403:
            print(f"HTTPError: 403 - Access Forbidden when checking directory exists at {url}")
            # Assume the directory exists if access is forbidden
            return True
        print(f"HTTPError: {e.code} - {e.reason} when checking directory exists at {url}")
        print(f"Server Response: {e.read().decode()}")
        raise e
    except urllib.error.URLError as e:
        print(f"URLError: {e.reason} when checking directory exists at {url}")
        raise e

def create_directory(url, auth):
    if directory_exists(url, auth):
        print(f"Directory already exists: {url}")
        return True
    
    req = urllib.request.Request(url, method='MKCOL')
    req.add_header('Authorization', 'Basic ' + base64.b64encode(f"{auth[0]}:{auth[1]}".encode()).decode())
    
    try:
        with urllib.request.urlopen(req) as response:
            return response.getcode() in (201, 204)
    except urllib.error.HTTPError as e:
        if e.code == 405:  # Method not allowed, directory might already exist
            return True
        print(f"Failed to create directory. HTTPError: {e.code} - {e.reason}")
        print(f"Server Response: {e.read().decode()}")
        return False
    except urllib.error.URLError as e:
        print(f"URLError: {e.reason} when creating directory at {url}")
        return False

def upload_file(local_path, remote_url, auth):
    with open(local_path, 'rb') as file:
        data = file.read()
        req = urllib.request.Request(remote_url, data=data, method='PUT')
        req.add_header('Authorization', 'Basic ' + base64.b64encode(f"{auth[0]}:{auth[1]}".encode()).decode())
        try:
            with urllib.request.urlopen(req) as response:
                return response.getcode() in (201, 204)
        except urllib.error.HTTPError:
            return False

def mask_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Mask sensitive information
    if 'global' in config:
        if 'labkey_api_key' in config['global']:
            config['global']['labkey_api_key'] = 'MASKED'
        if 'labkey_username' in config['global']:
            config['global']['labkey_username'] = 'MASKED'
        if 'labkey_password' in config['global']:
            config['global']['labkey_password'] = 'MASKED'
    
    return yaml.dump(config)

def main(snakemake):
    auth = (snakemake.params.username, snakemake.params.password)
    
    # Upload to LabKey/WebDAV if credentials are available
    if snakemake.params.webdav_url and snakemake.params.username and snakemake.params.password and snakemake.params.labkey_server and snakemake.params.project_name:
        # Create the main experiment directory
        main_dir = urllib.parse.urljoin(snakemake.params.webdav_url, f"{snakemake.params.experiment}/")
        if not create_directory(main_dir, auth):
            raise Exception(f"Failed to create main directory: {main_dir}")
        
        # Create the run ID subdirectory
        run_dir = urllib.parse.urljoin(main_dir, f"{snakemake.params.snakemake_run_id}/")
        if not create_directory(run_dir, auth):
            raise Exception(f"Failed to create run directory: {run_dir}")
        
        # Upload tarball
        remote_tarball_path = urllib.parse.urljoin(run_dir, f"{snakemake.wildcards.sample}.tar.zst")
        if not upload_file(snakemake.input.tarball, remote_tarball_path, auth):
            print(f"Failed to upload {snakemake.wildcards.sample}.tar.zst")
        
        # Upload masked config.yaml
        masked_config = mask_config(snakemake.input.config_file)
        remote_config_path = urllib.parse.urljoin(run_dir, "config.yaml")
        req = urllib.request.Request(remote_config_path, data=masked_config.encode(), method='PUT')
        req.add_header('Authorization', 'Basic ' + base64.b64encode(f"{auth[0]}:{auth[1]}".encode()).decode())
        try:
            with urllib.request.urlopen(req) as response:
                if response.getcode() not in (201, 204):
                    print("Failed to upload masked config.yaml")
        except urllib.error.HTTPError:
            print("Failed to upload masked config.yaml")

    # Ensure the tarball is always copied to final_output_path
    os.makedirs(f"{snakemake.params.final_output_path}/{snakemake.wildcards.sample}/{snakemake.params.snakemake_run_id}", exist_ok=True)
    os.system(f"cp {snakemake.input.tarball} {snakemake.params.final_output_path}/{snakemake.wildcards.sample}/{snakemake.params.snakemake_run_id}/")
    
    # Create done file
    with open(snakemake.output.done, 'w') as f:
        f.write("Process completed")

if __name__ == "__main__":
    main(snakemake)
