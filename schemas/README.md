# NVD JSON Schemas

This directory contains JSON Schema definitions for the NVD pipeline.

## Available Schemas

| Schema | Description |
|--------|-------------|
| `nvd-params.v2.3.0.schema.json` | Pipeline parameters schema (version 2.3.0) |
| `nvd-params.v2.4.0.schema.json` | Pipeline parameters schema (version 2.4.0) |
| `nvd-params.v2.5.0.schema.json` | Pipeline parameters schema (version 2.5.0) |
| `nvd-params.latest.schema.json` | Symlink to the current version |

## Usage

### VS Code (YAML)

Add a schema reference comment at the top of your params file:

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/dhoconno/nvd/main/schemas/nvd-params.latest.schema.json

samplesheet: samples.csv
experiment_id: 1001
tools: all
cutoff_percent: 0.001
```

### VS Code (JSON)

Include the `$schema` property in your JSON params file:

```json
{
  "$schema": "https://raw.githubusercontent.com/dhoconno/nvd/main/schemas/nvd-params.latest.schema.json",
  "samplesheet": "samples.csv",
  "experiment_id": 1001,
  "tools": "all",
  "cutoff_percent": 0.001
}
```

### Benefits

With the schema reference, your IDE provides:

- **Autocomplete** for parameter names
- **Type validation** (red squiggles for wrong types)
- **Hover documentation** (descriptions from schema)
- **Enum suggestions** (e.g., valid `tools` values)

### Using with Nextflow

Pass your params file to Nextflow with the `-params-file` option:

```bash
nextflow run dhoconno/nvd -params-file my-params.yaml -profile docker
```

Or use the NVD CLI wrapper:

```bash
nvd run -s samples.csv -e exp001 --config my-params.yaml
```

## Schema Versioning

Schema versions are tied to NVD pipeline versions. The `latest` symlink always
points to the current version for users who want to track updates automatically.

For reproducibility, you can reference a specific version:

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/dhoconno/nvd/main/schemas/nvd-params.v2.4.0.schema.json
```

## Validation

Validate a params file against the schema using the NVD CLI:

```bash
nvd params validate my-params.yaml
```

Or generate a template params file with the schema reference pre-configured:

```bash
nvd params init my-params.yaml
```
