"""Tests for py_nvd.models module, specifically NvdParams."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from py_nvd.models import (
    DEFAULT_HUMAN_VIRUS_FAMILIES,
    NvdParams,
    ParamSource,
    TracedParams,
    VALID_TOOLS,
    trace_merge,
)
from py_nvd.params import load_params_file


class TestNvdParamsInstantiation:
    """Tests for basic NvdParams instantiation."""

    def test_minimal_instantiation(self):
        """Can create NvdParams with no arguments (all defaults)."""
        p = NvdParams()
        assert p.tools is None
        assert p.preprocess is False
        assert p.cutoff_percent == 0.001

    def test_with_required_fields(self):
        """Can create with typical required fields."""
        p = NvdParams(
            samplesheet=Path("/tmp/samples.csv"),
            tools="blast",
            experiment_id="exp001",
        )
        assert p.samplesheet == Path("/tmp/samples.csv")
        assert p.tools == "blast"
        assert p.experiment_id == "exp001"

    def test_all_fields(self):
        """Can create with all fields specified."""
        p = NvdParams(
            # Core
            samplesheet=Path("/tmp/samples.csv"),
            results=Path("/tmp/results"),
            tools="blast,gottcha",
            experiment_id="exp001",
            max_concurrent_downloads=5,
            cleanup=True,
            work_dir=Path("/tmp/work"),
            # Database versions
            gottcha2_db_version="RefSeq-r220",
            blast_db_version="core-nt_2025-01-01",
            # Database paths
            gottcha2_db=Path("/db/gottcha2"),
            blast_db=Path("/db/blast"),
            blast_db_prefix="nt",
            # Preprocessing
            preprocess=True,
            dedup=True,
            trim_adapters=True,
            # Analysis
            cutoff_percent=0.01,
            entropy=0.85,
            tax_stringency=0.8,
            # LabKey
            labkey=True,
            labkey_server="example.com",
        )
        assert p.tools == "blast,gottcha"
        assert p.cleanup is True
        assert p.cutoff_percent == 0.01
        assert p.labkey is True


class TestNvdParamsToolsValidator:
    """Tests for tools field validation."""

    def test_valid_single_tool(self):
        """Single valid tool is accepted."""
        p = NvdParams(tools="blast")
        assert p.tools == "blast"

    def test_valid_multiple_tools(self):
        """Comma-separated valid tools are accepted."""
        p = NvdParams(tools="blast,gottcha,stat")
        assert p.tools == "blast,gottcha,stat"

    def test_valid_all_tools(self):
        """All valid tool names are accepted."""
        for tool in VALID_TOOLS:
            p = NvdParams(tools=tool)
            assert p.tools == tool

    def test_none_tools_allowed(self):
        """None is a valid value for tools."""
        p = NvdParams(tools=None)
        assert p.tools is None

    def test_invalid_tool_rejected(self):
        """Invalid tool name raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(tools="invalid_tool")
        assert "Invalid tools" in str(exc_info.value)

    def test_mixed_valid_invalid_rejected(self):
        """Mix of valid and invalid tools is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(tools="blast,invalid,gottcha")
        assert "Invalid tools" in str(exc_info.value)
        assert "invalid" in str(exc_info.value)


class TestNvdParamsRangeValidators:
    """Tests for 0-1 range validators."""

    def test_cutoff_percent_valid_range(self):
        """cutoff_percent accepts 0-1 range."""
        assert NvdParams(cutoff_percent=0.0).cutoff_percent == 0.0
        assert NvdParams(cutoff_percent=0.5).cutoff_percent == 0.5
        assert NvdParams(cutoff_percent=1.0).cutoff_percent == 1.0

    def test_cutoff_percent_above_one_rejected(self):
        """cutoff_percent > 1 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(cutoff_percent=1.5)
        assert "Must be between 0 and 1" in str(exc_info.value)

    def test_cutoff_percent_below_zero_rejected(self):
        """cutoff_percent < 0 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(cutoff_percent=-0.1)
        assert "Must be between 0 and 1" in str(exc_info.value)

    def test_entropy_valid_range(self):
        """entropy accepts 0-1 range."""
        assert NvdParams(entropy=0.0).entropy == 0.0
        assert NvdParams(entropy=0.9).entropy == 0.9
        assert NvdParams(entropy=1.0).entropy == 1.0

    def test_entropy_out_of_range_rejected(self):
        """entropy outside 0-1 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(entropy=1.1)
        with pytest.raises(ValidationError):
            NvdParams(entropy=-0.5)

    def test_tax_stringency_valid_range(self):
        """tax_stringency accepts 0-1 range."""
        assert NvdParams(tax_stringency=0.0).tax_stringency == 0.0
        assert NvdParams(tax_stringency=0.7).tax_stringency == 0.7
        assert NvdParams(tax_stringency=1.0).tax_stringency == 1.0

    def test_tax_stringency_out_of_range_rejected(self):
        """tax_stringency outside 0-1 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(tax_stringency=2.0)


class TestNvdParamsPositiveIntValidators:
    """Tests for positive integer validators."""

    def test_min_gottcha_reads_valid(self):
        """min_gottcha_reads accepts positive integers."""
        assert NvdParams(min_gottcha_reads=1).min_gottcha_reads == 1
        assert NvdParams(min_gottcha_reads=250).min_gottcha_reads == 250

    def test_min_gottcha_reads_zero_rejected(self):
        """min_gottcha_reads=0 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(min_gottcha_reads=0)
        assert "Must be >= 1" in str(exc_info.value)

    def test_min_gottcha_reads_negative_rejected(self):
        """Negative min_gottcha_reads raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(min_gottcha_reads=-5)

    def test_max_blast_targets_valid(self):
        """max_blast_targets accepts positive integers."""
        assert NvdParams(max_blast_targets=1).max_blast_targets == 1
        assert NvdParams(max_blast_targets=100).max_blast_targets == 100

    def test_max_blast_targets_zero_rejected(self):
        """max_blast_targets=0 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(max_blast_targets=0)

    def test_min_read_length_valid(self):
        """min_read_length accepts positive integers."""
        assert NvdParams(min_read_length=1).min_read_length == 1
        assert NvdParams(min_read_length=50).min_read_length == 50

    def test_min_read_length_zero_rejected(self):
        """min_read_length=0 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(min_read_length=0)


class TestNvdParamsNonNegativeIntValidators:
    """Tests for non-negative integer validators."""

    def test_min_read_quality_illumina_zero_allowed(self):
        """min_read_quality_illumina=0 is valid."""
        p = NvdParams(min_read_quality_illumina=0)
        assert p.min_read_quality_illumina == 0

    def test_min_read_quality_illumina_positive_allowed(self):
        """Positive min_read_quality_illumina is valid."""
        p = NvdParams(min_read_quality_illumina=20)
        assert p.min_read_quality_illumina == 20

    def test_min_read_quality_illumina_negative_rejected(self):
        """Negative min_read_quality_illumina raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(min_read_quality_illumina=-1)
        assert "Must be >= 0" in str(exc_info.value)

    def test_min_read_quality_nanopore_zero_allowed(self):
        """min_read_quality_nanopore=0 is valid."""
        p = NvdParams(min_read_quality_nanopore=0)
        assert p.min_read_quality_nanopore == 0

    def test_min_read_quality_nanopore_negative_rejected(self):
        """Negative min_read_quality_nanopore raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(min_read_quality_nanopore=-5)


class TestNvdParamsMaxReadLength:
    """Tests for max_read_length validator (optional positive int)."""

    def test_none_allowed(self):
        """max_read_length=None is valid (no limit)."""
        p = NvdParams(max_read_length=None)
        assert p.max_read_length is None

    def test_positive_allowed(self):
        """Positive max_read_length is valid."""
        p = NvdParams(max_read_length=500)
        assert p.max_read_length == 500

    def test_zero_rejected(self):
        """max_read_length=0 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(max_read_length=0)
        assert "Must be >= 1" in str(exc_info.value)

    def test_negative_rejected(self):
        """Negative max_read_length raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(max_read_length=-100)


class TestNvdParamsTypeCoercion:
    """Tests for Pydantic type coercion."""

    def test_string_to_path(self):
        """String paths are coerced to Path objects."""
        p = NvdParams(samplesheet="/tmp/samples.csv")
        assert isinstance(p.samplesheet, Path)
        assert p.samplesheet == Path("/tmp/samples.csv")

    def test_string_to_int(self):
        """String integers are coerced to int."""
        p = NvdParams(min_gottcha_reads="100")
        assert isinstance(p.min_gottcha_reads, int)
        assert p.min_gottcha_reads == 100

    def test_string_to_float(self):
        """String floats are coerced to float."""
        p = NvdParams(cutoff_percent="0.05")
        assert isinstance(p.cutoff_percent, float)
        assert p.cutoff_percent == 0.05

    def test_int_to_float(self):
        """Integers are coerced to float where expected."""
        p = NvdParams(entropy=1)
        assert isinstance(p.entropy, float)
        assert p.entropy == 1.0


class TestNvdParamsMerge:
    """Tests for NvdParams.merge() class method."""

    def test_merge_empty_sources(self):
        """Merge with no sources returns defaults."""
        p = NvdParams.merge()
        assert p.tools is None
        assert p.cutoff_percent == 0.001

    def test_merge_single_dict(self):
        """Merge with single dict source."""
        p = NvdParams.merge({"tools": "blast", "cutoff_percent": 0.01})
        assert p.tools == "blast"
        assert p.cutoff_percent == 0.01

    def test_merge_precedence_later_wins(self):
        """Later sources override earlier sources."""
        preset = {"tools": "all", "cutoff_percent": 0.01}
        cli = {"tools": "blast"}
        p = NvdParams.merge(preset, cli)
        assert p.tools == "blast"  # CLI wins
        assert p.cutoff_percent == 0.01  # From preset

    def test_merge_three_sources(self):
        """Three-way merge with correct precedence."""
        preset = {"tools": "all", "cutoff_percent": 0.01, "entropy": 0.8}
        params_file = {"tools": "blast", "entropy": 0.85}
        cli = {"cutoff_percent": 0.005}
        p = NvdParams.merge(preset, params_file, cli)
        assert p.tools == "blast"  # From params_file
        assert p.cutoff_percent == 0.005  # From CLI
        assert p.entropy == 0.85  # From params_file

    def test_merge_none_sources_skipped(self):
        """None sources are gracefully skipped."""
        p = NvdParams.merge(None, {"tools": "blast"}, None)
        assert p.tools == "blast"

    def test_merge_none_values_skipped(self):
        """None values within dicts don't override."""
        preset = {"tools": "all", "cutoff_percent": 0.01}
        cli = {"tools": None, "cutoff_percent": 0.005}
        p = NvdParams.merge(preset, cli)
        assert p.tools == "all"  # None didn't override
        assert p.cutoff_percent == 0.005  # Non-None did override

    def test_merge_with_nvdparams_instance(self):
        """Can merge NvdParams instances, not just dicts."""
        base = NvdParams(tools="all", cutoff_percent=0.01)
        override = {"tools": "blast"}
        p = NvdParams.merge(base, override)
        assert p.tools == "blast"
        assert p.cutoff_percent == 0.01

    def test_merge_validates_result(self):
        """Merged result is validated."""
        with pytest.raises(ValidationError):
            NvdParams.merge({"tools": "invalid_tool"})


class TestNvdParamsToNextflowArgs:
    """Tests for NvdParams.to_nextflow_args() method."""

    def test_basic_command_structure(self):
        """Command starts with nextflow run <pipeline>."""
        p = NvdParams()
        cmd = p.to_nextflow_args(Path("/path/to/pipeline"))
        assert cmd[0] == "nextflow"
        assert cmd[1] == "run"
        assert cmd[2] == "/path/to/pipeline"

    def test_params_as_double_dash_args(self):
        """Params are formatted as --param_name value (underscores for Nextflow)."""
        p = NvdParams(tools="blast", cutoff_percent=0.01)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--tools" in cmd
        assert "blast" in cmd
        assert "--cutoff_percent" in cmd
        assert "0.01" in cmd

    def test_underscores_preserved_for_nextflow(self):
        """Param names keep underscores (Nextflow/Groovy requires them)."""
        p = NvdParams(min_gottcha_reads=100)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--min_gottcha_reads" in cmd
        assert "--min-gottcha-reads" not in cmd

    def test_bool_to_string(self):
        """Booleans are converted to 'true'/'false' strings."""
        p = NvdParams(preprocess=True, labkey=False)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        # Find the value after --preprocess
        preprocess_idx = cmd.index("--preprocess")
        assert cmd[preprocess_idx + 1] == "true"
        labkey_idx = cmd.index("--labkey")
        assert cmd[labkey_idx + 1] == "false"

    def test_path_to_string(self):
        """Paths are converted to strings."""
        p = NvdParams(samplesheet=Path("/tmp/samples.csv"))
        cmd = p.to_nextflow_args(Path("/pipeline"))
        samplesheet_idx = cmd.index("--samplesheet")
        assert cmd[samplesheet_idx + 1] == "/tmp/samples.csv"

    def test_none_values_excluded(self):
        """None values are not included in command."""
        p = NvdParams(tools=None, samplesheet=None)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--tools" not in cmd
        assert "--samplesheet" not in cmd

    def test_list_to_comma_separated(self):
        """Lists are converted to comma-separated strings."""
        p = NvdParams(human_virus_families=["Adenoviridae", "Coronaviridae"])
        cmd = p.to_nextflow_args(Path("/pipeline"))
        families_idx = cmd.index("--human_virus_families")
        assert cmd[families_idx + 1] == "Adenoviridae,Coronaviridae"


class TestNvdParamsDefaults:
    """Tests for default values matching nextflow.config."""

    def test_default_cutoff_percent(self):
        """Default cutoff_percent matches nextflow.config."""
        assert NvdParams().cutoff_percent == 0.001

    def test_default_entropy(self):
        """Default entropy matches nextflow.config."""
        assert NvdParams().entropy == 0.9

    def test_default_tax_stringency(self):
        """Default tax_stringency matches nextflow.config."""
        assert NvdParams().tax_stringency == 0.7

    def test_default_min_gottcha_reads(self):
        """Default min_gottcha_reads matches nextflow.config."""
        assert NvdParams().min_gottcha_reads == 250

    def test_default_max_blast_targets(self):
        """Default max_blast_targets matches nextflow.config."""
        assert NvdParams().max_blast_targets == 100

    def test_default_preprocess(self):
        """Default preprocess matches nextflow.config."""
        assert NvdParams().preprocess is False

    def test_default_labkey(self):
        """Default labkey matches nextflow.config."""
        assert NvdParams().labkey is False

    def test_default_human_virus_families(self):
        """Default human_virus_families matches nextflow.config."""
        p = NvdParams()
        # Field creates a mutable list from the immutable tuple constant
        assert p.human_virus_families == list(DEFAULT_HUMAN_VIRUS_FAMILIES)
        # Verify it's a copy, not the same object
        assert p.human_virus_families is not DEFAULT_HUMAN_VIRUS_FAMILIES

    def test_default_hostile_index_name(self):
        """Default hostile_index_name matches nextflow.config."""
        assert (
            NvdParams().hostile_index_name
            == "human-t2t-hla.rs-viral-202401_ml-phage-202401"
        )

    def test_default_monoimage(self):
        """Default monoimage matches nextflow.config."""
        assert NvdParams().monoimage == "nrminor/nvd2:v0.1.0"


class TestLoadParamsFile:
    """Tests for load_params_file() from params module."""

    def test_load_yaml_file(self, tmp_path):
        """Can load params from YAML file."""
        yaml_file = tmp_path / "params.yaml"
        yaml_file.write_text("tools: blast\ncutoff_percent: 0.01\n")

        params = load_params_file(yaml_file)
        assert params["tools"] == "blast"
        assert params["cutoff_percent"] == 0.01

    def test_load_yml_extension(self, tmp_path):
        """Can load params from .yml file."""
        yml_file = tmp_path / "params.yml"
        yml_file.write_text("tools: gottcha\n")

        params = load_params_file(yml_file)
        assert params["tools"] == "gottcha"

    def test_load_json_file(self, tmp_path):
        """Can load params from JSON file."""
        json_file = tmp_path / "params.json"
        json_file.write_text('{"tools": "blast", "cutoff_percent": 0.01}')

        params = load_params_file(json_file)
        assert params["tools"] == "blast"
        assert params["cutoff_percent"] == 0.01

    def test_strips_schema_key(self, tmp_path):
        """$schema key is stripped from loaded params."""
        yaml_file = tmp_path / "params.yaml"
        yaml_file.write_text(
            '$schema: "https://example.com/schema.json"\ntools: blast\n'
        )

        params = load_params_file(yaml_file)
        assert "$schema" not in params
        assert params["tools"] == "blast"

    def test_empty_file_returns_empty_dict(self, tmp_path):
        """Empty file returns empty dict."""
        yaml_file = tmp_path / "empty.yaml"
        yaml_file.write_text("")

        params = load_params_file(yaml_file)
        assert params == {}

    def test_file_not_found_raises(self, tmp_path):
        """Missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_params_file(tmp_path / "nonexistent.yaml")

    def test_integration_with_nvdparams_merge(self, tmp_path):
        """load_params_file output works with NvdParams.merge()."""
        yaml_file = tmp_path / "params.yaml"
        yaml_file.write_text("tools: blast\ncutoff_percent: 0.05\n")

        file_params = load_params_file(yaml_file)
        cli_params = {"cutoff_percent": 0.01}

        merged = NvdParams.merge(file_params, cli_params)
        assert merged.tools == "blast"  # From file
        assert merged.cutoff_percent == 0.01  # CLI wins


class TestTraceMerge:
    """Tests for trace_merge() function."""

    def test_basic_tracing(self):
        """trace_merge returns TracedParams with source info."""
        traced = trace_merge(
            ("preset", {"tools": "all", "cutoff_percent": 0.01}),
            ("file", {"tools": "blast"}),
        )

        assert isinstance(traced, TracedParams)
        assert traced.params.tools == "blast"
        assert traced.params.cutoff_percent == 0.01

    def test_tracks_override(self):
        """Overridden values are tracked."""
        traced = trace_merge(
            ("preset", {"tools": "all"}),
            ("file", {"tools": "blast"}),
        )

        # Find the tools source
        tools_source = next(s for s in traced.sources if s.field == "tools")
        assert tools_source.value == "blast"
        assert tools_source.source == "file"
        assert tools_source.overridden_value == "all"
        assert tools_source.overridden_source == "preset"

    def test_tracks_non_overridden(self):
        """Non-overridden values show their source without override info."""
        traced = trace_merge(
            ("preset", {"cutoff_percent": 0.01}),
            ("file", {"tools": "blast"}),
        )

        # cutoff_percent came from preset, not overridden
        cutoff_source = next(s for s in traced.sources if s.field == "cutoff_percent")
        assert cutoff_source.value == 0.01
        assert cutoff_source.source == "preset"
        assert cutoff_source.overridden_value is None
        assert cutoff_source.overridden_source is None

    def test_tracks_defaults(self):
        """Default values are tracked as source='default'."""
        traced = trace_merge(
            ("file", {"tools": "blast"}),
        )

        # entropy should be default
        entropy_source = next(s for s in traced.sources if s.field == "entropy")
        assert entropy_source.value == 0.9
        assert entropy_source.source == "default"

    def test_counts_defaults(self):
        """defaults_used counts fields using default values."""
        traced = trace_merge(
            ("file", {"tools": "blast", "cutoff_percent": 0.01}),
        )

        # Should have many defaults (all the fields we didn't set)
        assert traced.defaults_used > 10

    def test_empty_sources(self):
        """trace_merge with no sources returns all defaults."""
        traced = trace_merge()

        assert traced.params.tools is None  # No default for tools
        assert traced.params.cutoff_percent == 0.001  # Default
        assert traced.defaults_used > 0

    def test_none_sources_skipped(self):
        """None sources are gracefully skipped."""
        traced = trace_merge(
            ("preset", None),
            ("file", {"tools": "blast"}),
        )

        assert traced.params.tools == "blast"
        tools_source = next(s for s in traced.sources if s.field == "tools")
        assert tools_source.source == "file"

    def test_three_way_override(self):
        """Three-way merge tracks the immediate predecessor."""
        traced = trace_merge(
            ("preset", {"tools": "all"}),
            ("file", {"tools": "gottcha"}),
            ("cli", {"tools": "blast"}),
        )

        tools_source = next(s for s in traced.sources if s.field == "tools")
        assert tools_source.value == "blast"
        assert tools_source.source == "cli"
        # Should show it overrode "file", not "preset"
        assert tools_source.overridden_value == "gottcha"
        assert tools_source.overridden_source == "file"

    def test_with_nvdparams_instance(self):
        """trace_merge accepts NvdParams instances as sources."""
        preset = NvdParams(tools="all", cutoff_percent=0.01)
        traced = trace_merge(
            ("preset", preset),
            ("file", {"tools": "blast"}),
        )

        assert traced.params.tools == "blast"
        assert traced.params.cutoff_percent == 0.01

    def test_validates_merged_result(self):
        """trace_merge validates the merged params."""
        with pytest.raises(Exception):  # ValidationError
            trace_merge(
                ("file", {"tools": "invalid_tool"}),
            )
