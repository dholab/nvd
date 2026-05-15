"""Tests for py_nvd.models module, specifically NvdParams."""

import logging
from pathlib import Path

import pytest
from pydantic import ValidationError

from py_nvd.models import (
    DEFAULT_HUMAN_VIRUS_FAMILIES,
    NvdParams,
    Run,
    TracedParams,
    trace_merge,
)
from py_nvd.params import load_params_file


class TestNvdParamsInstantiation:
    """Tests for basic NvdParams instantiation."""

    def test_minimal_instantiation(self) -> None:
        """Can create NvdParams with no arguments (all defaults)."""
        p = NvdParams()
        assert p.preprocess is False
        assert p.cutoff_percent == 0.001

    def test_with_required_fields(self) -> None:
        """Can create with typical required fields."""
        p = NvdParams(
            samplesheet=Path("/fake/samples.csv"),
            experiment_id="exp001",
        )
        assert p.samplesheet == Path("/fake/samples.csv")
        assert p.experiment_id == "exp001"

    def test_all_fields(self) -> None:
        """Can create with all fields specified."""
        p = NvdParams(
            # Core
            samplesheet=Path("/fake/samples.csv"),
            results=Path("/fake/results"),
            experiment_id="exp001",
            max_concurrent_downloads=5,
            cleanup=True,
            work_dir=Path("/fake/work"),
            # Database versions
            blast_db_version="core-nt_2025-01-01",
            # Database paths
            blast_db=Path("/db/blast"),
            blast_db_prefix="nt",
            # Preprocessing
            preprocess=True,
            dedup=True,
            dedup_seq=True,
            dedup_pos=True,
            trim_adapters=True,
            host_index=Path("/db/host.idx"),
            # Analysis
            cutoff_percent=0.01,
            entropy=0.85,
            tax_stringency=0.8,
            # LabKey
            labkey=True,
            labkey_server="example.com",
        )
        assert p.cleanup is True
        assert p.host_index == Path("/db/host.idx")
        assert p.cutoff_percent == 0.01
        assert p.labkey is True


class TestNvdParamsRangeValidators:
    """Tests for 0-1 range validators."""

    def test_cutoff_percent_valid_range(self) -> None:
        """cutoff_percent accepts 0-1 range."""
        assert NvdParams(cutoff_percent=0.0).cutoff_percent == 0.0
        assert NvdParams(cutoff_percent=0.5).cutoff_percent == 0.5
        assert NvdParams(cutoff_percent=1.0).cutoff_percent == 1.0

    def test_cutoff_percent_above_one_rejected(self) -> None:
        """cutoff_percent > 1 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(cutoff_percent=1.5)
        assert "Must be between 0 and 1" in str(exc_info.value)

    def test_cutoff_percent_below_zero_rejected(self) -> None:
        """cutoff_percent < 0 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(cutoff_percent=-0.1)
        assert "Must be between 0 and 1" in str(exc_info.value)

    def test_entropy_valid_range(self) -> None:
        """entropy accepts 0-1 range."""
        assert NvdParams(entropy=0.0).entropy == 0.0
        assert NvdParams(entropy=0.9).entropy == 0.9
        assert NvdParams(entropy=1.0).entropy == 1.0

    def test_entropy_out_of_range_rejected(self) -> None:
        """entropy outside 0-1 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(entropy=1.1)
        with pytest.raises(ValidationError):
            NvdParams(entropy=-0.5)

    def test_tax_stringency_valid_range(self) -> None:
        """tax_stringency accepts 0-1 range."""
        assert NvdParams(tax_stringency=0.0).tax_stringency == 0.0
        assert NvdParams(tax_stringency=0.7).tax_stringency == 0.7
        assert NvdParams(tax_stringency=1.0).tax_stringency == 1.0

    def test_tax_stringency_out_of_range_rejected(self) -> None:
        """tax_stringency outside 0-1 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(tax_stringency=2.0)


class TestNvdParamsPositiveIntValidators:
    """Tests for positive integer validators."""

    def test_host_index_build_params_valid(self) -> None:
        """Host index build parameters accept positive integers."""
        p = NvdParams(host_kmer_size=31, host_window_size=15, host_abs_threshold=2)
        assert p.host_kmer_size == 31
        assert p.host_window_size == 15
        assert p.host_abs_threshold == 2

    def test_virus_index_build_params_valid(self) -> None:
        """Virus index build parameters accept positive integers."""
        p = NvdParams(virus_kmer_size=31, virus_window_size=1, virus_abs_threshold=1)
        assert p.virus_kmer_size == 31
        assert p.virus_window_size == 1
        assert p.virus_abs_threshold == 1

    def test_host_kmer_size_zero_rejected(self) -> None:
        """host_kmer_size=0 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(host_kmer_size=0)
        assert "Must be >= 1" in str(exc_info.value)

    def test_host_window_size_negative_rejected(self) -> None:
        """Negative host_window_size raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(host_window_size=-5)

    def test_virus_window_size_zero_rejected(self) -> None:
        """virus_window_size=0 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(virus_window_size=0)

    def test_max_blast_targets_valid(self) -> None:
        """max_blast_targets accepts positive integers."""
        assert NvdParams(max_blast_targets=1).max_blast_targets == 1
        assert NvdParams(max_blast_targets=100).max_blast_targets == 100

    def test_max_blast_targets_zero_rejected(self) -> None:
        """max_blast_targets=0 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(max_blast_targets=0)

    def test_min_read_length_valid(self) -> None:
        """min_read_length accepts positive integers."""
        assert NvdParams(min_read_length=1).min_read_length == 1
        assert NvdParams(min_read_length=50).min_read_length == 50

    def test_min_read_length_zero_rejected(self) -> None:
        """min_read_length=0 raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(min_read_length=0)


class TestNvdParamsNonNegativeIntValidators:
    """Tests for non-negative integer validators."""

    def test_min_read_quality_illumina_zero_allowed(self) -> None:
        """min_read_quality_illumina=0 is valid."""
        p = NvdParams(min_read_quality_illumina=0)
        assert p.min_read_quality_illumina == 0

    def test_min_read_quality_illumina_positive_allowed(self) -> None:
        """Positive min_read_quality_illumina is valid."""
        p = NvdParams(min_read_quality_illumina=20)
        assert p.min_read_quality_illumina == 20

    def test_min_read_quality_illumina_negative_rejected(self) -> None:
        """Negative min_read_quality_illumina raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(min_read_quality_illumina=-1)
        assert "Must be >= 0" in str(exc_info.value)

    def test_min_read_quality_nanopore_zero_allowed(self) -> None:
        """min_read_quality_nanopore=0 is valid."""
        p = NvdParams(min_read_quality_nanopore=0)
        assert p.min_read_quality_nanopore == 0

    def test_min_read_quality_nanopore_negative_rejected(self) -> None:
        """Negative min_read_quality_nanopore raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(min_read_quality_nanopore=-5)


class TestNvdParamsMaxReadLength:
    """Tests for max_read_length validator (optional positive int)."""

    def test_none_allowed(self) -> None:
        """max_read_length=None is valid (no limit)."""
        p = NvdParams(max_read_length=None)
        assert p.max_read_length is None

    def test_positive_allowed(self) -> None:
        """Positive max_read_length is valid."""
        p = NvdParams(max_read_length=500)
        assert p.max_read_length == 500

    def test_zero_rejected(self) -> None:
        """max_read_length=0 raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(max_read_length=0)
        assert "Must be >= 1" in str(exc_info.value)

    def test_negative_rejected(self) -> None:
        """Negative max_read_length raises ValidationError."""
        with pytest.raises(ValidationError):
            NvdParams(max_read_length=-100)


class TestNvdParamsTypeCoercion:
    """Tests for Pydantic type coercion."""

    def test_string_to_path(self) -> None:
        """String paths are coerced to Path objects."""
        p = NvdParams(samplesheet="/fake/samples.csv")  # ty: ignore[invalid-argument-type]
        assert isinstance(p.samplesheet, Path)
        assert p.samplesheet == Path("/fake/samples.csv")

    def test_string_to_int(self) -> None:
        """String integers are coerced to int."""
        p = NvdParams(max_blast_targets="100")  # ty: ignore[invalid-argument-type]
        assert isinstance(p.max_blast_targets, int)
        assert p.max_blast_targets == 100

    def test_string_to_float(self) -> None:
        """String floats are coerced to float."""
        p = NvdParams(cutoff_percent="0.05")  # ty: ignore[invalid-argument-type]
        assert isinstance(p.cutoff_percent, float)
        assert p.cutoff_percent == 0.05

    def test_int_to_float(self) -> None:
        """Integers are coerced to float where expected."""
        p = NvdParams(entropy=1)
        assert isinstance(p.entropy, float)
        assert p.entropy == 1.0


class TestNvdParamsMerge:
    """Tests for NvdParams.merge() class method."""

    def test_merge_empty_sources(self) -> None:
        """Merge with no sources returns defaults."""
        p = NvdParams.merge()
        assert p.cutoff_percent == 0.001

    def test_merge_single_dict(self) -> None:
        """Merge with single dict source."""
        p = NvdParams.merge({"dedup": True, "cutoff_percent": 0.01})
        assert p.dedup is True
        assert p.cutoff_percent == 0.01

    def test_merge_precedence_later_wins(self) -> None:
        """Later sources override earlier sources."""
        preset = {"dedup": False, "cutoff_percent": 0.01}
        cli = {"dedup": True}
        p = NvdParams.merge(preset, cli)
        assert p.dedup is True  # CLI wins
        assert p.cutoff_percent == 0.01  # From preset

    def test_merge_three_sources(self) -> None:
        """Three-way merge with correct precedence."""
        preset = {"dedup": True, "cutoff_percent": 0.01, "entropy": 0.8}
        params_file = {"dedup_seq": True, "entropy": 0.85}
        cli = {"cutoff_percent": 0.005}
        p = NvdParams.merge(preset, params_file, cli)
        assert p.dedup_seq is True  # From params_file
        assert p.cutoff_percent == 0.005  # From CLI
        assert p.entropy == 0.85  # From params_file

    def test_merge_none_sources_skipped(self) -> None:
        """None sources are gracefully skipped."""
        p = NvdParams.merge(None, {"dedup": True}, None)
        assert p.dedup is True

    def test_merge_none_values_skipped(self) -> None:
        """None values within dicts don't override."""
        preset = {"dedup": True, "cutoff_percent": 0.01}
        cli = {"dedup": None, "cutoff_percent": 0.005}
        p = NvdParams.merge(preset, cli)
        assert p.dedup is True  # None didn't override
        assert p.cutoff_percent == 0.005  # Non-None did override

    def test_merge_with_nvdparams_instance(self) -> None:
        """Can merge NvdParams instances, not just dicts."""
        base = NvdParams(dedup=True, cutoff_percent=0.01)
        override = {"dedup_pos": True}
        p = NvdParams.merge(base, override)
        assert p.dedup is True
        assert p.dedup_pos is True
        assert p.cutoff_percent == 0.01

    def test_merge_validates_result(self) -> None:
        """Merged result is validated."""
        with pytest.raises(ValidationError):
            NvdParams.merge({"cutoff_percent": -0.1})

    def test_merge_stateless_mode_precedence(self) -> None:
        """CLI stateless/taxonomy_dir override preset values."""
        preset = {
            "stateless": False,
            "taxonomy_dir": "/preset/taxonomy",
        }
        cli = {
            "stateless": True,
            "taxonomy_dir": "/cli/taxonomy",
        }
        p = NvdParams.merge(preset, cli)
        assert p.stateless is True  # CLI wins
        assert p.taxonomy_dir == Path("/cli/taxonomy")  # CLI wins

    def test_merge_stateless_cli_none_preserves_preset(self) -> None:
        """CLI None values don't override preset stateless settings."""
        preset = {
            "stateless": True,
            "taxonomy_dir": "/preset/taxonomy",
        }
        cli = {
            "stateless": None,  # User didn't specify
            "taxonomy_dir": None,  # User didn't specify
        }
        p = NvdParams.merge(preset, cli)
        assert p.stateless is True  # Preset preserved
        assert p.taxonomy_dir == Path("/preset/taxonomy")  # Preset preserved


class TestNvdParamsToNextflowArgs:
    """Tests for NvdParams.to_nextflow_args() method."""

    def test_basic_command_structure(self) -> None:
        """Command starts with nextflow run <pipeline>."""
        p = NvdParams()
        cmd = p.to_nextflow_args(Path("/path/to/pipeline"))
        assert cmd[0] == "nextflow"
        assert cmd[1] == "run"
        assert cmd[2] == "/path/to/pipeline"

    def test_params_as_double_dash_args(self) -> None:
        """Params are formatted as --param_name value (underscores for Nextflow)."""
        p = NvdParams(dedup=True, cutoff_percent=0.01)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--dedup" in cmd
        assert "true" in cmd
        assert "--cutoff_percent" in cmd
        assert "0.01" in cmd

    def test_underscores_preserved_for_nextflow(self) -> None:
        """Param names keep underscores (Nextflow/Groovy requires them)."""
        p = NvdParams(host_kmer_size=31)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--host_kmer_size" in cmd
        assert "--host-kmer-size" not in cmd

    def test_bool_to_string(self) -> None:
        """Booleans are converted to 'true'/'false' strings."""
        p = NvdParams(preprocess=True, labkey=False)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        # Find the value after --preprocess
        preprocess_idx = cmd.index("--preprocess")
        assert cmd[preprocess_idx + 1] == "true"
        labkey_idx = cmd.index("--labkey")
        assert cmd[labkey_idx + 1] == "false"

    def test_path_to_string(self) -> None:
        """Paths are converted to strings."""
        p = NvdParams(samplesheet=Path("/fake/samples.csv"))
        cmd = p.to_nextflow_args(Path("/pipeline"))
        samplesheet_idx = cmd.index("--samplesheet")
        assert cmd[samplesheet_idx + 1] == "/fake/samples.csv"

    def test_none_values_excluded(self) -> None:
        """None values are not included in command."""
        p = NvdParams(samplesheet=None, host_index=None)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--samplesheet" not in cmd
        assert "--host_index" not in cmd

    def test_list_to_comma_separated(self) -> None:
        """Lists are converted to comma-separated strings."""
        p = NvdParams(human_virus_families=["Adenoviridae", "Coronaviridae"])
        cmd = p.to_nextflow_args(Path("/pipeline"))
        families_idx = cmd.index("--human_virus_families")
        assert cmd[families_idx + 1] == "Adenoviridae,Coronaviridae"

    def test_stateless_mode_params(self) -> None:
        """Stateless mode params are correctly propagated."""
        p = NvdParams(stateless=True, taxonomy_dir=Path("/shared/taxonomy"))
        cmd = p.to_nextflow_args(Path("/pipeline"))

        # stateless should be "true"
        assert "--stateless" in cmd
        stateless_idx = cmd.index("--stateless")
        assert cmd[stateless_idx + 1] == "true"

        # taxonomy_dir should be the path string
        assert "--taxonomy_dir" in cmd
        taxonomy_idx = cmd.index("--taxonomy_dir")
        assert cmd[taxonomy_idx + 1] == "/shared/taxonomy"

    def test_stateless_false_still_included(self) -> None:
        """stateless=False is included (explicit false, not omitted)."""
        p = NvdParams(stateless=False)
        cmd = p.to_nextflow_args(Path("/pipeline"))

        # stateless should be "false" (not omitted)
        assert "--stateless" in cmd
        stateless_idx = cmd.index("--stateless")
        assert cmd[stateless_idx + 1] == "false"

    def test_taxonomy_dir_none_excluded(self) -> None:
        """taxonomy_dir=None is excluded from command."""
        p = NvdParams(taxonomy_dir=None)
        cmd = p.to_nextflow_args(Path("/pipeline"))
        assert "--taxonomy_dir" not in cmd


class TestNvdParamsDefaults:
    """Tests for default values matching nextflow.config."""

    def test_default_cutoff_percent(self) -> None:
        """Default cutoff_percent matches nextflow.config."""
        assert NvdParams().cutoff_percent == 0.001

    def test_default_entropy(self) -> None:
        """Default entropy matches nextflow.config."""
        assert NvdParams().entropy == 0.9

    def test_default_tax_stringency(self) -> None:
        """Default tax_stringency matches nextflow.config."""
        assert NvdParams().tax_stringency == 0.7

    def test_default_host_index_url(self) -> None:
        """Host depletion is off by default."""
        assert NvdParams().host_index is None
        assert NvdParams().host_index_url is None
        assert NvdParams().host_contaminants_fasta is None

    def test_default_virus_index_sources(self) -> None:
        """Virus enrichment has no index source by default."""
        assert NvdParams().virus_index is None
        assert NvdParams().virus_index_url is None
        assert NvdParams().virus_reference_fasta is None
        assert NvdParams().virus_kmer_size == 31
        assert NvdParams().virus_window_size == 1
        assert NvdParams().virus_abs_threshold == 1
        assert NvdParams().virus_rel_threshold == 0.0

    def test_default_max_blast_targets(self) -> None:
        """Default max_blast_targets matches nextflow.config."""
        assert NvdParams().max_blast_targets == 100

    def test_default_preprocess(self) -> None:
        """Default preprocess matches nextflow.config."""
        assert NvdParams().preprocess is False

    def test_default_labkey(self) -> None:
        """Default labkey matches nextflow.config."""
        assert NvdParams().labkey is False

    def test_default_human_virus_families(self) -> None:
        """Default human_virus_families matches nextflow.config."""
        p = NvdParams()
        # Field creates a mutable list from the immutable tuple constant
        assert p.human_virus_families == list(DEFAULT_HUMAN_VIRUS_FAMILIES)
        # Verify it's a copy, not the same object
        assert p.human_virus_families is not DEFAULT_HUMAN_VIRUS_FAMILIES

    def test_default_monoimage(self) -> None:
        """Default monoimage matches nextflow.config."""
        assert NvdParams().monoimage == "nrminor/nvd:v3.0.0"


class TestLoadParamsFile:
    """Tests for load_params_file() from params module."""

    def test_load_yaml_file(self, tmp_path: Path) -> None:
        """Can load params from YAML file."""
        yaml_file = tmp_path / "params.yaml"
        yaml_file.write_text("dedup: true\ncutoff_percent: 0.01\n")

        params = load_params_file(yaml_file)
        assert params["dedup"] is True
        assert params["cutoff_percent"] == 0.01

    def test_load_yml_extension(self, tmp_path: Path) -> None:
        """Can load params from .yml file."""
        yml_file = tmp_path / "params.yml"
        yml_file.write_text("dedup_seq: true\n")

        params = load_params_file(yml_file)
        assert params["dedup_seq"] is True

    def test_load_json_file(self, tmp_path: Path) -> None:
        """Can load params from JSON file."""
        json_file = tmp_path / "params.json"
        json_file.write_text('{"dedup_pos": true, "cutoff_percent": 0.01}')

        params = load_params_file(json_file)
        assert params["dedup_pos"] is True
        assert params["cutoff_percent"] == 0.01

    def test_strips_schema_key(self, tmp_path: Path) -> None:
        """$schema key is stripped from loaded params."""
        yaml_file = tmp_path / "params.yaml"
        yaml_file.write_text(
            '$schema: "https://example.com/schema.json"\ndedup: true\n',
        )

        params = load_params_file(yaml_file)
        assert "$schema" not in params
        assert params["dedup"] is True

    def test_empty_file_returns_empty_dict(self, tmp_path: Path) -> None:
        """Empty file returns empty dict."""
        yaml_file = tmp_path / "empty.yaml"
        yaml_file.write_text("")

        params = load_params_file(yaml_file)
        assert params == {}

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        """Missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_params_file(tmp_path / "nonexistent.yaml")

    def test_integration_with_nvdparams_merge(self, tmp_path: Path) -> None:
        """load_params_file output works with NvdParams.merge()."""
        yaml_file = tmp_path / "params.yaml"
        yaml_file.write_text("dedup: true\ncutoff_percent: 0.05\n")

        file_params = load_params_file(yaml_file)
        cli_params = {"cutoff_percent": 0.01}

        merged = NvdParams.merge(file_params, cli_params)
        assert merged.dedup is True  # From file
        assert merged.cutoff_percent == 0.01  # CLI wins


class TestTraceMerge:
    """Tests for trace_merge() function."""

    def test_basic_tracing(self) -> None:
        """trace_merge returns TracedParams with source info."""
        traced = trace_merge(
            ("preset", {"dedup_seq": True, "cutoff_percent": 0.01}),
            ("file", {"dedup": True}),
        )

        assert isinstance(traced, TracedParams)
        assert traced.params.dedup is True
        assert traced.params.dedup_seq is True
        assert traced.params.cutoff_percent == 0.01

    def test_tracks_override(self) -> None:
        """Overridden values are tracked."""
        traced = trace_merge(
            ("preset", {"cutoff_percent": 0.01}),
            ("file", {"cutoff_percent": 0.005}),
        )

        cutoff_source = next(s for s in traced.sources if s.field == "cutoff_percent")
        assert cutoff_source.value == 0.005
        assert cutoff_source.source == "file"
        assert cutoff_source.overridden_value == 0.01
        assert cutoff_source.overridden_source == "preset"

    def test_tracks_non_overridden(self) -> None:
        """Non-overridden values show their source without override info."""
        traced = trace_merge(
            ("preset", {"cutoff_percent": 0.01}),
            ("file", {"dedup": True}),
        )

        # cutoff_percent came from preset, not overridden
        cutoff_source = next(s for s in traced.sources if s.field == "cutoff_percent")
        assert cutoff_source.value == 0.01
        assert cutoff_source.source == "preset"
        assert cutoff_source.overridden_value is None
        assert cutoff_source.overridden_source is None

    def test_tracks_defaults(self) -> None:
        """Default values are tracked as source='default'."""
        traced = trace_merge(
            ("file", {"dedup": True}),
        )

        # entropy should be default
        entropy_source = next(s for s in traced.sources if s.field == "entropy")
        assert entropy_source.value == 0.9
        assert entropy_source.source == "default"

    def test_counts_defaults(self) -> None:
        """defaults_used counts fields using default values."""
        traced = trace_merge(
            ("file", {"dedup": True, "cutoff_percent": 0.01}),
        )

        # Should have many defaults (all the fields we didn't set)
        assert traced.defaults_used > 10

    def test_empty_sources(self) -> None:
        """trace_merge with no sources returns all defaults."""
        traced = trace_merge()

        assert traced.params.dedup is False  # Default
        assert traced.params.cutoff_percent == 0.001  # Default
        assert traced.defaults_used > 0

    def test_none_sources_skipped(self) -> None:
        """None sources are gracefully skipped."""
        traced = trace_merge(
            ("preset", None),
            ("file", {"dedup": True}),
        )

        assert traced.params.dedup is True
        dedup_source = next(s for s in traced.sources if s.field == "dedup")
        assert dedup_source.source == "file"

    def test_three_way_override(self) -> None:
        """Three-way merge tracks the immediate predecessor."""
        traced = trace_merge(
            ("preset", {"dedup_seq": True}),
            ("file", {"dedup_pos": True}),
            ("cli", {"dedup": True}),
        )

        dedup_source = next(s for s in traced.sources if s.field == "dedup")
        assert dedup_source.value is True
        assert dedup_source.source == "cli"
        # Should show it overrode "file", not "preset"
        assert dedup_source.overridden_value is None
        assert dedup_source.overridden_source is None

    def test_with_nvdparams_instance(self) -> None:
        """trace_merge accepts NvdParams instances as sources."""
        preset = NvdParams(dedup=True, cutoff_percent=0.01)
        traced = trace_merge(
            ("preset", preset),
            ("file", {"dedup_pos": True}),
        )

        assert traced.params.dedup is True
        assert traced.params.dedup_pos is True
        assert traced.params.cutoff_percent == 0.01

    def test_validates_merged_result(self) -> None:
        """trace_merge validates the merged params."""
        with pytest.raises(ValidationError):
            trace_merge(
                ("file", {"cutoff_percent": -0.1}),
            )


class TestNvdParamsSlackChannelValidator:
    """Tests for slack_channel field validation."""

    def test_valid_channel_id(self) -> None:
        """Valid Slack channel IDs are accepted."""
        p = NvdParams(slack_channel="C0123456789")
        assert p.slack_channel == "C0123456789"

    def test_valid_channel_id_short(self) -> None:
        """Short channel IDs are accepted."""
        p = NvdParams(slack_channel="C123")
        assert p.slack_channel == "C123"

    def test_valid_channel_id_long(self) -> None:
        """Long channel IDs are accepted."""
        p = NvdParams(slack_channel="C0123456789ABCDEF")
        assert p.slack_channel == "C0123456789ABCDEF"

    def test_none_allowed(self) -> None:
        """None is a valid value for slack_channel."""
        p = NvdParams(slack_channel=None)
        assert p.slack_channel is None

    def test_invalid_missing_c_prefix(self) -> None:
        """Channel ID without C prefix is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="0123456789")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_lowercase_c(self) -> None:
        """Channel ID with lowercase c prefix is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="c0123456789")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_lowercase_letters(self) -> None:
        """Channel ID with lowercase letters is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="C0123abcdef")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_channel_name(self) -> None:
        """Channel name (not ID) is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="#general")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_empty_after_c(self) -> None:
        """Channel ID with only C is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="C")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_empty_string(self) -> None:
        """Empty string is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_whitespace(self) -> None:
        """Channel ID with whitespace is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="C0123 456789")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_special_characters(self) -> None:
        """Channel ID with special characters is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="C0123-456789")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_underscore(self) -> None:
        """Channel ID with underscore is rejected."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="C0123_456789")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_invalid_hash_prefix(self) -> None:
        """Channel name with # prefix is rejected even with valid ID after."""
        with pytest.raises(ValidationError) as exc_info:
            NvdParams(slack_channel="#C0123456789")
        assert "Invalid Slack channel ID" in str(exc_info.value)

    def test_slack_enabled_default(self) -> None:
        """Default slack_enabled is False."""
        p = NvdParams()
        assert p.slack_enabled is False

    def test_slack_enabled_true(self) -> None:
        """slack_enabled can be set to True."""
        p = NvdParams(slack_enabled=True)
        assert p.slack_enabled is True


class TestRunDurationSeconds:
    """Tests for Run.duration_seconds property."""

    def _make_completed_run(self, started_at: str, completed_at: str) -> Run:
        """Helper to create a completed Run with the given timestamps."""
        return Run(
            run_id="test_run",
            sample_set_id="abc123",
            started_at=started_at,
            status="completed",
            completed_at=completed_at,
        )

    def _make_running_run(self, started_at: str) -> Run:
        """Helper to create a running Run (no completed_at)."""
        return Run(
            run_id="test_run",
            sample_set_id="abc123",
            started_at=started_at,
            status="running",
            completed_at=None,
        )

    def test_completed_run_returns_duration(self) -> None:
        """Completed run with valid timestamps returns duration in seconds."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00Z",
            completed_at="2024-01-01T10:05:00Z",
        )
        assert run.duration_seconds == 300.0

    def test_in_progress_run_returns_none(self) -> None:
        """Run without completed_at returns None."""
        run = self._make_running_run(started_at="2024-01-01T10:00:00Z")
        assert run.duration_seconds is None

    def test_zero_duration(self) -> None:
        """Same start and end time returns 0."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00Z",
            completed_at="2024-01-01T10:00:00Z",
        )
        assert run.duration_seconds == 0.0

    def test_z_suffix_timestamps(self) -> None:
        """Timestamps with Z suffix are parsed correctly."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00Z",
            completed_at="2024-01-01T11:00:00Z",
        )
        assert run.duration_seconds == 3600.0

    def test_explicit_utc_offset(self) -> None:
        """Timestamps with +00:00 offset are parsed correctly."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00+00:00",
            completed_at="2024-01-01T11:00:00+00:00",
        )
        assert run.duration_seconds == 3600.0

    def test_naive_timestamps(self) -> None:
        """Naive timestamps (no timezone) are handled correctly."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00",
            completed_at="2024-01-01T11:00:00",
        )
        assert run.duration_seconds == 3600.0

    def test_mixed_timezone_awareness(self) -> None:
        """Mixed naive and aware timestamps are handled correctly."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00Z",
            completed_at="2024-01-01T11:00:00",  # naive
        )
        assert run.duration_seconds == 3600.0

    def test_different_timezone_offsets(self) -> None:
        """Different timezone offsets compute correct duration."""
        # 10:00 UTC and 12:00+02:00 are the same instant
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00+00:00",
            completed_at="2024-01-01T12:00:00+02:00",
        )
        assert run.duration_seconds == 0.0

    def test_timestamps_with_milliseconds(self) -> None:
        """Timestamps with milliseconds are parsed correctly."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00.000Z",
            completed_at="2024-01-01T10:00:01.500Z",
        )
        assert run.duration_seconds == 1.5

    def test_timestamps_with_microseconds(self) -> None:
        """Timestamps with microseconds are parsed correctly."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00.000000+00:00",
            completed_at="2024-01-01T10:00:00.500000+00:00",
        )
        assert run.duration_seconds == 0.5

    def test_invalid_started_at_returns_none(self) -> None:
        """Invalid started_at timestamp returns None."""
        run = self._make_completed_run(
            started_at="not-a-date",
            completed_at="2024-01-01T10:00:00Z",
        )
        assert run.duration_seconds is None

    def test_invalid_completed_at_returns_none(self) -> None:
        """Invalid completed_at timestamp returns None."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00Z",
            completed_at="not-a-date",
        )
        assert run.duration_seconds is None

    def test_empty_started_at_returns_none(self) -> None:
        """Empty started_at string returns None."""
        run = self._make_completed_run(
            started_at="",
            completed_at="2024-01-01T10:00:00Z",
        )
        assert run.duration_seconds is None

    def test_empty_completed_at_returns_none(self) -> None:
        """Empty completed_at string returns None."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:00:00Z",
            completed_at="",
        )
        assert run.duration_seconds is None

    def test_negative_duration_returns_none(self) -> None:
        """Negative duration (end before start) returns None."""
        run = self._make_completed_run(
            started_at="2024-01-01T10:05:00Z",
            completed_at="2024-01-01T10:00:00Z",  # Before start
        )
        assert run.duration_seconds is None

    def test_negative_duration_logs_warning(
        self,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Negative duration logs a warning."""
        with caplog.at_level(logging.WARNING):
            run = self._make_completed_run(
                started_at="2024-01-01T10:05:00Z",
                completed_at="2024-01-01T10:00:00Z",
            )
            _ = run.duration_seconds

        assert "negative duration" in caplog.text.lower()

    def test_invalid_timestamp_logs_warning(
        self,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Invalid timestamp logs a warning."""
        with caplog.at_level(logging.WARNING):
            run = self._make_completed_run(
                started_at="not-a-date",
                completed_at="2024-01-01T10:00:00Z",
            )
            _ = run.duration_seconds

        assert "failed to parse" in caplog.text.lower()
