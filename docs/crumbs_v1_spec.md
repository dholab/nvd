# CRUMBS v1 Specification

CRUMBS, short for Coverage-Robust Unified Mapped Base Support, is NVD's experimental abundance-proxy statistic for classified contigs. It summarizes filtered read mapback support while accounting for two important patterns in NVD viral contigs: gaps in mapback coverage on contigs, and extreme localized pileups. Both of these patterns can confound existing coverage statistics, making contigs contribute to much or too little to abundance profiles.

However, please note: while CRUMBS is intended to be used as a relative abundance proxy, it should not be interpreted as a genome-copy abundance. Its purpose is to associate contigs with a measure of abundance mass that can be useful in downstream benchmarks like OPAL.

## Inputs

For each classified contig, CRUMBS v1 requires a full-contig per-base depth vector from filtered read mapback:

```text
d = [d1, d2, ..., dn]
```

where `n` is the contig length and every `di` is the depth at one contig base. Positions with zero mapped reads must be included. In NVD, this vector is currently produced by streaming `samtools depth -aa` over the filtered mapback BAM.

Each contig must also have one taxonomic assignment from the final BLAST/LCA reporting path. Repeated BLAST hit rows are not independent abundance evidence; they collapse to one assignment per contig before taxon-level aggregation.

## Percentiles and caps

CRUMBS v1 uses full-contig empirical nearest-rank percentiles over the complete depth vector, including zero-depth bases. After sorting depths as:

```text
x1 <= x2 <= ... <= xn
```

the percentile cap for percentile `p` is:

```text
cap_p = xceil(p * n)
```

The p95 cap uses `p = 0.95`; the p99 cap uses `p = 0.99`. These caps are observed integer depths, not interpolated quantiles. Medians reported alongside CRUMBS use the conventional median definition rather than the nearest-rank convention.

Including zeros is part of the statistic. A contig with very low breadth can have a p95 or p99 cap of zero, in which case it contributes zero CRUMBS mass even if a small covered island exists elsewhere on the contig. Will invisible to CRUMBS, that covered island will remain visible in evidence fields such as covered bases, breadth, raw aligned bases, and maximum depth.

## Contig CRUMBS

For a contig with depth vector `d` and percentile cap `cap_p`, the percentile-specific Contig CRUMBS score is:

```text
CRUMBS_p(contig) = Σ min(di, cap_p)
```

Equivalently:

```text
CRUMBS_p(contig) = contig_length × winsorized_full_length_mean_depth_p
```

NVD emits both p95- and p99-capped values as `crumbs_p95` and `crumbs_p99`. The default CRUMBS v1 score used for taxonomic profiles is `crumbs_p99`.

The p99 default is intentionally less severe than p95 for sparse but nonzero evidence. A p95 cap becomes zero when nonzero coverage falls below about five percent of the contig; a p99 cap has the same behavior only when nonzero coverage falls below about one percent. Both caps dampen localized extreme pileups, but p99 is less likely to erase low-breadth evidence for big contigs.

## Taxon CRUMBS

For each assigned taxon, Taxon CRUMBS is the length-normalized aggregate of the default Contig CRUMBS scores for contigs assigned to that taxon:

```text
Taxon CRUMBS(taxon) = Σ CRUMBS_p99(contig) / Σ contig_length
```

The numerator and denominator include contigs assigned to that taxon by NVD's final BLAST/LCA or dominant-hit annotation. The denominator is the total assigned contig span for that taxon, not a reference genome length.

This length normalization is deliberate. Normalizing a summed Contig CRUMBS mass directly would reward a taxon for being recovered as more same-depth fragments or for having more assembled contig span. Taxon CRUMBS instead asks: among the contig bases NVD assigned to this taxon, what was the robust capped mapback support per base?

If a contig is assigned to a higher-rank LCA, its CRUMBS mass stops at that assigned taxon and its ancestors as represented in the profile path. CRUMBS v1 does not push that mass down into genus, species, or strain descendants.

## Relative profile percentage

When emitting a relative taxonomic profile, NVD normalizes positive Taxon CRUMBS values across emitted taxa:

```text
PERCENTAGE(taxon) = 100 × Taxon CRUMBS(taxon) / Σ Taxon CRUMBS(all emitted taxa)
```

Taxa with classified contigs but zero Taxon CRUMBS remain in sidecar evidence tables. They do not receive epsilon mass, and they are omitted from the BioBoxes-style profile because zero-percentage rows do not contribute abundance mass and are typically dropped by OPAL internally.

If all classified taxa have zero Taxon CRUMBS, the BioBoxes-style profile has no positive data rows. The contig-level and taxon-level evidence tables are still the source of truth for diagnosing why no abundance mass was emitted.

## Interpretation

The claim made by CRUMBS v1 is:

```text
normalized robust mapback support per assigned contig base
```

It is not:

```text
NVD genome-copy abundance
```

CRUMBS v1 exposes, rather than corrects away, NVD's assembly and classification behavior. Fragmentation, low-breadth coverage, ambiguous BLAST placement, accessory-region hits, enrichment bias, and localized pileups can all affect the profile. Those effects should be interpreted alongside the sidecar contig and taxon evidence tables.

The statistic is intentionally simple. It does not cluster contigs into genome bins, infer genome-copy abundance, de-fragment assemblies, divide by reference genome lengths, or model taxon-specific recovery bias. Those may be useful future analyses, but they are outside CRUMBS v1.
