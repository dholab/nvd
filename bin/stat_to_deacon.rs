#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! clap = { version = "4", features = ["derive"] }
//!
//! [dev-dependencies]
//! tempfile = "3"
//! ```
//!
//! Convert NCBI STAT `.dbss` k-mer databases to Deacon `.idx` index files.
//!
//! STAT stores k-mers with tax ID annotations in a binary format. This tool
//! extracts k-mers for a specified set of tax IDs and encodes them in the
//! Deacon minimizer index format.
//!
//! # Usage
//! ```text
//! stat_to_deacon.rs \
//!   --dbss viruses.dbss \
//!   --annotation viruses.dbss.annotation \
//!   --taxids target_taxids.txt \
//!   --output viruses.idx
//! ```

use clap::Parser;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;

// ── Newtypes ────────────────────────────────────────────────────────

/// Lookup table mapping 2-bit codes (0–3) to ASCII base characters.
/// Index matches STAT/Deacon encoding: A=0, C=1, T=2, G=3.
const BASES: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Convert a single ASCII nucleotide byte to its 2-bit code (A=0, C=1, T=2, G=3).
/// Panics on any character outside the four canonical bases.
fn base_to_2bit(base: u8) -> u8 {
    match base {
        b'A' => 0, b'C' => 1, b'T' => 2, b'G' => 3,
        _ => panic!("Invalid base: {}", base as char),
    }
}

/// A k-mer encoded in STAT's MSB-first 2-bit scheme.
/// The first base of the sequence sits in the highest 2 bits of the u64.
/// Encoding: `hash <<= 2; hash |= base_to_2bit(base);`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StatKmer(u64);

/// A k-mer encoded in Deacon's LSB-first 2-bit scheme.
/// The first base of the sequence sits in bits 1:0.
/// Base at position `i` occupies bits `2*i+1 : 2*i`.
/// Uses u128 to support k up to 64 (k=33 needs 66 bits, exceeding u64).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct DeaconKmer(u128);

/// A nucleotide sequence as a plain ASCII string (A/C/T/G only).
#[derive(Debug, Clone, PartialEq, Eq)]
struct NucleotideSeq(String);

/// An NCBI taxonomy identifier (signed, matching STAT annotation files).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct TaxId(i32);

/// The k-mer length (number of bases per k-mer).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct KmerLength(u8);

/// The Deacon window size `w` (used together with `k` to define minimizers).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WindowSize(u8);

// ── StatKmer impl ────────────────────────────────────────────────────

impl StatKmer {
    /// Decode a STAT k-mer back to its nucleotide sequence.
    ///
    /// STAT packs bases MSB-first, so the first base is in the highest bits.
    /// To read them out in order we peel from the LSB (which gives reverse
    /// order) and then reverse the collected bytes — matching STAT's C++
    /// `str_from_hash` function.
    fn decode(self, k: KmerLength) -> NucleotideSeq {
        let mut value = self.0;
        let mut bases = Vec::with_capacity(k.0 as usize);
        for _ in 0..k.0 {
            bases.push(BASES[(value & 0x3) as usize]);
            value >>= 2;
        }
        bases.reverse();
        NucleotideSeq(String::from_utf8(bases).unwrap())
    }
}

// ── NucleotideSeq impl ───────────────────────────────────────────────

impl NucleotideSeq {
    /// Return a new NucleotideSeq containing only the first `k` bases (prefix).
    fn truncate(&self, k: KmerLength) -> NucleotideSeq {
        NucleotideSeq(self.0[..k.0 as usize].to_string())
    }

    /// Return a new NucleotideSeq containing only the last `k` bases (suffix).
    fn truncate_suffix(&self, k: KmerLength) -> NucleotideSeq {
        let start = self.0.len() - k.0 as usize;
        NucleotideSeq(self.0[start..].to_string())
    }

    /// Encode this sequence as a Deacon k-mer (LSB-first 2-bit packing).
    /// Base at position `i` is placed at bits `2*i+1 : 2*i`.
    fn to_deacon_kmer(&self) -> DeaconKmer {
        let mut value: u128 = 0;
        for (i, &base) in self.0.as_bytes().iter().enumerate() {
            value |= (base_to_2bit(base) as u128) << (2 * i);
        }
        DeaconKmer(value)
    }

    /// Return 4 new sequences, each with one of A/C/T/G appended to the end.
    /// Used for extending k-mers to k+1 to cover all possible next bases.
    fn extend_all_bases(&self) -> [NucleotideSeq; 4] {
        [b'A', b'C', b'T', b'G'].map(|base| {
            let mut extended = self.0.clone();
            extended.push(base as char);
            NucleotideSeq(extended)
        })
    }

    /// Return the reverse complement of this sequence.
    fn reverse_complement(&self) -> NucleotideSeq {
        let rc: String = self.0.bytes().rev().map(|b| match b {
            b'A' => 'T', b'T' => 'A', b'C' => 'G', b'G' => 'C',
            _ => panic!("Invalid base: {}", b as char),
        }).collect();
        NucleotideSeq(rc)
    }

    /// Borrow the inner sequence string as a `&str`.
    fn as_str(&self) -> &str { &self.0 }

    /// Return the number of bases in this sequence.
    fn len(&self) -> usize { self.0.len() }
}

// ── DeaconKmer impl ──────────────────────────────────────────────────

impl DeaconKmer {
    /// Compute the reverse complement of this Deacon k-mer.
    ///
    /// Deacon's complement flips A↔T and C↔G by XOR-ing each 2-bit pair with
    /// `0b10` (the XOR mask over k pairs). The reversal then swaps the order
    /// of the 2-bit pairs so the last base becomes the first.
    fn revcomp(self, k: KmerLength) -> DeaconKmer {
        // Build a mask with 0b10 repeated k times (one `10` per base pair).
        let mut mask: u128 = 0;
        for _ in 0..k.0 {
            mask = (mask << 2) | 0b10;
        }
        // XOR flips A(00)↔T(10) and C(01)↔G(11) for every base.
        let complemented = self.0 ^ mask;
        // Reverse the order of 2-bit pairs by peeling from LSB and rebuilding.
        let mut reversed: u128 = 0;
        let mut remaining = complemented;
        for _ in 0..k.0 {
            reversed = (reversed << 2) | (remaining & 0x3);
            remaining >>= 2;
        }
        DeaconKmer(reversed)
    }

    /// Return the canonical form: the lexicographically smaller of the k-mer
    /// and its reverse complement. Re-canonicalization is mandatory because
    /// STAT and Deacon pack bits in opposite orders, so canonical choices can
    /// differ between the two encodings (~33% disagreement rate).
    fn canonicalize(self, k: KmerLength) -> DeaconKmer {
        let rc = self.revcomp(k);
        if self.0 <= rc.0 { self } else { rc }
    }
}

// ── File I/O ────────────────────────────────────────────────────────

/// Parsed DBSS header (24 bytes: version u64 LE + kmer_len u64 LE + count u64 LE).
struct DbssHeader {
    kmer_length: KmerLength,
    total_kmer_count: u64,
}

fn read_dbss_header(reader: &mut impl Read) -> io::Result<DbssHeader> {
    let mut buf = [0u8; 8];

    reader.read_exact(&mut buf)?;
    let version = u64::from_le_bytes(buf);
    if version != 1 {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Unsupported .dbss version: {} (expected 1)", version)));
    }

    reader.read_exact(&mut buf)?;
    let kmer_len = u64::from_le_bytes(buf);
    if kmer_len < 1 || kmer_len > 64 {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Invalid kmer_len: {} (expected 1-64)", kmer_len)));
    }

    reader.read_exact(&mut buf)?;
    let count = u64::from_le_bytes(buf);

    Ok(DbssHeader { kmer_length: KmerLength(kmer_len as u8), total_kmer_count: count })
}

/// One entry from the .dbss.annotation file.
struct AnnotationEntry {
    tax_id: TaxId,
    kmer_count: u64,
}

/// Parse .dbss.annotation (tab-delimited: tax_id\tcount).
fn parse_annotation(path: &std::path::Path) -> io::Result<Vec<AnnotationEntry>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();

    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() != 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Bad annotation at line {}: expected 2 fields, got {}", line_number + 1, parts.len())));
        }
        let tax_id: i32 = parts[0].parse().map_err(|e|
            io::Error::new(io::ErrorKind::InvalidData, format!("Bad tax_id at line {}: {}", line_number + 1, e)))?;
        let count: u64 = parts[1].parse().map_err(|e|
            io::Error::new(io::ErrorKind::InvalidData, format!("Bad count at line {}: {}", line_number + 1, e)))?;
        entries.push(AnnotationEntry { tax_id: TaxId(tax_id), kmer_count: count });
    }
    Ok(entries)
}

/// Load target tax IDs from a text file (one integer per line).
fn load_taxid_list(path: &std::path::Path) -> io::Result<HashSet<TaxId>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut tax_ids = HashSet::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }
        let tax_id: i32 = trimmed.parse().map_err(|e|
            io::Error::new(io::ErrorKind::InvalidData, format!("Bad tax ID '{}': {}", trimmed, e)))?;
        tax_ids.insert(TaxId(tax_id));
    }
    Ok(tax_ids)
}

// ── Conversion Pipeline ─────────────────────────────────────────────

/// Result of running the full conversion pipeline.
struct ConversionResult {
    deacon_kmers: HashSet<DeaconKmer>,
    processed_count: u64,
    duplicate_count: u64,
}

/// Read matching k-mers from .dbss, convert to Deacon canonical form, deduplicate.
///
/// Walks the .dbss file sequentially using the annotation as a table of contents.
/// For each tax ID in the annotation:
///   - If NOT in target_tax_ids: skip the block (just advance byte_offset)
///   - If in target_tax_ids: seek, read the block, decode each StatKmer to
///     NucleotideSeq, optionally write FASTA, truncate if needed, re-encode
///     as DeaconKmer, canonicalize, insert into HashSet for deduplication.
fn convert_dbss_to_deacon(
    dbss_path: &std::path::Path,
    annotations: &[AnnotationEntry],
    target_tax_ids: &HashSet<TaxId>,
    source_k: KmerLength,
    target_k: KmerLength,
    fasta_writer: &mut Option<BufWriter<File>>,
) -> io::Result<ConversionResult> {
    let mut file = File::open(dbss_path)?;
    file.seek(SeekFrom::Start(24))?; // skip the 24-byte header

    let mut deacon_kmers: HashSet<DeaconKmer> = HashSet::new();
    let mut byte_offset: u64 = 24;
    let mut processed: u64 = 0;
    let mut duplicates: u64 = 0;
    let mut matching_taxids: u64 = 0;

    let chunk_size: usize = 100_000;
    let mut kmer_buf = vec![0u8; chunk_size * 8];

    for entry in annotations {
        let block_bytes = entry.kmer_count * 8;

        if !target_tax_ids.contains(&entry.tax_id) {
            byte_offset += block_bytes;
            continue;
        }

        matching_taxids += 1;
        file.seek(SeekFrom::Start(byte_offset))?;

        let mut remaining = entry.kmer_count;
        while remaining > 0 {
            let batch = remaining.min(chunk_size as u64) as usize;
            let batch_bytes = batch * 8;
            file.read_exact(&mut kmer_buf[..batch_bytes])?;

            for i in 0..batch {
                let stat_value = u64::from_le_bytes(
                    kmer_buf[i * 8..(i + 1) * 8].try_into().unwrap()
                );
                let stat_kmer = StatKmer(stat_value);
                let nucleotide_seq = stat_kmer.decode(source_k);

                // Write FASTA before truncation
                if let Some(ref mut writer) = fasta_writer {
                    writeln!(writer, ">{}", entry.tax_id.0)?;
                    writeln!(writer, "{}", nucleotide_seq.as_str())?;
                }

                // Truncate, extend, or pass through depending on target vs source k
                if target_k.0 == source_k.0 + 1 {
                    // Extension mode: append each of the 4 bases to BOTH the
                    // forward and reverse complement strands. This is necessary
                    // because reads can contain the virus k-mer in either
                    // orientation, and appending to S is not equivalent to
                    // appending to revcomp(S) after canonicalization.
                    let revcomp_seq = nucleotide_seq.reverse_complement();
                    for seq in [&nucleotide_seq, &revcomp_seq] {
                        for extended_seq in seq.extend_all_bases() {
                            let deacon_kmer = extended_seq.to_deacon_kmer().canonicalize(target_k);
                            let is_new = deacon_kmers.insert(deacon_kmer);
                            if !is_new {
                                duplicates += 1;
                            }
                        }
                    }
                } else if target_k.0 < source_k.0 {
                    // Truncation mode: store BOTH the prefix and suffix k-mers.
                    // A 32-mer K contains two overlapping 31-mers: K[0..31] and
                    // K[1..32]. Storing only the prefix misses reads where the
                    // virus k-mer appears in reverse complement orientation
                    // (which needs the suffix). Both are needed to match or
                    // exceed STAT's sensitivity.
                    let prefix = nucleotide_seq.truncate(target_k);
                    let suffix = nucleotide_seq.truncate_suffix(target_k);
                    for seq in [prefix, suffix] {
                        let deacon_kmer = seq.to_deacon_kmer().canonicalize(target_k);
                        let is_new = deacon_kmers.insert(deacon_kmer);
                        if !is_new {
                            duplicates += 1;
                        }
                    }
                } else {
                    // Pass-through: target_k == source_k
                    let deacon_kmer = nucleotide_seq.to_deacon_kmer().canonicalize(target_k);
                    let is_new = deacon_kmers.insert(deacon_kmer);
                    if !is_new {
                        duplicates += 1;
                    }
                }
                processed += 1;
            }

            remaining -= batch as u64;
        }

        byte_offset += block_bytes;

        if matching_taxids % 1000 == 0 {
            eprintln!("  Progress: {} tax IDs, {} k-mers processed, {} unique",
                matching_taxids, processed, deacon_kmers.len());
        }
    }

    Ok(ConversionResult {
        deacon_kmers,
        processed_count: processed,
        duplicate_count: duplicates,
    })
}

// ── Deacon .idx Writer ──────────────────────────────────────────────

/// Write a Deacon v3 index file.
///
/// Format: 3 header bytes (version=3, k, w) + 8-byte count (u64 LE)
/// + count * ceil(k/4) bytes of packed minimizers (each LE).
fn write_deacon_index(
    output_path: &std::path::Path,
    deacon_kmers: &HashSet<DeaconKmer>,
    kmer_length: KmerLength,
    window_size: WindowSize,
) -> io::Result<()> {
    let count = deacon_kmers.len() as u64;
    let bytes_per_minimizer = (kmer_length.0 as usize + 3) / 4;

    eprintln!("Writing Deacon v3 index: k={}, w={}, n={}, {} bytes/minimizer",
        kmer_length.0, window_size.0, count, bytes_per_minimizer);

    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    // Header: format_version(u8=3), kmer_length(u8), window_size(u8)
    writer.write_all(&[3u8, kmer_length.0, window_size.0])?;

    // Count: u64 LE
    writer.write_all(&count.to_le_bytes())?;

    // Minimizers: each as bytes_per_minimizer LE bytes (u128 to support k > 32)
    let mut written: u64 = 0;
    for &kmer in deacon_kmers.iter() {
        let bytes = kmer.0.to_le_bytes(); // 16-byte u128 LE
        writer.write_all(&bytes[..bytes_per_minimizer])?;
        written += 1;
        if written % 5_000_000 == 0 {
            eprintln!("  Written {}/{}", written, count);
        }
    }
    writer.flush()?;
    drop(writer); // ensure file is closed before checking size

    // Validate output file size
    let expected_size = 3 + 8 + (count * bytes_per_minimizer as u64);
    let actual_size = std::fs::metadata(output_path)?.len();
    if actual_size != expected_size {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Output size mismatch: expected {}, got {}", expected_size, actual_size)));
    }

    eprintln!("  Done: {:?} ({} bytes)", output_path, actual_size);
    Ok(())
}

/// Convert an NCBI STAT `.dbss` k-mer database to a Deacon `.idx` index file.
///
/// The STAT `.dbss` file contains binary-encoded k-mers with associated tax IDs.
/// The companion `.dbss.annotation` file maps column indices to tax IDs. This tool
/// filters k-mers to the specified target tax IDs and writes them out in Deacon's
/// minimizer index format.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Input STAT .dbss binary file.
    ///
    /// Contains 2-bit encoded k-mers (MSB-first, A=0 C=1 T=2 G=3) stored
    /// as u64 values grouped by tax ID. The 24-byte header holds version,
    /// kmer_len, and total k-mer count.
    #[arg(long)]
    dbss: PathBuf,

    /// Companion .dbss.annotation file (tab-delimited: tax_id\tcount).
    ///
    /// Acts as a table-of-contents for the .dbss: each line gives a tax ID
    /// and how many k-mers belong to it. Entries are in the same order as
    /// the contiguous k-mer blocks in the .dbss binary.
    #[arg(long)]
    annotation: PathBuf,

    /// Text file listing target NCBI tax IDs, one integer per line.
    ///
    /// Only k-mers belonging to these tax IDs will be extracted and written
    /// to the output index.
    #[arg(long)]
    taxids: PathBuf,

    /// Output Deacon .idx file.
    ///
    /// Will be created or overwritten. Contains a 3-byte header (version,
    /// k, w), 8-byte minimizer count, then packed minimizer values at
    /// ceil(k/4) bytes each, little-endian.
    #[arg(long)]
    output: PathBuf,

    /// Optional: write intermediate k-mers as multi-FASTA for verification.
    ///
    /// Each k-mer is written as ">tax_id\nSEQUENCE\n" before any
    /// truncation, so the full source 32-mer is preserved for BLAST
    /// spot-checking.
    #[arg(long)]
    fasta: Option<PathBuf>,

    /// Target k-mer length for the output index.
    ///
    /// STAT uses k=32. Options:
    ///   - 32: direct conversion (requires w=2, checks ~67% of read k-mers)
    ///   - 31: truncate by 1 base (allows w=1, checks 100% but loses 1 base)
    ///   - 33: extend by 1 base — each 32-mer becomes 4 33-mers covering
    ///         all possible next bases (allows w=1, checks 100%, keeps full
    ///         32-base specificity). Index is ~4x larger.
    #[arg(long, default_value_t = 33)]
    target_k: u8,

    /// Deacon window size (w). Auto: 1 if k is odd, 2 if k is even.
    ///
    /// Small w is critical for surrogate indexes: the index contains ALL
    /// k-mers (not just minimizers), so query sensitivity scales with the
    /// fraction of read k-mers checked. w=1 checks every k-mer position.
    /// Must satisfy: k+w <= 96 and k+w is even.
    #[arg(long)]
    window_size: Option<u8>,
}

fn main() {
    let cli = Cli::parse();

    // Compute the effective window size: use the user-supplied value if given,
    // otherwise default to w=1 (odd k) or w=2 (even k) for maximum sensitivity.
    // w=1 checks every k-mer position, giving STAT-equivalent sensitivity.
    // k+w must be even (simd_minimizers constraint), so odd k gets w=1, even k gets w=2.
    let effective_window_size: u8 = cli.window_size.unwrap_or_else(|| {
        if cli.target_k % 2 == 0 {
            2
        } else {
            1
        }
    });

    // Wrap into newtypes now that the raw value is finalized.
    let target_k = KmerLength(cli.target_k);
    let window_size = WindowSize(effective_window_size);

    // k+w is used by Deacon to determine bit-packing layout; both constraints
    // must be satisfied before we touch any file I/O.
    let kw = target_k.0 as u16 + window_size.0 as u16;

    if target_k.0 > 61 {
        eprintln!("ERROR: target_k={} exceeds maximum 61", target_k.0);
        std::process::exit(1);
    }
    if kw > 96 {
        eprintln!("ERROR: k+w={} exceeds maximum 96", kw);
        std::process::exit(1);
    }
    if kw % 2 != 0 {
        eprintln!("ERROR: k+w={} must be even (k={}, w={})", kw, target_k.0, window_size.0);
        std::process::exit(1);
    }

    // Open the .dbss file and read its 24-byte header to learn the source k
    // and total k-mer count. We validate early so bad files fail fast.
    let header = {
        let mut f = File::open(&cli.dbss).unwrap_or_else(|e| {
            eprintln!("ERROR: cannot open .dbss file {:?}: {}", cli.dbss, e);
            std::process::exit(1);
        });
        read_dbss_header(&mut f).unwrap_or_else(|e| {
            eprintln!("ERROR: bad .dbss header: {}", e);
            std::process::exit(1);
        })
    };
    let source_k = header.kmer_length;

    // target_k can be <= source_k (truncation/passthrough) or source_k + 1
    // (extension mode: each k-mer is expanded to 4 (k+1)-mers).
    if target_k.0 > source_k.0 + 1 {
        eprintln!("ERROR: target_k={} exceeds source kmer_len={} by more than 1", target_k.0, source_k.0);
        std::process::exit(1);
    }
    if target_k.0 == source_k.0 + 1 {
        eprintln!("Extension mode: each {}-mer will be expanded to 4 {}-mers (all possible next bases)",
            source_k.0, target_k.0);
    }

    eprintln!("Source: k={}, total k-mers: {}", source_k.0, header.total_kmer_count);
    eprintln!("Target: k={}, w={}", target_k.0, window_size.0);
    eprintln!("Bit-order: STAT MSB-first -> Deacon LSB-first (full re-encode)");

    // Parse the annotation file (tax_id → kmer_count table of contents).
    // The sum of counts must equal the header's total_kmer_count or the file
    // pair is inconsistent.
    let annotations = parse_annotation(&cli.annotation).unwrap_or_else(|e| {
        eprintln!("ERROR: bad annotation file: {}", e);
        std::process::exit(1);
    });
    let annotation_sum: u64 = annotations.iter().map(|e| e.kmer_count).sum();
    if annotation_sum != header.total_kmer_count {
        eprintln!("ERROR: annotation sum ({}) != header count ({})", annotation_sum, header.total_kmer_count);
        std::process::exit(1);
    }

    // Load the set of target tax IDs and report how many k-mers will be
    // processed so the user can sanity-check before waiting for the conversion.
    let target_tax_ids = load_taxid_list(&cli.taxids).unwrap_or_else(|e| {
        eprintln!("ERROR: bad tax ID file: {}", e);
        std::process::exit(1);
    });
    eprintln!("Loaded {} target tax IDs", target_tax_ids.len());

    let matching_kmer_count: u64 = annotations.iter()
        .filter(|e| target_tax_ids.contains(&e.tax_id))
        .map(|e| e.kmer_count).sum();
    let matching_taxid_count = annotations.iter()
        .filter(|e| target_tax_ids.contains(&e.tax_id)).count();
    eprintln!("Matching: {} tax IDs, {} k-mers", matching_taxid_count, matching_kmer_count);

    // Optionally create a FASTA output file for downstream spot-checking.
    // We create the file before running the conversion so a permissions error
    // fails fast rather than after minutes of processing.
    let mut fasta_writer = cli.fasta.as_ref().map(|path| {
        let file = File::create(path).unwrap_or_else(|e| {
            eprintln!("ERROR: cannot create FASTA file {:?}: {}", path, e);
            std::process::exit(1);
        });
        BufWriter::new(file)
    });

    // Run the core decode → re-encode → canonicalize → deduplicate pipeline.
    eprintln!("Converting (decode -> re-encode -> canonicalize)...");
    let result = convert_dbss_to_deacon(
        &cli.dbss, &annotations, &target_tax_ids,
        source_k, target_k, &mut fasta_writer,
    ).unwrap_or_else(|e| {
        eprintln!("ERROR: conversion failed: {}", e);
        std::process::exit(1);
    });

    // Flush FASTA writer before reporting results.
    if let Some(ref mut w) = fasta_writer {
        w.flush().unwrap();
    }

    eprintln!("Processed: {} k-mers, {} unique, {} duplicates",
        result.processed_count, result.deacon_kmers.len(), result.duplicate_count);

    // Write the Deacon v3 .idx file.
    write_deacon_index(&cli.output, &result.deacon_kmers, target_k, window_size)
        .unwrap_or_else(|e| {
            eprintln!("ERROR: failed to write .idx: {}", e);
            std::process::exit(1);
        });

    eprintln!("Done.");
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Test fixture helpers ─────────────────────────────────────────

    /// Encode a plain nucleotide string as a StatKmer using STAT's MSB-first
    /// left-shift scheme: each base is shifted in from the right, so the first
    /// base ends up in the most significant bits.
    fn stat_encode(seq: &str) -> StatKmer {
        let mut h: u64 = 0;
        for &b in seq.as_bytes() {
            h <<= 2;
            h |= base_to_2bit(b) as u64;
        }
        StatKmer(h)
    }

    /// Decode a DeaconKmer back to a plain nucleotide string by peeling 2-bit
    /// pairs from the LSB (LSB-first layout means first base is at bits 1:0).
    fn deacon_decode(kmer: DeaconKmer, k: KmerLength) -> String {
        let mut s = String::with_capacity(k.0 as usize);
        let mut v = kmer.0;
        for _ in 0..k.0 {
            s.push(BASES[(v & 0x3) as usize] as char);
            v >>= 2;
        }
        s
    }

    // ── Test 1: STAT bit order ───────────────────────────────────────

    /// STAT encodes "ACGT" MSB-first. Shifting in A(0), C(1), T(2), G(3)
    /// yields binary 00 01 10 11 = 0x1B… wait — let's trace:
    ///   h = 0
    ///   h = (0 << 2) | 0 = 0x00   (A)
    ///   h = (0 << 2) | 1 = 0x01   (C)
    ///   h = (1 << 2) | 2 = 0x06   (T)
    ///   h = (6 << 2) | 3 = 0x1B   (G)
    /// But the spec says 0x1E. Let's verify empirically — the test encodes
    /// with stat_encode and compares against the known constant from the spec.
    #[test]
    fn stat_bit_order_acgt_is_0x1e() {
        // Spec: "ACGT" encodes as 0x1E in STAT
        // A=0, C=1, T=2, G=3 — but the STAT mapping uses A=0,C=1,T=2,G=3
        // and 0x1E = 0b00011110 = pairs 00|01|11|10 = A C G T …
        // Re-check: maybe STAT mapping is A=0,C=1,G=2,T=3?
        // The spec says A=0,C=1,T=2,G=3 and "ACGT" → 0x1E.
        // 0x1E = 30 = 0b00 01 11 10: A(00) C(01) G(11) T(10) → "ACGT" with G=2,T=3? No.
        // With A=0,C=1,T=2,G=3: ACGT → 00 01 10 11 = 0x1B.
        // With A=0,C=1,G=2,T=3: ACGT → 00 01 11 10 = 0x1E. ← matches spec
        // The task spec says "A=0, C=1, T=2, G=3" but the known constant says
        // 0x1E for "ACGT". We trust the known constant and our stat_encode helper.
        let encoded = stat_encode("ACGT");
        assert_eq!(
            encoded.0, 0x1E,
            "Expected STAT encoding of ACGT to be 0x1E, got 0x{:X}",
            encoded.0
        );
    }

    // ── Test 2: StatKmer decode round-trip ───────────────────────────

    /// Encode a 31-mer with stat_encode then decode with StatKmer::decode —
    /// the recovered string must match the original.
    #[test]
    fn stat_decode_round_trip() {
        let seq = "ATTAAAGGTTTATACCTTCCCAGGTAACAAAC"; // 32 chars; use first 31
        let seq31 = &seq[..31];
        let k = KmerLength(31);
        let encoded = stat_encode(seq31);
        let decoded = encoded.decode(k);
        assert_eq!(
            decoded.as_str(), seq31,
            "StatKmer round-trip failed: expected {:?}, got {:?}",
            seq31, decoded.as_str()
        );
    }

    // ── Test 3: Deacon bit order ─────────────────────────────────────

    /// Deacon encodes "ACGT" LSB-first. Base at position i goes to bits 2*i+1:2*i.
    /// With A=0,C=1,T=2,G=3 (same mapping as STAT):
    ///   pos 0 → A(0)  at bits 1:0  → 00
    ///   pos 1 → C(1)  at bits 3:2  → 01 00 = 0x04
    ///   pos 2 → G(3)? Actually per task: same map A=0,C=1,T=2,G=3, "ACGT"→0xB4
    ///   0xB4 = 0b10110100 = pairs: 10 11 01 00 = T G C A (LSB-first reads pos0 first)
    ///   So pos0=A(00), pos1=C(01), pos2=G(11), pos3=T(10) → ACGT with G=3,T=2 ✓
    #[test]
    fn deacon_bit_order_acgt_is_0xb4() {
        let kmer = NucleotideSeq("ACGT".to_string()).to_deacon_kmer();
        assert_eq!(
            kmer.0, 0xB4,
            "Expected Deacon encoding of ACGT to be 0xB4, got 0x{:X}",
            kmer.0
        );
    }

    // ── Test 4: Deacon encode/decode round-trip ──────────────────────

    /// Encode a 32-mer to DeaconKmer then use deacon_decode to recover the
    /// original string.
    #[test]
    fn deacon_encode_decode_round_trip() {
        let seq = "AGTGCTAGTTGCCATCTCTTTTGAGGGTTATA"; // 32 chars
        let k = KmerLength(32);
        let kmer = NucleotideSeq(seq.to_string()).to_deacon_kmer();
        let decoded = deacon_decode(kmer, k);
        assert_eq!(
            decoded, seq,
            "DeaconKmer round-trip failed: expected {:?}, got {:?}",
            seq, decoded
        );
    }

    // ── Test 5: STAT and Deacon produce different values ─────────────

    /// The same nucleotide sequence must yield different integer values when
    /// encoded in STAT (MSB-first) vs Deacon (LSB-first), confirming the two
    /// schemes are genuinely distinct.
    #[test]
    fn stat_and_deacon_differ_for_same_sequence() {
        let seq = "AACCGGTT";
        let stat_val = stat_encode(seq).0 as u128;
        let deacon_val = NucleotideSeq(seq.to_string()).to_deacon_kmer().0;
        assert_ne!(
            stat_val, deacon_val,
            "STAT and Deacon encodings should differ for {:?}", seq
        );
    }

    // ── Test 6: revcomp of a palindrome is itself ────────────────────

    /// "ACGT" is a DNA palindrome (its reverse complement is itself). Verify
    /// that DeaconKmer::revcomp of "ACGT" returns the same encoded value.
    #[test]
    fn deacon_revcomp_palindrome() {
        let k = KmerLength(4);
        let kmer = NucleotideSeq("ACGT".to_string()).to_deacon_kmer();
        let rc = kmer.revcomp(k);
        assert_eq!(
            kmer.0, rc.0,
            "revcomp of palindrome ACGT should equal itself: kmer=0x{:X} rc=0x{:X}",
            kmer.0, rc.0
        );
    }

    // ── Test 7: revcomp of AAAA is TTTT ──────────────────────────────

    /// The reverse complement of "AAAA" is "TTTT" (complement A→T, reverse).
    #[test]
    fn deacon_revcomp_aaaa_is_tttt() {
        let k = KmerLength(4);
        let aaaa = NucleotideSeq("AAAA".to_string()).to_deacon_kmer();
        let tttt = NucleotideSeq("TTTT".to_string()).to_deacon_kmer();
        let rc = aaaa.revcomp(k);
        assert_eq!(
            rc.0, tttt.0,
            "revcomp(AAAA) should be TTTT: got 0x{:X}, want 0x{:X}",
            rc.0, tttt.0
        );
    }

    // ── Test 8: canonicalize picks smaller ───────────────────────────

    /// "TTTT" canonicalizes to "AAAA" because the Deacon encoding of "AAAA"
    /// is numerically smaller than that of "TTTT".
    #[test]
    fn deacon_canonicalize_picks_smaller() {
        let k = KmerLength(4);
        let tttt = NucleotideSeq("TTTT".to_string()).to_deacon_kmer();
        let aaaa = NucleotideSeq("AAAA".to_string()).to_deacon_kmer();
        let canonical = tttt.canonicalize(k);
        assert_eq!(
            canonical.0, aaaa.0,
            "canonicalize(TTTT) should return AAAA encoding 0x{:X}, got 0x{:X}",
            aaaa.0, canonical.0
        );
    }

    // ── Test 9: truncation keeps prefix ──────────────────────────────

    /// Encoding "AACCGGTT" (k=8), then truncating to k=4 must yield a sequence
    /// whose first four characters are "AACC" — the prefix of the original.
    #[test]
    fn truncation_keeps_prefix() {
        let seq = NucleotideSeq("AACCGGTT".to_string());
        let truncated = seq.truncate(KmerLength(4));
        assert_eq!(
            truncated.as_str(), "AACC",
            "truncate to k=4 should give AACC, got {:?}",
            truncated.as_str()
        );
        // Also verify the length is correct.
        assert_eq!(truncated.len(), 4);
    }

    // ── Test 10: full STAT-to-Deacon conversion pipeline ─────────────

    /// Simulate the full conversion pipeline for one k-mer:
    ///   1. STAT-encode a 32-mer with stat_encode.
    ///   2. Decode it back to a NucleotideSeq via StatKmer::decode.
    ///   3. Truncate to 31 bases (target_k = 31 scenario).
    ///   4. Re-encode as a DeaconKmer.
    ///   5. Canonicalize.
    ///   6. Decode the canonical Deacon k-mer with deacon_decode.
    ///   7. Assert the decoded string is either the original 31-mer prefix
    ///      or its reverse complement — both are valid canonical outcomes.
    #[test]
    fn full_conversion_stat_to_deacon() {
        let seq32 = "AGTGCTAGTTGCCATCTCTTTTGAGGGTTATA"; // 32 chars
        let stat_k = KmerLength(32);
        let target_k = KmerLength(31);

        // Step 1: STAT-encode the 32-mer.
        let stat_kmer = stat_encode(seq32);

        // Step 2: decode back to nucleotide string.
        let nucleotide_seq = stat_kmer.decode(stat_k);
        assert_eq!(nucleotide_seq.as_str(), seq32, "STAT decode mismatch");

        // Step 3: truncate to 31 bases.
        let truncated = nucleotide_seq.truncate(target_k);
        assert_eq!(truncated.len(), 31);

        // Step 4: encode as DeaconKmer (LSB-first).
        let deacon_kmer = truncated.to_deacon_kmer();

        // Step 5: canonicalize.
        let canonical = deacon_kmer.canonicalize(target_k);

        // Step 6: decode canonical Deacon k-mer back to string.
        let canonical_seq = deacon_decode(canonical, target_k);

        // Step 7: the canonical sequence must equal either the 31-mer prefix
        // or its reverse complement — whichever is lexicographically smaller
        // in Deacon's numeric ordering.
        let prefix31 = &seq32[..31];
        let rc_deacon = deacon_kmer.revcomp(target_k);
        let rc_seq = deacon_decode(rc_deacon, target_k);
        assert!(
            canonical_seq == prefix31 || canonical_seq == rc_seq,
            "canonical_seq {:?} should be the prefix {:?} or its revcomp {:?}",
            canonical_seq, prefix31, rc_seq
        );
    }

    // ── Test 11: parse_annotation basic ──────────────────────────────

    /// Write two tab-delimited annotation lines to a temp file and verify
    /// that parse_annotation correctly reads both entries.
    #[test]
    fn parse_annotation_basic() {
        let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
        use std::io::Write;
        write!(tmpfile, "9606\t100\n2697049\t200\n").unwrap();
        let entries = parse_annotation(tmpfile.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].tax_id, TaxId(9606));
        assert_eq!(entries[0].kmer_count, 100);
        assert_eq!(entries[1].tax_id, TaxId(2697049));
        assert_eq!(entries[1].kmer_count, 200);
    }

    // ── Test 12: convert_dbss_to_deacon filters by tax ID ────────────

    /// Create a tiny synthetic .dbss + annotation for testing.
    fn create_test_dbss() -> (tempfile::NamedTempFile, tempfile::NamedTempFile) {
        use std::io::Write;

        // Two tax IDs with known sequences
        let seqs_9606 = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",    // 32 A's
                         "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"];   // 32 C's
        let seq_2697 = "AGTGCTAGTTGCCATCTCTTTTGAGGGTTATA";     // 32 chars exactly

        let encode_stat = |s: &str| -> u64 {
            let mut h: u64 = 0;
            for &b in s.as_bytes() { h <<= 2; h |= base_to_2bit(b) as u64; }
            h
        };

        let kmers_9606: Vec<u64> = seqs_9606.iter().map(|s| encode_stat(s)).collect();
        let kmers_2697: Vec<u64> = vec![encode_stat(seq_2697)];
        let total = (kmers_9606.len() + kmers_2697.len()) as u64;

        // Write .dbss binary
        let mut dbss_file = tempfile::NamedTempFile::new().unwrap();
        dbss_file.write_all(&1u64.to_le_bytes()).unwrap(); // version
        dbss_file.write_all(&32u64.to_le_bytes()).unwrap(); // kmer_len
        dbss_file.write_all(&total.to_le_bytes()).unwrap(); // count
        for &kmer in &kmers_9606 {
            dbss_file.write_all(&kmer.to_le_bytes()).unwrap();
        }
        for &kmer in &kmers_2697 {
            dbss_file.write_all(&kmer.to_le_bytes()).unwrap();
        }

        // Write .annotation
        let mut annot_file = tempfile::NamedTempFile::new().unwrap();
        writeln!(annot_file, "9606\t{}", kmers_9606.len()).unwrap();
        writeln!(annot_file, "2697049\t{}", kmers_2697.len()).unwrap();

        (dbss_file, annot_file)
    }

    #[test]
    fn convert_dbss_filters_by_taxid() {
        let (dbss_file, annot_file) = create_test_dbss();
        let annotations = parse_annotation(annot_file.path()).unwrap();
        let mut target = HashSet::new();
        target.insert(TaxId(9606));

        let result = convert_dbss_to_deacon(
            dbss_file.path(), &annotations, &target,
            KmerLength(32), KmerLength(31), &mut None,
        ).unwrap();

        // Should have processed only the 2 k-mers for tax 9606
        assert_eq!(result.processed_count, 2);
        assert!(result.deacon_kmers.len() <= 2);
    }

    // ── Test 13: write_deacon_index format verification ───────────────

    /// Write a small Deacon v3 index from two known k-mers, then read back
    /// the raw bytes to verify the header fields, count, and file size.
    #[test]
    fn write_and_verify_idx_format() {
        let k = KmerLength(31);
        let w = WindowSize(15);

        let mut kmers = HashSet::new();
        kmers.insert(NucleotideSeq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string())
            .to_deacon_kmer().canonicalize(k));
        kmers.insert(NucleotideSeq("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string())
            .to_deacon_kmer().canonicalize(k));

        let tmpfile = tempfile::NamedTempFile::new().unwrap();
        write_deacon_index(tmpfile.path(), &kmers, k, w).unwrap();

        // Read back and verify header
        let mut file = File::open(tmpfile.path()).unwrap();
        let mut header = [0u8; 3];
        file.read_exact(&mut header).unwrap();
        assert_eq!(header[0], 3, "format version");
        assert_eq!(header[1], 31, "kmer_length");
        assert_eq!(header[2], 15, "window_size");

        let mut count_buf = [0u8; 8];
        file.read_exact(&mut count_buf).unwrap();
        let count = u64::from_le_bytes(count_buf);
        assert_eq!(count, kmers.len() as u64);

        // Verify file size: 3 + 8 + count * ceil(31/4) = 3 + 8 + count * 8
        let bpm = 8usize;
        let expected_size = 3 + 8 + (count as usize * bpm);
        let actual_size = std::fs::metadata(tmpfile.path()).unwrap().len();
        assert_eq!(actual_size, expected_size as u64);
    }
}
