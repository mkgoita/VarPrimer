#!/usr/bin/env python3
"""
primer_design.py
================
Design PCR primers around genomic variant positions.

Input : CSV file with columns  gene, position  (e.g. BRCA1, chr17:43044295)
Output: CSV file with columns  Name, Sequence, TM, GC_percent, ProductSize, Exon

Usage example
-------------
python primer_design.py \
    --input  variants.csv \
    --output primers_output.csv \
    --email  your@email.com \
    --build  hg38
"""

import argparse
import csv
import logging
import os
import sys
import time
import warnings
import urllib.request
import urllib.parse
from math import log as math_log
from typing import Optional, Tuple, List, Dict

import pandas as pd
import requests
from Bio import Entrez, SeqIO

# ── Optional primer3-py ────────────────────────────────────────────────────────
try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False
    warnings.warn(
        "primer3-py not found – using built-in pure-Python primer picker. "
        "Install primer3-py for better primer design.",
        ImportWarning,
        stacklevel=2,
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────
NCBI_FETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_SUMMARY_URL= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ENSEMBL_REST    = "https://rest.ensembl.org"
NCBI_SLEEP      = 0.4          # seconds between NCBI requests (≤3/s without API key)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Sequence / Thermodynamics helpers
# ─────────────────────────────────────────────────────────────────────────────

def gc_percent(seq: str) -> float:
    """Return GC% of a DNA sequence (0-100)."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return round(100.0 * gc / len(seq), 1) if seq else 0.0


def tm_basic(seq: str) -> float:
    """
    Basic Tm calculation.
    Uses the Wallace rule for short oligos (< 14 nt),
    the nearest-neighbour approximation for longer ones.
    """
    seq = seq.upper().replace("U", "T")
    n   = len(seq)
    if n == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")

    if n < 14:
        return float(4 * gc + 2 * at)

    # SantaLucia 1998 nearest-neighbour unified parameters (kcal/mol, cal/mol/K)
    NN_DH = {
        "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2,
        "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,
        "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
        "CG": -10.6,"GC": -9.8, "GG": -8.0, "CC": -8.0,
    }
    NN_DS = {
        "AA": -22.2,"TT": -22.2,"AT": -20.4,"TA": -21.3,
        "CA": -22.7,"TG": -22.7,"GT": -22.4,"AC": -22.4,
        "CT": -21.0,"AG": -21.0,"GA": -22.2,"TC": -22.2,
        "CG": -27.2,"GC": -24.4,"GG": -19.9,"CC": -19.9,
    }
    dH = sum(NN_DH.get(seq[i:i+2], -8.0) for i in range(n - 1))
    dS = sum(NN_DS.get(seq[i:i+2], -21.0) for i in range(n - 1))
    # Add initiation parameters
    if seq[0] in "GC":
        dH += 0.1;  dS += -2.8
    else:
        dH += 2.3;  dS += 4.1
    if seq[-1] in "GC":
        dH += 0.1;  dS += -2.8
    else:
        dH += 2.3;  dS += 4.1

    R        = 1.987          # cal/(mol·K)
    oligo_nM = 250            # primer concentration nM
    dS_total = dS - R * math_log(1 / (oligo_nM * 1e-9))
    # Tm = ΔH / ΔS  (convert dH kcal→cal, dS already in cal)
    tm = (dH * 1000) / dS_total - 273.15
    return round(tm, 1)


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Exon annotation via Ensembl REST
# ─────────────────────────────────────────────────────────────────────────────

def get_exons_ensembl(gene_name: str, build: str) -> Tuple[List[Dict], str]:
    """
    Return (exons, transcript_id) where exons is a sorted list of:
        { exon_number, chrom, start, end, strand }
    and transcript_id is the Ensembl ID of the selected transcript
    (MANE Select preferred).
    """
    server = "https://grch37.rest.ensembl.org" if build in ("hg19","GRCh37") \
             else "https://rest.ensembl.org"

    # 1. Look up the gene ID
    headers = {"Content-Type": "application/json"}
    r = requests.get(
        f"{server}/xrefs/symbol/homo_sapiens/{gene_name}",
        headers=headers
    )
    r.raise_for_status()
    hits = [h for h in r.json() if h.get("type") == "gene"]
    if not hits:
        raise ValueError(f"Gene '{gene_name}' not found in Ensembl ({build}).")
    gene_id = hits[0]["id"]

    # 2. Fetch transcripts with MANE annotation
    time.sleep(0.2)
    r2 = requests.get(
        f"{server}/lookup/id/{gene_id}",
        headers=headers,
        params={"expand": 1, "mane": 1, "content-type": "application/json"},
    )
    r2.raise_for_status()
    gene_data = r2.json()

    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        raise ValueError(f"No transcripts found for gene '{gene_name}'.")

    # Transcript selection priority:
    #   1. MANE Select  (gold-standard clinical transcript — most clinically relevant)
    #   2. Ensembl canonical (is_canonical == 1)
    #   3. Longest transcript (fallback)
    canon = None

    # Priority 1 — MANE Select
    for t in transcripts:
        mane_list = t.get("MANE", [])
        if any(m.get("type") == "MANE_Select" for m in mane_list):
            canon = t
            log.info("Using MANE Select transcript for %s: %s (%s)",
                     gene_name, t.get("id"), t.get("display_name", ""))
            break

    # Priority 2 — Ensembl canonical
    if canon is None:
        for t in transcripts:
            if t.get("is_canonical") == 1:
                canon = t
                log.info("Using canonical transcript for %s (no MANE Select): %s (%s)",
                         gene_name, t.get("id"), t.get("display_name", ""))
                break

    # Priority 3 — longest
    if canon is None:
        canon = max(transcripts, key=lambda t: t.get("length", 0))
        log.warning("No MANE Select or canonical transcript for %s; using longest: %s",
                    gene_name, canon.get("id"))

    transcript_id = canon.get("id", "")
    strand = canon.get("strand", 1)
    chrom  = canon.get("seq_region_name", "")

    exons_raw = canon.get("Exon", [])
    if strand == 1:
        exons_sorted = sorted(exons_raw, key=lambda e: e["start"])
    else:
        exons_sorted = sorted(exons_raw, key=lambda e: e["start"], reverse=True)

    exons = []
    for i, ex in enumerate(exons_sorted, start=1):
        exons.append({
            "exon_number": i,
            "chrom":       chrom,
            "start":       int(ex["start"]),
            "end":         int(ex["end"]),
            "strand":      strand,
        })

    return exons, transcript_id


def find_exon(position: int, exons: List[Dict]) -> int:
    """
    Return the exon number that contains the position.
    If the position is intronic, return the number of the nearest exon.
    """
    # Exact match
    for ex in exons:
        if ex["start"] <= position <= ex["end"]:
            return ex["exon_number"]

    # Nearest exon (by minimum distance to any exon boundary)
    best_num  = exons[0]["exon_number"]
    best_dist = float("inf")
    for ex in exons:
        dist = min(abs(position - ex["start"]), abs(position - ex["end"]))
        if dist < best_dist:
            best_dist = dist
            best_num  = ex["exon_number"]
    log.warning(
        "Position %d is intronic; assigning the nearest exon (%d)", position, best_num
    )
    return best_num


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Fetch genomic sequence — primary: Ensembl REST; fallback: NCBI Entrez
# ─────────────────────────────────────────────────────────────────────────────
# Using Ensembl REST guarantees the sequence is byte-for-byte identical to
# what the Ensembl genome browser displays, so primers will be findable there.

# NCBI accession map kept as fallback only
CHROM_ACC_HG38 = {
    "chr1":  "NC_000001.11", "chr2":  "NC_000002.12", "chr3":  "NC_000003.12",
    "chr4":  "NC_000004.12", "chr5":  "NC_000005.10", "chr6":  "NC_000006.12",
    "chr7":  "NC_000007.14", "chr8":  "NC_000008.11", "chr9":  "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9",  "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX":  "NC_000023.11", "chrY":  "NC_000024.10",
}
CHROM_ACC_HG19 = {
    "chr1":  "NC_000001.10", "chr2":  "NC_000002.11", "chr3":  "NC_000003.11",
    "chr4":  "NC_000004.11", "chr5":  "NC_000005.9",  "chr6":  "NC_000006.11",
    "chr7":  "NC_000007.13", "chr8":  "NC_000008.10", "chr9":  "NC_000009.11",
    "chr10": "NC_000010.10", "chr11": "NC_000011.9",  "chr12": "NC_000012.11",
    "chr13": "NC_000013.10", "chr14": "NC_000014.8",  "chr15": "NC_000015.9",
    "chr16": "NC_000016.9",  "chr17": "NC_000017.10", "chr18": "NC_000018.9",
    "chr19": "NC_000019.9",  "chr20": "NC_000020.10", "chr21": "NC_000021.8",
    "chr22": "NC_000022.10", "chrX":  "NC_000023.10", "chrY":  "NC_000024.9",
}


def chrom_to_accession(chrom: str, build: str) -> str:
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    acc_map = CHROM_ACC_HG19 if build in ("hg19", "GRCh37") else CHROM_ACC_HG38
    acc = acc_map.get(chrom)
    if not acc:
        raise ValueError(f"Unknown chromosome '{chrom}' for build '{build}'.")
    return acc


def fetch_sequence_ensembl(chrom: str, start: int, end: int, build: str) -> str:
    """
    Fetch genomic sequence via Ensembl REST (positive-strand, 1-based inclusive).
    Returns uppercase sequence string.
    This is the primary sequence source — identical to what Ensembl shows.
    """
    server = "https://grch37.rest.ensembl.org" if build in ("hg19", "GRCh37") \
             else "https://rest.ensembl.org"
    chrom_bare = chrom.lstrip("chr")   # Ensembl uses bare numbers: "2" not "chr2"
    url = f"{server}/sequence/region/human/{chrom_bare}:{start}..{end}:1"
    r = requests.get(url, headers={"Content-Type": "text/plain"})
    r.raise_for_status()
    seq = r.text.strip().upper()
    if not seq:
        raise ValueError(f"Ensembl returned empty sequence for {chrom}:{start}-{end}.")
    return seq


def fetch_sequence_ncbi(chrom: str, start: int, end: int, build: str,
                        email: str, api_key: Optional[str] = None) -> str:
    """NCBI Entrez fallback — used if Ensembl REST fails."""
    acc = chrom_to_accession(chrom, build)
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.efetch(
        db        = "nucleotide",
        id        = acc,
        rettype   = "fasta",
        retmode   = "text",
        seq_start = start,
        seq_stop  = end,
    )
    record = SeqIO.read(handle, "fasta")
    handle.close()
    time.sleep(NCBI_SLEEP)
    return str(record.seq).upper()


def fetch_sequence(chrom: str, start: int, end: int, build: str, email: str,
                   api_key: Optional[str] = None) -> str:
    """
    Fetch genomic sequence [start, end] (1-based inclusive).
    Tries Ensembl REST first (same source as Ensembl genome browser),
    falls back to NCBI Entrez on failure.
    Returns uppercase sequence string.
    """
    try:
        time.sleep(0.1)                  # be polite to Ensembl
        seq = fetch_sequence_ensembl(chrom, start, end, build)
        log.debug("Sequence fetched from Ensembl REST (%d bp).", len(seq))
        return seq
    except Exception as ens_exc:
        log.warning("Ensembl sequence fetch failed (%s); falling back to NCBI.", ens_exc)

    # NCBI fallback
    seq = fetch_sequence_ncbi(chrom, start, end, build, email, api_key)
    log.debug("Sequence fetched from NCBI (%d bp).", len(seq))
    return seq


def check_ucsc_insilico_pcr(fwd_seq: str, rev_seq: str, build: str = "hg38") -> int:
    """
    Submits primers to UCSC In-Silico PCR (hgPcr) and returns the number of expected amplicons.
    Uses 'genome-test.gi.ucsc.edu' since it often avoids the Cloudflare blocks on the primary server.
    Returns 0 if no match or error, 1 if unique, >1 if there are off-target amplicons.
    Returns -1 if the query fails completely.
    """
    import subprocess
    import bs4
    
    # We use 'genome-test' mirror which is typically more accessible for the CGI tools
    # We still use 'Submit=Submit' (capital S) as it is standard for hgPcr
    url_base = f"https://genome-test.gi.ucsc.edu/cgi-bin/hgPcr?db={build}"
    url_params = (
        f"&wp_target=genome&wp_f={fwd_seq}&wp_r={rev_seq}"
        f"&Submit=Submit&wp_size=4000&wp_perfect=15&wp_good=15"
        f"&boolshad.wp_flipReverse=0&wp_append=on&boolshad.wp_append=0"
    )
    
    try:
        # Using curl with a generic browser User-Agent
        cmd = ["curl", "-s", "-L", "-A", "Mozilla/5.0", url_base + url_params]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        html = result.stdout
        
        # Check for Turnstile/Cloudflare challenge specifically
        # Challenge pages are short and contain specific script markers
        is_blocked = len(html) < 2000 and ('turnstile' in html.lower() or 'cf-browser-verification' in html.lower())
        
        if not html or is_blocked:
            # Fallback to primary server if mirror is down or blocked
            url_fallback = f"https://genome.ucsc.edu/cgi-bin/hgPcr?db={build}" + url_params
            cmd[-1] = url_fallback
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            html = result.stdout
            is_blocked = len(html) < 2000 and ('turnstile' in html.lower() or 'cf-browser-verification' in html.lower())
            
        if not html or is_blocked:
            log.warning("UCSC: Blocked by Cloudflare on both primary and mirror.")
            return -1

        # If "No matches" is returned, amplicons = 0
        if 'No matches to' in html:
            return 0
            
        soup = bs4.BeautifulSoup(html, 'html.parser')
        pre_tags = soup.find_all(['pre', 'PRE'])
            
        if pre_tags:
            # Each hit is an anchor link within the PRE block
            hits = pre_tags[0].find_all('a')
            if hits:
                return len(hits)
            
            # Fallback: Count lines starting with '>' in the PRE block
            text = pre_tags[0].get_text()
            lines = [l for l in text.split('\n') if l.strip().startswith('>')]
            if lines:
                return len(lines)
            
        return 0
    except Exception as exc:
        log.warning("UCSC In-Silico PCR check failed: %s", exc)
        return -1


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Primer design
# ─────────────────────────────────────────────────────────────────────────────

def design_primers_primer3(
    sequence: str,
    target_offset: int,
    target_len: int = 1,
    product_size_range: Tuple[int, int] = (100, 500),
    tm_min: float = 57.0,
    tm_opt: float = 62.0,
    tm_max: float = 67.0,
    gc_min: float = 40.0,
    gc_max: float = 65.0,
    primer_len_min: int = 18,
    primer_len_opt: int = 20,
    primer_len_max: int = 25,
) -> List[Dict]:
    """
    Use primer3-py to design primer pairs.
    Returns a list of dicts (up to 5) with left/right sequences, TMs, GCs, product size.
    """
    result = primer3.design_primers(
        seq_args={
            "SEQUENCE_ID":      "target",
            "SEQUENCE_TEMPLATE": sequence,
            "SEQUENCE_TARGET":  [target_offset, target_len],
        },
        global_args={
            "PRIMER_OPT_SIZE":            primer_len_opt,
            "PRIMER_MIN_SIZE":            primer_len_min,
            "PRIMER_MAX_SIZE":            primer_len_max,
            "PRIMER_OPT_TM":              tm_opt,
            "PRIMER_MIN_TM":              tm_min,
            "PRIMER_MAX_TM":              tm_max,
            "PRIMER_MIN_GC":              gc_min,
            "PRIMER_MAX_GC":              gc_max,
            "PRIMER_PRODUCT_SIZE_RANGE":  [list(product_size_range)],
            "PRIMER_NUM_RETURN":          5,
            "PRIMER_EXPLAIN_FLAG":        1,
        },
    )

    pairs = []
    for i in range(5):
        key_l = f"PRIMER_LEFT_{i}_SEQUENCE"
        key_r = f"PRIMER_RIGHT_{i}_SEQUENCE"
        if key_l in result and key_r in result:
            left_seq  = result[key_l]
            right_seq = result[key_r]
            prod_size = result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)
            left_pos  = result.get(f"PRIMER_LEFT_{i}",  [0, 0])   # [start, length]
            right_pos = result.get(f"PRIMER_RIGHT_{i}", [0, 0])   # [3'-end, length]
            pairs.append({
                "left_seq":    left_seq,
                "right_seq":   right_seq,
                "left_tm":     round(result.get(f"PRIMER_LEFT_{i}_TM",  tm_basic(left_seq)), 1),
                "right_tm":    round(result.get(f"PRIMER_RIGHT_{i}_TM", tm_basic(right_seq)), 1),
                "left_gc":     round(result.get(f"PRIMER_LEFT_{i}_GC_PERCENT",  gc_percent(left_seq)), 1),
                "right_gc":    round(result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", gc_percent(right_seq)), 1),
                "product_size": prod_size,
                # Offsets within the sequence window (0-based)
                "left_start":  left_pos[0],
                "left_end":    left_pos[0] + left_pos[1] - 1,
                "right_end":   right_pos[0],                          # 3'-end (rightmost)
                "right_start": right_pos[0] - right_pos[1] + 1,      # 5'-end of rev primer
            })
    return pairs


def design_primers_builtin(
    sequence: str,
    target_offset: int,
    product_size_range: Tuple[int, int] = (100, 500),
    tm_min: float = 55.0,
    tm_max: float = 68.0,
    gc_min: float = 40.0,
    gc_max: float = 65.0,
    primer_len_min: int = 18,
    primer_len_max: int = 25,
) -> List[Dict]:
    """
    Pure-Python greedy primer picker.
    Slides a window on the left of the target and looks for a complementary
    window on the right that meets TM / GC / product-size constraints.
    Returns a list of candidate pairs sorted by score (TM distance from 62 °C).
    """
    min_prod, max_prod = product_size_range

    def is_valid(seq: str) -> bool:
        tm  = tm_basic(seq)
        gc  = gc_percent(seq)
        return (tm_min <= tm <= tm_max) and (gc_min <= gc <= gc_max)

    candidates = []

    for fwd_end in range(primer_len_min - 1, min(target_offset, len(sequence))):
        for fwd_start in range(max(0, fwd_end - primer_len_max + 1), fwd_end - primer_len_min + 2):
            fwd_seq = sequence[fwd_start: fwd_end + 1]
            if not is_valid(fwd_seq):
                continue

            # Right primer: its 3' end must be ≥ target_offset
            rev_start_min = fwd_start + min_prod - 1
            rev_start_max = fwd_start + max_prod - 1

            for rev_end in range(
                min(rev_start_min + primer_len_min - 1, len(sequence) - 1),
                min(rev_start_max + primer_len_max,     len(sequence)),
            ):
                for rev_start in range(
                    max(rev_end - primer_len_max + 1, target_offset),
                    rev_end - primer_len_min + 2,
                ):
                    rev_template = sequence[rev_start: rev_end + 1]
                    rev_seq = reverse_complement(rev_template)
                    if not is_valid(rev_seq):
                        continue

                    prod = rev_end - fwd_start + 1
                    if not (min_prod <= prod <= max_prod):
                        continue

                    # Score: distance of both TMs from 62 °C
                    score = abs(tm_basic(fwd_seq) - 62) + abs(tm_basic(rev_seq) - 62)
                    candidates.append((score, {
                        "left_seq":    fwd_seq,
                        "right_seq":   rev_seq,
                        "left_tm":     tm_basic(fwd_seq),
                        "right_tm":    tm_basic(rev_seq),
                        "left_gc":     gc_percent(fwd_seq),
                        "right_gc":    gc_percent(rev_seq),
                        "product_size": prod,
                        "left_start":  fwd_start,
                        "left_end":    fwd_end,
                        "right_start": rev_start,
                        "right_end":   rev_end,
                    }))
                    
                    # Keep only the top 10 candidates to save memory
                    if len(candidates) > 50:
                        candidates.sort(key=lambda x: x[0])
                        candidates = candidates[:10]

    candidates.sort(key=lambda x: x[0])
    return [c[1] for c in candidates[:10]]


def design_primers(
    sequence: str,
    target_offset: int,
    **kwargs,
) -> List[Dict]:
    """Dispatch to primer3-py or built-in picker."""
    if HAS_PRIMER3:
        return design_primers_primer3(sequence, target_offset, **kwargs)
    else:
        # Remove primer3-specific args that the built-in doesn't accept
        kwargs.pop("target_len", None)
        kwargs.pop("tm_opt", None)
        kwargs.pop("primer_len_opt", None)
        return design_primers_builtin(sequence, target_offset, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Main per-variant workflow
# ─────────────────────────────────────────────────────────────────────────────

# cDNA position pattern: c.617, c.617G>A, c.617+2A>G, c.-14G>A, c.*5A>G etc.
import re as _re
_CDNA_RE = _re.compile(r'^c\.([\-\*]?\d+)', _re.IGNORECASE)


def is_cdna_position(pos_str: str) -> bool:
    """Return True if pos_str looks like a cDNA position (e.g. c.617G>A, c.617)."""
    return bool(_CDNA_RE.match(pos_str.strip()))


def parse_cdna_pos(pos_str: str) -> int:
    """
    Extract the numeric cDNA coordinate from an HGVS-like string.
    Examples:
        c.617G>A  -> 617
        c.617     -> 617
        c.*5A>G   -> 5   (UTR/stop-adjacent variants)
        c.-14G>A  -> -14 (5'UTR)
    Returns the absolute value for intronic positions like c.617+2 -> 617.
    """
    m = _CDNA_RE.match(pos_str.strip())
    if not m:
        raise ValueError(f"Cannot parse cDNA position from '{pos_str}'.")
    offset_str = m.group(1)
    # Drop the '*' prefix (stop-adjacent) — treat as positive
    return abs(int(offset_str.replace('*', '')))


def nm_to_ensembl_transcript(nm_id: str, build: str) -> str:
    """
    Convert a RefSeq transcript accession (NM_ / NR_) to an Ensembl ENST ID
    via the Ensembl xrefs REST API.  Returns the ENST ID as a string.
    Raises ValueError if no mapping is found.
    """
    server = "https://grch37.rest.ensembl.org" if build in ("hg19", "GRCh37") \
             else "https://rest.ensembl.org"
    # Strip version (NM_001105.5  →  NM_001105)
    nm_base = nm_id.split(".")[0]
    url = f"{server}/xrefs/symbol/homo_sapiens/{nm_base}"
    r = requests.get(url, params={"external_db": "RefSeq_mRNA"},
                     headers={"Content-Type": "application/json"})
    r.raise_for_status()
    hits = [h for h in r.json() if h.get("type") == "transcript"]
    if not hits:
        raise ValueError(f"No Ensembl transcript found for RefSeq ID '{nm_id}'.")
    enst = hits[0]["id"]
    log.info("Resolved %s → %s", nm_id, enst)
    return enst


def cdna_to_genomic(
    transcript_id: str,
    cdna_pos: int,
    build: str,
) -> Tuple[str, int, int]:
    """
    Convert a cDNA position to genomic coordinates via Ensembl REST.
    transcript_id may be an ENST ID or a RefSeq NM_/NR_ accession
    (it will be resolved automatically).
    Returns (chrom, genomic_start, genomic_end) — 1-based.
    """
    server = "https://grch37.rest.ensembl.org" if build in ("hg19", "GRCh37") \
             else "https://rest.ensembl.org"

    # Resolve NM_ / NR_ accessions to ENST
    enst_id = transcript_id
    if transcript_id.upper().startswith(("NM_", "NR_", "XM_")):
        try:
            enst_id = nm_to_ensembl_transcript(transcript_id, build)
        except Exception as exc:
            raise ValueError(
                f"Could not resolve RefSeq ID '{transcript_id}' to Ensembl: {exc}"
            ) from exc

    url = f"{server}/map/cdna/{enst_id}/{cdna_pos}..{cdna_pos}"
    r = requests.get(url, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    data = r.json()

    mappings = data.get("mappings", [])
    if not mappings:
        raise ValueError(
            f"Ensembl could not map cDNA position {cdna_pos} on transcript {enst_id}."
        )
    m = mappings[0]
    chrom     = m["seq_region_name"]
    gen_start = int(m["start"])
    gen_end   = int(m["end"])
    log.info("Mapped cDNA pos %d on %s → chr%s:%d-%d",
             cdna_pos, enst_id, chrom, gen_start, gen_end)
    return chrom, gen_start, gen_end


def parse_position(pos_str: str) -> Tuple[str, int]:
    """
    Parse GENOMIC position strings such as:
        chr17:43044295  |  17:43044295  |  chr17:43,044,295
    Returns (chrom, pos_1based).
    For cDNA input, use cdna_to_genomic() instead.
    """
    pos_str = pos_str.replace(",", "").strip()
    if ":" not in pos_str:
        raise ValueError(f"Cannot parse position '{pos_str}'. Expected chr:pos or c.NNN.")
    chrom, pos = pos_str.split(":", 1)
    return chrom.strip(), int(pos.strip())


def process_variant(
    gene: str,
    position_str: str,
    build: str = "hg38",
    email: str = "",
    api_key: Optional[str] = None,
    flank: int = 400,
    product_size_range: Tuple[int, int] = (100, 500),
    tm_min: float = 57.0,
    tm_opt: float = 62.0,
    tm_max: float = 67.0,
    gc_min: float = 40.0,
    gc_max: float = 65.0,
    user_transcript: str = "",
    cdna_label_override: str = "",
) -> Tuple[List[Dict], str]:
    """
    Main workflow for a single variant.
    Returns (rows, skip_reason). rows is a list of primer dicts, skip_reason is empty if successful.
    """
    gene          = gene.strip().upper()
    pos_input     = position_str.strip()
    skip_reason   = ""
    # If a cDNA label was pre-supplied (both columns were filled in the CSV),
    # we skip cDNA→genomic conversion and use the label purely for output.
    cdna_label    = cdna_label_override.strip()
    transcript_id = ""        # MANE Select transcript, filled by get_exons_ensembl

    # ── 1. Exon annotation (always needed; also gives MANE transcript ID) ─────
    try:
        exons, transcript_id = get_exons_ensembl(gene, build)
    except Exception as exc:
        log.warning("Ensembl lookup failed for %s: %s. Exon label = '?'.", gene, exc)
        exons         = []
        transcript_id = ""

    # ── 2. Resolve position ────────────────────────────────────────────────────
    # If cdna_label_override was set, position_str is already genomic — skip conversion.
    if cdna_label_override:
        # Both columns were filled: use genomic position directly.
        try:
            chrom, var_pos = parse_position(pos_input)
        except ValueError as exc:
            reason = f"Position parse error: {exc}"
            log.error("%s for %s", reason, gene)
            return [], reason
        log.info("Processing %s @ %s:%d (genomic; cDNA label=%s, build=%s)",
                 gene, chrom, var_pos, cdna_label, build)
    elif is_cdna_position(pos_input):
        # Only cDNA provided — convert to genomic via Ensembl map/cdna
        cdna_num   = parse_cdna_pos(pos_input)

        # Which transcript to use for coordinate mapping?
        # Priority: user-supplied transcript > MANE Select
        mapping_transcript = user_transcript.strip() if user_transcript.strip() else transcript_id
        if not mapping_transcript:
            reason = "Cannot convert cDNA position: no transcript found"
            log.error("%s for %s", reason, gene)
            return [], reason
        if user_transcript.strip():
            log.info("Using user-specified transcript for cDNA mapping: %s", mapping_transcript)
        try:
            chrom, gen_start, _gen_end = cdna_to_genomic(mapping_transcript, cdna_num, build)
            var_pos = gen_start
        except Exception as exc:
            reason = f"cDNA→genomic conversion failed: {exc}"
            log.error("%s for %s %s", reason, gene, pos_input)
            return [], reason
        log.info("Processing %s @ %s (cDNA %s → chr%s:%d, build=%s)",
                 gene, pos_input, pos_input, chrom, var_pos, build)
    else:
        # Plain genomic coordinate
        try:
            chrom, var_pos = parse_position(pos_input)
        except ValueError as exc:
            reason = f"Position parse error: {exc}"
            log.error("%s for %s", reason, gene)
            return [], reason
        log.info("Processing %s @ %s:%d (build=%s)", gene, chrom, var_pos, build)

    # ── 3. Determine exon number ────────────────────────────────────────────────
    exon_num = find_exon(var_pos, exons) if exons else "?"

    # ── 4. Fetch flanking sequence ─────────────────────────────────────────────
    seq_start = max(1, var_pos - flank)
    seq_end   = var_pos + flank
    try:
        sequence = fetch_sequence(chrom, seq_start, seq_end, build, email, api_key)
    except Exception as exc:
        reason = f"Failed to fetch sequence: {exc}"
        log.error("%s for %s %s", reason, gene, pos_input)
        return [], reason

    target_offset = var_pos - seq_start    # 0-based offset of variant inside window

    # ── 5. Design primers ──────────────────────────────────────────────────────
    candidates = design_primers(
        sequence           = sequence,
        target_offset      = target_offset,
        product_size_range = product_size_range,
        tm_min  = tm_min,
        tm_opt  = tm_opt,
        tm_max  = tm_max,
        gc_min  = gc_min,
        gc_max  = gc_max,
    )

    if not candidates:
        reason = "No primer candidates found meeting basic criteria (Tm/GC/Size)"
        log.warning("%s for %s @ %s.", reason, gene, pos_input)
        return [], reason

    # ── 6. Validation (In-Silico PCR) Loop ─────────────────────────────────────
    # We test candidates until we find one that is unique (1 amplicon total)
    pair = None
    amplicons = 0
    
    for i, cand in enumerate(candidates):
        log.info("Checking candidate %d/%d for uniqueness...", i + 1, len(candidates))
        amps = check_ucsc_insilico_pcr(cand["left_seq"], cand["right_seq"], build)
        if amps == 1:
            log.info("Candidate %d is unique! Selecting this pair.", i + 1)
            pair = cand
            amplicons = amps
            break
        elif amps > 1:
            log.info("Candidate %d is non-unique (%d amplicons).", i + 1, amps)
        elif amps == 0:
            log.info("Candidate %d returned no hits (check primer mapping).", i + 1)
        else:
            log.info("Candidate %d check failed.", i + 1)

    if pair is None:
        reason = f"No unique primer pair found among {len(candidates)} candidates"
        log.warning("%s for %s @ %s. Skipping variant.", reason, gene, pos_input)
        return [], reason

    # ── 7. Build output rows ───────────────────────────────────────────────────
    exon_label = exon_num if exon_num != "?" else "1"
    fwd_name   = f"{gene}_{exon_label}F"
    rev_name   = f"{gene}_{exon_label}R"

    # Compute genomic coordinates of each primer from their position in the window
    chrom_bare = chrom.lstrip("chr")
    def _primer_coords(win_start: int, seq_offset_start: int, seq_offset_end: int) -> str:
        """Convert 0-based window offsets to genomic coordinate string."""
        g_start = win_start + seq_offset_start
        g_end   = win_start + seq_offset_end
        return f"chr{chrom_bare}:{g_start}-{g_end}"

    fwd_coords = _primer_coords(seq_start, pair.get("left_start", 0), pair.get("left_end", 0))
    rev_coords = _primer_coords(seq_start, pair.get("right_start", 0), pair.get("right_end", 0))

    base_row = {
        "Gene":             gene,
        "Exon":             exon_num,
        "Genomic_Position": f"chr{chrom_bare}:{var_pos}",
        "cDNA_Position":    cdna_label,
        "Transcript":       user_transcript.strip() or transcript_id,
        "Build":            build,
        "ProductSize":      pair["product_size"],
        "Off_Target_Amplicons": amplicons - 1 if amplicons > 0 else (0 if amplicons == 0 else "Error"),
    }
    rows = [
        {"Name": fwd_name, "Sequence": pair["left_seq"],
         "TM": pair["left_tm"], "GC_percent": pair["left_gc"],
         "Primer_Coords": fwd_coords, **base_row},
        {"Name": rev_name, "Sequence": pair["right_seq"],
         "TM": pair["right_tm"], "GC_percent": pair["right_gc"],
         "Primer_Coords": rev_coords, **base_row},
    ]
    return rows, ""


# ─────────────────────────────────────────────────────────────────────────────
# 6.  CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Design PCR primers from variant genomic positions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Input CSV. Required columns: 'gene'. "
                        "Supply either 'cdna_position' (e.g. c.617G>A) or "
                        "'position' (genomic, e.g. chr2:157774114). "
                        "Optional column: 'transcript' (NM_ or ENST ID used for cDNA mapping).")
    p.add_argument("--output", "-o", default="primers_output.csv",
                   help="Output CSV file path.")
    p.add_argument("--email",  "-e", required=True,
                   help="Email address for NCBI Entrez (required by NCBI).")
    p.add_argument("--api-key", default=None,
                   help="NCBI API key (optional; allows >3 requests/sec).")
    p.add_argument("--build", choices=["hg38","GRCh38","hg19","GRCh37"],
                   default="hg38",
                   help="Reference genome build.")
    p.add_argument("--flank", type=int, default=400,
                   help="Flanking bases on each side of the variant.")
    p.add_argument("--product-size", default="300-500",
                   help="Acceptable PCR product size range, e.g. 100-500.")
    p.add_argument("--tm-min", type=float, default=57.0, help="Minimum primer Tm (°C).")
    p.add_argument("--tm-opt", type=float, default=62.0, help="Optimal primer Tm (°C).")
    p.add_argument("--tm-max", type=float, default=67.0, help="Maximum primer Tm (°C).")
    p.add_argument("--gc-min", type=float, default=40.0, help="Minimum GC%.")
    p.add_argument("--gc-max", type=float, default=65.0, help="Maximum GC%.")
    return p.parse_args()


def main():
    args = parse_args()

    # Parse product size range
    try:
        lo, hi = (int(x) for x in args.product_size.split("-"))
        product_size_range = (lo, hi)
    except ValueError:
        sys.exit(f"[ERROR] --product-size must be in the form MIN-MAX, e.g. 100-500.")

    # Read input CSV
    try:
        df_input = pd.read_csv(args.input)
    except FileNotFoundError:
        sys.exit(f"[ERROR] Input file not found: {args.input}")

    # Normalize column names
    df_input.columns = [c.strip().lower().replace(" ", "_") for c in df_input.columns]
    if "gene" not in df_input.columns:
        sys.exit("[ERROR] Input CSV must contain a 'gene' column.")
    has_cdna = "cdna_position" in df_input.columns
    has_geom = "position" in df_input.columns
    if not has_cdna and not has_geom:
        sys.exit("[ERROR] Input CSV must contain a 'cdna_position' column "
                 "(e.g. c.617G>A) and/or a 'position' column (e.g. chr2:157774114).")

    all_rows = []
    skipped_variants = []
    n_total  = len(df_input)

    for idx, row in df_input.iterrows():
        gene = str(row["gene"]).strip()

        # ── Determine which position to use ─────────────────────────────────
        def _col(name: str) -> str:
            val = str(row.get(name, "")).strip()
            return "" if val.lower() in ("", "nan", "none") else val

        cdna_pos_raw = _col("cdna_position") if has_cdna else ""
        geom_pos_raw = _col("position")      if has_geom else ""
        user_tx      = _col("transcript") if "transcript" in df_input.columns else ""

        if cdna_pos_raw and geom_pos_raw:
            position          = geom_pos_raw
            cdna_label_arg    = cdna_pos_raw
        elif cdna_pos_raw:
            position       = cdna_pos_raw
            cdna_label_arg = ""
        else:
            position       = geom_pos_raw
            cdna_label_arg = ""

        log.info("── Variant %d/%d: %s %s ──", idx + 1, n_total, gene, position)

        try:
            rows, reason = process_variant(
                gene                = gene,
                position_str        = position,
                build               = args.build,
                email               = args.email,
                api_key             = args.api_key,
                flank               = args.flank,
                product_size_range  = product_size_range,
                tm_min              = args.tm_min,
                tm_opt              = args.tm_opt,
                tm_max              = args.tm_max,
                gc_min              = args.gc_min,
                gc_max              = args.gc_max,
                user_transcript     = user_tx,
                cdna_label_override = cdna_label_arg,
            )
            if not rows and reason:
                skipped_variants.append({
                    "Gene": gene,
                    "Target_Position": position,
                    "Skip_Reason": reason
                })
            else:
                all_rows.extend(rows)
        except Exception as exc:
            log.error("Unexpected error for %s %s: %s", gene, position, exc)
            skipped_variants.append({
                "Gene": gene, 
                "Target_Position": position, 
                "Skip_Reason": f"Crash: {exc}"
            })

    # Write output CSV
    if all_rows:
        out_df = pd.DataFrame(all_rows, columns=[
            "Name", "Sequence", "TM", "GC_percent",
            "ProductSize", "Off_Target_Amplicons", "Primer_Coords",
            "Gene", "Exon", "Genomic_Position", "cDNA_Position", "Transcript", "Build",
        ])
        out_df.to_csv(args.output, index=False)
        log.info("Saved %d primer records to %s", len(all_rows), args.output)
        print(f"\nDone! Primers written to: {args.output}")
    else:
        log.warning("No primers were designed.")
        print("\nNo primers were designed. Check warnings above.")

    # Write skipped variants CSV
    if skipped_variants:
        skipped_file = args.output.replace(".csv", "_skipped.csv")
        if not skipped_file.endswith("_skipped.csv"):
            skipped_file += "_skipped.csv"
        pd.DataFrame(skipped_variants).to_csv(skipped_file, index=False)
        log.info("Saved %d skipped variants to %s", len(skipped_variants), skipped_file)
        print(f"Skipped variants recorded in: {skipped_file}")


if __name__ == "__main__":
    main()
