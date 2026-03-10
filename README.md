# VarPrimer: Primer Design Tool

Design PCR primers from variant positions. Reads a CSV of genes and variant positions
(genomic or cDNA), fetches exon structure from Ensembl using the **MANE Select** transcript,
retrieves genomic sequence from **Ensembl REST**, designs primers with **primer3** (or a
built-in fallback), validates primer uniqueness via **UCSC In-Silico PCR**, and writes
results to a CSV file.

---

## Installation

```bash
pip install -r requirements.txt
```

> **Windows note:** `primer3-py` requires C++ Build Tools. If installation fails, the tool
> automatically falls back to a built-in pure-Python primer picker.

---

## Usage

```bash
python varprimer.py \
    --input  variants.csv \
    --output primers_output.csv \
    --email  your@email.com
```

### All options

| Argument | Default | Description |
|---|---|---|
| `--input` / `-i` | *(required)* | Input CSV file |
| `--output` / `-o` | `primers_output.csv` | Output CSV path |
| `--email` / `-e` | *(required)* | Email for NCBI Entrez (policy requirement) |
| `--api-key` | — | NCBI API key (raises rate limit to 10 req/sec) |
| `--build` | `hg38` | Genome build: `hg38` / `GRCh38` / `hg19` / `GRCh37` |
| `--flank` | `300` | Flanking bases on each side of the variant |
| `--product-size` | `100-500` | PCR product size range (bp) |
| `--tm-min/opt/max` | `57/62/67` | Primer melting temperature range (°C) |
| `--gc-min/max` | `40/65` | GC% range |

---

## Input CSV format

The CSV must have a `gene` column plus **at least one** of `cdna_position` or `position`:

| Column | Required? | Description |
|---|---|---|
| `gene` | ✅ | Gene symbol (e.g. `ACVR1`) |
| `cdna_position` | ✅ or `position` | cDNA notation (e.g. `c.617G>A`, `c.639_642del`) |
| `position` | ✅ or `cdna_position` | Genomic coordinate (e.g. `chr2:157774114`) |
| `transcript` | Optional | RefSeq (`NM_XXXXX`) or Ensembl (`ENST…`) ID — used only for cDNA-to-genomic conversion. Omit to use MANE Select automatically. |

### Priority rules

| What's filled | Behaviour |
|---|---|
| Both `cdna_position` **and** `position` | Genomic position is used as anchor; cDNA label is kept for output only |
| Only `cdna_position` | Converted to genomic coords via Ensembl `map/cdna` using the specified or MANE Select transcript |
| Only `position` | Used directly as the genomic anchor |

### Example

```csv
gene,cdna_position,transcript,position
ACVR1,c.617G>A,NM_001111067.4,chr2:157774114
MESD,c.639_642del,NM_015154.3,chr15:80979281
LAMC2,c.343C>T,NM_005562.3,chr1:183215527
BRCA1,,,chr17:43044295
```

- Commas in positions (e.g. `chr17:43,044,295`) are stripped automatically.
- The `chr` prefix is optional.

---

## Output CSV format

| Column | Description |
|---|---|
| `Name` | Primer name: `GENE_ExonF` / `GENE_ExonR` |
| `Sequence` | Primer sequence (5′ → 3′, positive strand) |
| `TM` | Melting temperature (°C) |
| `GC_percent` | GC content (%) |
| `ProductSize` | Expected PCR amplicon size (bp) |
| `Off_Target_Amplicons` | Number of unintended genomic amplicons found via UCSC In-Silico PCR (`0` = highly specific) |
| `Primer_Coords` | Exact chromosomal coordinates of this primer (e.g. `chr2:157773816-157773837`) |
| `Gene` | Gene symbol |
| `Exon` | Exon number (MANE Select transcript) |
| `Genomic_Position` | Genomic anchor position of the variant |
| `cDNA_Position` | cDNA notation from input (blank if only genomic was given) |
| `Transcript` | Transcript used (user-specified or MANE Select) |
| `Build` | Genome build |

### Naming convention

- `ACVR1_6F` — forward primer for ACVR1, variant in exon 6
- `ACVR1_6R` — reverse primer for ACVR1, variant in exon 6

If the variant is **intronic**, the nearest exon number is used and a warning is logged.

---

## Example output

```
Name,Sequence,TM,GC_percent,ProductSize,Primer_Coords,Gene,Exon,...
ACVR1_6F,TGGGCATTCTCTCATCATCCCA,61.2,50.0,496,chr2:157773816-157773837,ACVR1,6,...
ACVR1_6R,ACTAACAGGCCACGTGTCCC,62.1,60.0,496,chr2:157774292-157774311,ACVR1,6,...
```

## How exon numbers are assigned

1. Ensembl is queried for all transcripts of the gene.
2. The **MANE Select** transcript is preferred (clinical gold standard); falls back to Ensembl canonical, then longest.
3. Exons are numbered 1 → N in transcription order (5′ → 3′).
4. If the variant falls in an intron, the nearest exon number is assigned (with a log warning).

> **cDNA position accuracy:** if your variant report uses a specific RefSeq transcript
> (e.g. `NM_001105.5`), enter it in the `transcript` column. Omitting it will use the
> MANE Select transcript, which may number exons differently.

---

## Notes

- Requires an internet connection (Ensembl REST + NCBI Entrez APIs).
- Sequence is fetched from **Ensembl REST** (same source as Ensembl genome browser); NCBI is a fallback.
- NCBI rate-limits unauthenticated requests to 3 req/sec. Use `--api-key` for faster runs.
- TM is calculated with the **SantaLucia 1998 nearest-neighbour model** (or via primer3 if installed).
