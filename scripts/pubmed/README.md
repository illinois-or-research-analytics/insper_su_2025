# PubMed Sampling and Author Country Extraction Scripts

This repository contains two Python scripts designed to interact with PubMed through the NCBI Entrez API (via Biopython) for bibliometric studies:

1. **`sampling_pubmed.py`** - Samples PubMed journal articles randomly within a range of PMIDs, counting their reference out-degree (number of references cited), producing data suitable for ABM experiments.

2. **`check_author_pubmed.py`** - Retrieves PubMed articles related to exosomes and extracts country affiliations of authors along with publication year, to explore the geographic spread of research.

---

# Script 1: sampling_pubmed.py

## What it does
- Randomly samples PMIDs between two given IDs.
- Fetches each article, checks if it's a journal article with at least 5 references.
- Saves a CSV file listing each PMID and its number of references (out-degree).

## How to run
Run with click options:
```bash
python sampling_pubmed.py \
  --id_2020 31900000 \
  --id_2025 37500000 \
  --sample_size 1000 \
  --output_dir ../outputs \
  --filename pubmed_references
```
- `id_2020`: PMID roughly from year 2020  
- `id_2025`: PMID roughly from year 2025  
- `sample_size`: how many valid articles to sample  
- `output_dir`: where to save the output  
- `filename`: CSV file base name (it will append `_1000.csv` for 1000 samples)

## Output
A CSV in your `outputs/` folder like:

#node_id    out_degree  
34500012    12  
35200134    6  
...         ...

---

# Script 2: check_author_pubmed.py

## What it does
- Searches PubMed for articles tagged with `"exosomes"[MeSH Terms]` from 1983 to 2025.
- Fetches author affiliations and publication years.
- Parses the last segment after the last comma in the affiliation string as the country.
- Saves a CSV listing PMID, year, and country.

## How to run

python check_author_pubmed.py

It automatically:
- fetches up to 100,000 articles,
- extracts affiliation countries,
- and writes to `../outputs/exosome_affiliations.csv`.

## Output
A CSV in your `outputs/` folder like:

PMID    Year    Country  
33410001    2021    USA  
33410001    2021    Germany  
33410234    2020    Unknown  
...         ...     ...

---

# Notes & Recommendations

**Entrez Usage:**  
Make sure to personalize `Entrez.email` in both scripts. This is required by NCBI.

**Outputs:**  
Place your `outputs/` folder outside of your scripts to keep things organized.

**Scaling:**  
For large samples (10,000+), consider running with `nohup` to keep processes alive after closing your terminal.


---
