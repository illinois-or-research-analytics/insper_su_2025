# Co‑Citation Strength

Compute co‑citation strength from an edgelist CSV file. The co‑citation strength between nodes *i* and *j* is the number of source nodes that cite both *i* and *j*.

## Output Columns

| Column Name       | Description                                                      |
|-------------------|------------------------------------------------------------------|
| `source1`         | ID of the first target node (*i*)                                |
| `source2`         | ID of the second target node (*j*)                               |
| `Co-strength`     | Co‑citation strength               |
| `type1`           | `type` metadata for `source1` (if `--output_aux_path` provided)  |
| `year1`           | `year` metadata for `source1` (if `--output_aux_path` provided)  |
| `type2`           | `type` metadata for `source2` (if `--output_aux_path` provided)  |
| `year2`           | `year` metadata for `source2` (if `--output_aux_path` provided)  |

---

## Options and Usage

```bash
module load python3/3.11.0
source venv/bin/activate
pip install -r requirements.txt

python co_citation.py \
  --edgelist_path         /path/to/edgelist.csv \
  --csv_path              /path/to/output.csv \
  [--output_aux_path      /path/to/output.aux] \
  [--in_degree_threshold  0] \
  [--n_nodes              None] \
  [--threshold-co-strength 0]

--edgelist_path (required): Path to the edgelist CSV file.

--csv_path (required): Path where the co‑citation CSV will be saved.

--output_aux_path (optional): output.aux CSV with columns node_id, type, year for metadata merging.

--in_degree_threshold (optional): Minimum in‑degree threshold for target nodes (default: 0).

--n_nodes (optional): Number of top target nodes by in‑degree to consider (default: all).

--threshold-co-strength (optional): Minimum co‑citation strength to include (default: 0).
