# Bibliographic Coupling

**Description**  
Compute bibliographic coupling from an edgelist CSV file. The bibliographic coupling between node *i* and *j* is the number of references they both cite in common.

**Formula**  
Let E be the set of directed edges (i, k) in the citation graph. For each source node i, define:

T_i = { k | (i, k) in E }

Then the bibliographic coupling strength between nodes i and j is:

Biblio-strength(i, j) = |T_i ∩ T_j|

---

## Output Columns

| Column Name        | Description                                                               |
|--------------------|---------------------------------------------------------------------------|
| `source1`          | ID of the first source node                                               |
| `source2`          | ID of the second source node                                              |
| `Biblio-strength`  | Between node *i* and *j*                                                  |
| `type1`            | Type for `source1` (if auxiliary file provided)                           |
| `year1`            | Year for `source1` (if auxiliary file provided)                           |
| `type2`            | Type for `source2` (if auxiliary file provided)                           |
| `year2`            | Year for `source2` (if auxiliary file provided)                           |

---

## Options and Usage

```bash
module load python3/3.11.0
source venv/bin/activate
pip install -r requirements.txt

python bibliographic.py \
  --edgelist_path /path/to/edgelist.csv \
  --csv_path      /path/to/output.csv \
  [--output_aux_path /path/to/aux.csv] \
  [--out_degree_threshold 0] \
  [--n_nodes None] \
  [--threshold-biblio-strength 0]

--edgelist_path (required): Path to the edgelist CSV file.

--csv_path (required): Path where the bibliographic coupling CSV will be saved.

--output_aux_path (optional): Auxiliary CSV with seed metadata (node_id, type, year).

--out_degree_threshold (optional): Minimum out-degree threshold for source nodes (default: 0).

--n_nodes (optional): Number of top source nodes by out-degree to consider (default: all).

--threshold-biblio-strength (optional): Minimum bibliographic strength to include (default: 0).
