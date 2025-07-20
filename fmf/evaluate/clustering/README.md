## Leiden Clustering Pipeline

**Description**  
Generates clustering information for ABM experiments using the Leiden algorithm with a specified resolution parameter. Optionally, computes additional metrics such as global and average local clustering coefficients (GCC and ALCC), node coverage, and structural properties of each cluster. If superstars are present, their cluster statistics are also extracted.

---

### Output Files

Each valid experiment folder produces:

- `leidenClustering.csv` — Leiden clustering file (node id, cluster id).
- `clustering_metrics.csv` — Metrics per cluster (size, edges, mincut, etc.).
- `summary.csv` — One summary row per experiment (written to path given via `--output`).
- `superstars_gamma_X.csv` — Rows for each superstar cluster (only if superstars exist).

---

### Output Columns

#### `summary.csv`

| Column Name       | Description                                                             |
|-------------------|-------------------------------------------------------------------------|
| `exp_name`        | Name of the experiment folder                                            |
| `node_coverage`   | Fraction of nodes assigned to non-singleton clusters                    |
| `min`             | Minimum cluster size (excluding singletons)                             |
| `q1`              | First quartile of cluster sizes                                          |
| `median`          | Median cluster size                                                     |
| `q3`              | Third quartile of cluster sizes                                          |
| `p90`             | 90th percentile of cluster sizes                                         |
| `max`             | Maximum cluster size                                                    |
| `singletons`      | Number of singleton clusters (size = 1)                                 |
| `gcc`             | Global Clustering Coefficient (Networkit)                               |
| `alcc`            | Average Local Clustering Coefficient (Networkit)                        |

#### `superstars_gamma_X.csv` (if applicable)

| Column Name          | Description                                                 |
|----------------------|-------------------------------------------------------------|
| `ss_id`              | Node ID of the superstar                                    |
| `ss_fit`             | Peak fitness value of the superstar                         |
| `cluster_id`         | ID of the cluster the superstar belongs to                  |
| `size`               | Size of the cluster                                         |
| `intra_edges`        | Number of internal edges within the cluster                 |
| `boundary_edges`     | Number of edges connecting to nodes outside the cluster     |
| `normalized_density` | Normalized density of the cluster                           |
| `mincut`             | Mincut value computed on the induced subgraph               |
| `experiment`         | Name of the experiment folder                               |
| `gamma`              | Resolution parameter used in Leiden clustering              |

---

### Usage on Valhalla

```bash
module load python3/3.11.0
source venv/bin/activate
pip install -r requirements.txt

python3 leiden_pipeline.py \
  --experiments-dirs /path/to/exp1 --experiments-dirs /path/to/exp2 ... \
  --output /path/to/summary.csv \
  --ss-out /path/to/superstars.csv \
  --full-pipeline \
  --gamma 0.001
```

---

### Command-Line Options

| Option               | Description                                                                 |
|----------------------|-----------------------------------------------------------------------------|
| `--experiments-dirs` | **(Required)** List of paths to experiment folders                          |
| `--output`           | **(Required)** Path to save the experiment summary CSV (`summary.csv`)      |
| `--ss-out`           | **(Required)** Path to save superstar cluster CSV                           |
| `--full-pipeline`    | Run Leiden clustering + full metrics pipeline (optional but recommended)    |
| `--gamma`            | Resolution parameter for Leiden clustering (default: `0.001`)               |
| `--clustering-name`  | Name of clustering file to be generated (default: `leidenClustering.csv`)   |
| `--metrics-name`     | Name of clustering metrics file (default: `clustering_metrics.csv`)         |
