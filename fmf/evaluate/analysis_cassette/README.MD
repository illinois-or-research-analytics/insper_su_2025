## Analysis Cassette

**Description**  
Recursively evaluates ABM experiments, merging metrics from each experiment into a single CSV file. Each row corresponds to one experiment.

---

### Output Columns

| Column Name                      | Description                                                             |
|----------------------------------|-------------------------------------------------------------------------|
| `name_of_experiment`             | Identifier of the experiment (folder name)                              |
| `date_of_experiment`             | Date when the experiment was run (YYYY-MM-DD)                           |
| `input_nodelist`                 | Path to the input nodelist file                                         |
| `input_edgelist`                 | Path to the input edgelist file                                         |
| `count_of_seeds`                 | Number of seed nodes used                                               |
| `count_of_agents`                | Total number of agents simulated                                        |
| `count_of_nodes`                 | Total number of nodes in the network                                    |
| `num_cycles`                     | Number of simulation cycles                                             |
| `growth_rate`                    | Growth rate parameter                                                   |
| `agent_environment`              | Agent environment setting (e.g., RA or SA)                              |
| `alpha`                          | Alpha parameter                                                         |
| `count_of_superstars`            | Number of superstars introduced                                         |
| `seeds_in_degree_distribution`   | In‑degree distribution metrics for seed nodes                           |
| `agents_in_degree_distribution`  | In‑degree distribution metrics for all agents                           |
| `out_degree_bag`                 | Out‑degree distribution metrics                                         |
| `run_time`                       | Total runtime of the experiment (e.g., in seconds)                      |
| `memory_requested`               | Amount of memory requested (e.g., in MB)                                |
| `memory_used`                    | Actual memory used during execution (e.g., in MB)                       |

---

### Usage on Valhalla

```bash

module load python3/3.11.0
source venv/bin/activate
pip install -r requirements.txt

python3 cpp_script.py \
  --abm-outputs /path/to/abm_outputs \
  [--csv-path    /path/to/output.csv] \
  [--full-pipe]

--abm-outputs (required): Root folder containing all experiment subfolders.

--csv-path (optional): Where to write the merged CSV (defaults to abm-outputs).

--full-pipe (optional): Run GCC & ALCC calculations (these take longer).

Example Folder Structure

abm_outputs/
├── exp1/
├── exp2/
├── random.csv
└── other_experiments/
    ├── exp3/
    ├── exp4/
    └── random2.csv

In this layout, exp1, exp2, exp3 and exp4 will be evaluated; files like random.csv & random2.csv will be ignored.

Experiment Recognition Criteria
To be recognized as a valid experiment, each subfolder must contain:

errors/
└── abm.err

scripts/
└── submit_job.sbatch

output/
├── output.aux
└── output.edgelist
