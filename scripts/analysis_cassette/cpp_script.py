from pathlib import Path
import time
from typing import Iterator
import pandas as pd
import os
import click
import networkit as nk

def find_experiment_root(start_path: str = None) -> Path:
    """
    Walks up from start_path (or current working directory if None) until it finds
    a directory that contains any of the marker entries: 'bak', 'errors', 'output',
    'scripts'. Returns the Path to that directory, or None if no such directory is
    found before reaching the filesystem root.
    """
    if start_path is None:
        
        start_path = Path(__file__).resolve().parent

    fixed_markers = {"bak", "errors", "output", "scripts"}

    for directory in [start_path] + list(start_path.parents):
        try:
            entries = {entry.name for entry in directory.iterdir()}
        except PermissionError:
            continue

        if fixed_markers.intersection(entries):
            return directory

    return None

def count_of_nodes_seeds_and_agents(aux_path: Path) -> dict:
    """
    Reads the given output.aux CSV file and returns a dictionary with:
      - "seeds": count of rows where type == "seed"
      - "agents": count of rows where type == "agent"
      - "total": sum of seeds and agents
    """
    df = pd.read_csv(aux_path, low_memory=False)
    seed_count = (df["type"] == "seed").sum()
    agent_count = (df["type"] == "agent").sum()
    return {
        "seeds": int(seed_count),
        "agents": int(agent_count),
        "total": int(seed_count + agent_count)
    }

def extract_num_cycles_from_err(err_path: Path) -> int:
    """
    Reads the abm.err file, splits each line by whitespace, 
    finds '--num-cycles', then returns the next token as int.
    """
    with open(err_path, "r") as f:
        for line in f:
            if "--num-cycles" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "--num-cycles" and i + 1 < len(parts):
                        return int(parts[i + 1])
    raise ValueError(f"'--num-cycles' not found in {err_path}")

def extract_growth_rate_from_err(err_path: Path) -> int:
    """
    Reads the abm.err file, splits each line by whitespace, finds '--growth-rate',
    then returns the next token as a float.

    Raises:
        FileNotFoundError: If err_path does not exist.
        ValueError: If '--growth-rate' flag is not found in the file.
    """
    with open(err_path, "r") as f:
        for line in f:
            if "--growth-rate" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "--growth-rate" and i + 1 < len(parts):
                        return float(parts[i + 1])
    raise ValueError(f"'--growth-rate' not found in {err_path}")

def extract_agent_environment(abm_err: Path) -> str:
    """
    Reads the abm.err file and determines the agent environment based on the
    presence of specific flags. Returns "sa" or "ra" based on the flags found.
    """
    background_flags = ["--recency-weight", "--fitness-weight", "--preferential-weight"]
    with open(abm_err, "r") as f:
        for line in f:
            for flag in background_flags:
                if flag in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == flag and i + 1 < len(parts):
                            if parts[i+1] == "0.33":
                                return "sa"
    return "ra"

def count_number_of_superstars(aux_path: Path) -> int:
    """
    Reads the given output.aux CSV file and returns the count of rows where
    fit_peak_value > 1000.
    """
    df = pd.read_csv(aux_path, low_memory=False)
    superstar_count = (df["fit_peak_value"] > 1000).sum()
    return superstar_count

def in_degree_distribution(aux_path: Path) -> dict:
    """
    Reads output.aux and returns a dict with five-number summaries of in_degree
    for each type (“seed” and “agent”):
      {
        "seed":  {"min": ..., "Q1": ..., "median": ..., "Q3": ..., "max": ...},
        "agent": {"min": ..., "Q1": ..., "median": ..., "Q3": ..., "max": ...}
      }
    """
    df = pd.read_csv(aux_path, low_memory=False)
    result = {}
    for node_type, group in df.groupby("type"):
        s = group["in_degree"]
        result[node_type] = {
            "min": int(s.min()),
            "Q1":  float(s.quantile(0.25)),
            "median": float(s.quantile(0.50)),
            "Q3":  float(s.quantile(0.75)),
            "max": int(s.max())
        }
    return result

def extract_out_degree_distribution(abm_err: Path) -> str:
    """
    Reads the abm.err file, finds the '--out-degree-bag' flag, and returns the
    filename part of the provided path (the basename after the last '/').

    Raises:
        FileNotFoundError: If abm_err does not exist.
        ValueError: If '--out-degree-bag' flag is not found in the file.
    """
    with open(abm_err, "r") as f:
        for line in f:
            if "--out-degree-bag" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "--out-degree-bag":
                        path = parts[i + 1]
                        path = path.split("/")[-1]
                        return path
    raise ValueError(f"'--out-degree-bag' not found in {abm_err}")

def extract_run_time_from_err(err_path: Path) -> str:
    """
    Reads the given abm.err file and returns the wall-clock run time as a string
    in "h:mm:ss" or "m:ss" format, extracted from the line starting with 
    'Elapsed (wall clock) time'.
    """
    with open(err_path, "r") as f:
        for line in f:
            if "Elapsed (wall clock) time" in line:
                parts = line.split(": ", 1)
                if len(parts) == 2:
                    return parts[1].strip()
    raise ValueError(f"'Elapsed (wall clock) time' not found in {err_path}")

def extract_memory_requested_from_sbatch(sbatch_path: Path) -> str:
    """
    Reads submit_job.sbatch and returns the memory requested (e.g. "200GB"),
    found on the line starting with '#SBATCH --mem='.
    """
    with open(sbatch_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#SBATCH --mem="):
                return line.split("=", 1)[1]
    raise ValueError(f"'--mem=' not found in {sbatch_path}")

def extract_memory_used_from_err(err_path: Path) -> int:
    """
    Reads the abm.err file and returns the maximum resident set size as a string
    in gigabytes (e.g. "7.49GB").

    Raises:
        FileNotFoundError: If err_path does not exist.
        ValueError: If 'Maximum resident set size' is not found in the file.
    """
    with open(err_path, "r") as f:
        for line in f:
            if "Maximum resident set size" in line:
                parts = line.split(":")
                if len(parts) >= 2:
                    kb = int(parts[-1].strip())
                    gb = kb / (1024 ** 2)
                    gb = str(round(gb, 2))
                    return gb + "GB"
    raise ValueError(f"'Maximum resident set size' not found in {err_path}")

def extract_input_files_from_amb_err(abm_err: Path) -> dict:
    """
    Reads the abm.err file and extracts the input nodelist and edgelist files.
    Returns a dictionary with keys 'nodelist' and 'edgelist'.
    """
    input_files = {"nodelist": None, "edgelist": None}
    with open(abm_err, "r") as f:
        for line in f:
            if "--edgelist" in line and "--nodelist" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "--nodelist" and i + 1 < len(parts):
                        input_files["nodelist"] = parts[i + 1].split("/")[-1]
                    elif part == "--edgelist" and i + 1 < len(parts):
                        input_files["edgelist"] = parts[i + 1].split("/")[-1]
    if not input_files["nodelist"] or not input_files["edgelist"]:
        raise ValueError(f"'--nodelist' or '--edgelist' not found in {abm_err}")
    return input_files

def extract_alpha_from_abm_err(abm_err: Path) -> float:
    """
    Reads amb_err and returns the alpha value as a float (from '--alpha' flag).;
    if no alpha line is found, returns the string 'Random'.

    Raises:
        FileNotFoundError: If abm_err does not exist.
    """
    with open(abm_err, "r") as f:
        for line in f:
            if "--alpha" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "--alpha" and i + 1 < len(parts):
                        raw = parts[i + 1].strip('"').strip("'")
                        return float(raw)
    return "Random"

def extract_date_of_last_modification(abm_err: Path) -> str:
    """
    Returns the last modification date of the abm.err file as a string.
    The date is formatted as a string in the local time zone.
    """
    #from: https://www.geeksforgeeks.org/python-os-path-getctime-method/
    c_time = os.path.getmtime(abm_err)
    local_time = time.ctime(c_time)
    return local_time

def Extract_clustering_coefficient(edge_list: Path) -> dict:
    """
    Reads the edge list file and calculates the clustering coefficient
    using Networkit. Returns a dictionary with the GCC and ALCC values.
    Note: This function assumes the edge list is in CSV format with columns
    '#source' and 'target as the edgelist generated by the abm script'.
    """
    edges = pd.read_csv(edge_list)

    # first_col = edges.columns[0]
    # second_col = edges.columns[1]
    # # Rename columns to acept any edgelist format:
    # edges = edges.rename(columns={first_col: '#source', second_col: 'target'})
    
    min_source = edges["#source"].min()
    min_target = edges["target"].min()

    min_id = min_source if min_source < min_target else min_target

    reader = nk.graphio.EdgeListReader(
            separator=",",
            firstNode=min_id,
            commentPrefix="#",
            continuous=False,
            directed=False
        )
    G = reader.read(str(edge_list))
    cluster_coeff = nk.globals.ClusteringCoefficient()

    return {"gcc":cluster_coeff.exactGlobal(G),
            "alcc": cluster_coeff.sequentialAvgLocal(G)}

def generate_metrics(abm_err: Path, output_aux: Path, submit_job_sbatch: Path, edge_list: Path) -> dict:
    """
    Main function to extract various metrics from the provided files.
    Returns a dictionary with the following and generates a csv file.
    """
    experiment_root = find_experiment_root(Path(abm_err).resolve().parent)
    dict_in_degree_distribution = in_degree_distribution(output_aux)
    dict_caunt_of_nodes = count_of_nodes_seeds_and_agents(output_aux)
    clustering_coeff = Extract_clustering_coefficient(edge_list)
    metrics = {
        "name_of_experiment": experiment_root.name,
        "date_of_experiment": extract_date_of_last_modification(abm_err),
        "input_nodelist": extract_input_files_from_amb_err(abm_err).get("nodelist"),
        "input_edgelist": extract_input_files_from_amb_err(abm_err).get("edgelist"),
        "count_of_seeds": dict_caunt_of_nodes.get("seeds"),
        "count_of_agents": dict_caunt_of_nodes.get("agents"),
        "count_of_nodes": dict_caunt_of_nodes.get("total"),
        "num_cycles": extract_num_cycles_from_err(abm_err),
        "growth_rate": extract_growth_rate_from_err(abm_err),
        "agent_environment": extract_agent_environment(abm_err),
        "alpha": extract_alpha_from_abm_err(abm_err),
        "GCC": clustering_coeff.get("gcc"),
        "ALCC": clustering_coeff.get("alcc"),
        "count_of_superstars": count_number_of_superstars(output_aux),
        "seed_indeg_min": dict_in_degree_distribution.get("seed").get("min"),
        "seed_indeg_Q1": dict_in_degree_distribution.get("seed").get("Q1"),
        "seed_indeg_median": dict_in_degree_distribution.get("seed").get("median"),
        "seed_indeg_Q3": dict_in_degree_distribution.get("seed").get("Q3"),
        "seed_indeg_max": dict_in_degree_distribution.get("seed").get("max"),
        "agent_indeg_min": dict_in_degree_distribution.get("agent").get("min"),
        "agent_indeg_Q1": dict_in_degree_distribution.get("agent").get("Q1"),
        "agent_indeg_median": dict_in_degree_distribution.get("agent").get("median"),
        "agent_indeg_Q3": dict_in_degree_distribution.get("agent").get("Q3"),
        "agent_indeg_max": dict_in_degree_distribution.get("agent").get("max"),
        "out_degree_distribution": extract_out_degree_distribution(abm_err),
        "run_time": extract_run_time_from_err(abm_err),
        "memory_requested": extract_memory_requested_from_sbatch(submit_job_sbatch),
        "memory_used": extract_memory_used_from_err(abm_err),
    }
    return metrics


def extract_metrics_V2(abm_err: Path, output_aux: Path, submit_job_sbatch: Path, 
                edge_list: Path, df_metrics: pd.DataFrame = None) -> pd.DataFrame:
    """
    Runs the second version of the cassette, extractig metrics from the files provided by function:
    main_V2 and returns a pd.DataFrame with the metrics.
    """
    metrics = generate_metrics(abm_err, output_aux, submit_job_sbatch, edge_list)

    if df_metrics is None:
        df_metrics = pd.DataFrame([metrics])
        return df_metrics

    name = metrics["name_of_experiment"]
    df_metrics = df_metrics[df_metrics["name_of_experiment"] != name]

    df_metrics = df_metrics.reset_index(drop=True)

    df_metrics.loc[len(df_metrics)] = metrics

    return df_metrics

def discover_experiments(root: Path) -> Iterator[Path]:
    """
    Recursively walk the directory tree under `root` and yield every folder
    that looks like a complete ABM experiment, i.e. it contains all three
    subdirectories AND each of those contains the required file.

    An “experiment folder” must have:
      • errors/abm.err
      • scripts/submit_job.sbatch
      • output/output.aux

    Yields:
        Path to each experiment directory found anywhere under `root`.
    """
    root = Path(root)
    for dirpath, dirnames, filenames in os.walk(root):
        d = Path(dirpath)

        # check that all three sub-folders exist
        if {"errors", "scripts", "output"}.issubset(dirnames):
            # now verify the four required files inside them
            err_file   = d / "errors"  / "abm.err"
            sbatch_file= d / "scripts" / "submit_job.sbatch"
            aux_file   = d / "output"  / "output.aux"
            edge_list  = d / "output"  / "output.edgelist"
            if err_file.exists() and edge_list.exists()\
               and sbatch_file.exists() and aux_file.exists():
                yield d
            else:
                print(f"Skipping {d.name} not indentified as a complete experiment folder")
                print(f"{d.name} Missing one or more files: 'abm.err', 'edge_list' " \
                      "'submit_job.sbatch', 'output.aux'\n")

@click.command(help=
    """
    Recursively evaluate ABM experiment folders and update a master CSV.

    Usage:

    cpp_script.py --abm-outputs .../abm_outputs [--csv-path PATH]

    Options:

    --abm-outputs -> PATH Path to the abm_outputs directory (obligatory).

    --csv-path -> PATH CSV file to read/write metrics (default: .../abm_outputs/experiment_metrics.csv).
    """
    )

@click.option(
    "--abm-outputs", "abm_outputs",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help="-> Path to the abm_outputs directory containing experiment folders" \
    ""
)
@click.option(
    "--csv-path", "csv_path",
    required=False,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="-> Path to the master CSV file for metrics. If not provided," \
    "it will default to 'experiment_metrics.csv' in the abm_outputs directory." \
    ""
)
def cli(abm_outputs: Path, csv_path: Path = None):

    abm_outputs = Path(abm_outputs)

    metrics_df = None
    if csv_path:
        csv_path = Path(csv_path)
    else:
        csv_path = abm_outputs / "experiment_metrics.csv"

    if csv_path.exists():
        metrics_df = pd.read_csv(csv_path)
        processed = set(metrics_df["name_of_experiment"])
    else:
        metrics_df = None
        processed = set()

    evaluated_exp = []
    for exp_folder in discover_experiments(abm_outputs):

        if exp_folder.name in processed:
            print(f"Skipping already processed experiment: {exp_folder.name}\n")
            continue

        print(f"Evaluating experiment: {exp_folder.name}\n")

        abm_err           = exp_folder / "errors"  / "abm.err"
        output_aux        = exp_folder / "output"  / "output.aux"
        submit_job_sbatch = exp_folder / "scripts" / "submit_job.sbatch"
        edge_list         = exp_folder / "output" / "output.edgelist"

        metrics_df = extract_metrics_V2(
            abm_err=abm_err,
            output_aux=output_aux,
            submit_job_sbatch=submit_job_sbatch,
            edge_list=edge_list,
            df_metrics=metrics_df
        )
        evaluated_exp.append(exp_folder.name)
    print(f"Evaluated experiments: {evaluated_exp}\n")

    if evaluated_exp:
        print(f"{csv_path.name} file saved to {csv_path.parent}")
        metrics_df.to_csv(csv_path, index=False)
    else:
        print("No new experiments to evaluate or no experiments found. No changes made to the CSV file.")

if __name__ == "__main__":
    cli()