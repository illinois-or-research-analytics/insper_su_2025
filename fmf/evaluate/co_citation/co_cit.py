import time
from typing import Dict, Sequence, Set, Tuple
import pandas as pd
from pathlib import Path
from itertools import combinations
import click
import csv

def preprocess_edgelist(
    edgelist_path: Path,
    # output_aux_path: Path = None,
    in_degree_threshold: int = 0,
    n_nodes: int = None
) -> Tuple[Sequence[int], Dict[int, Set[int]]]:
    """
    Returns:
      top_ids:     list of the top-n target IDs with in_degree > threshold
      sources_by_target: mapping target ID → set of source IDs
    """
    df = pd.read_csv(edgelist_path)
    first_col = df.columns[0]
    second_col = df.columns[1]

    df = df.rename(columns={first_col: '#source', second_col: 'target'})

    # if output_aux_path is not None:
    #     df_aux = pd.read_csv(output_aux_path)
    #     df_aux = df_aux[df_aux["type"] == "seed"]
    #     ids_seeds = df_aux["node_id"].value_counts().index.to_list()

    #     df = df[
    #         df['#source'].isin(ids_seeds) &
    #         df['target'].isin(ids_seeds)
    #     ]

    vc = df['target'].value_counts()
    filtered = vc[vc > in_degree_threshold]
    if n_nodes is None:
        top_ids = filtered.index.tolist()
    else:
        top_ids = filtered.head(n_nodes).index.tolist()
    sources_by_target = (
        df[df['target'].isin(top_ids)]
          .groupby('target')['#source']
          .agg(set)
          .to_dict()
    )
    return top_ids, sources_by_target

def merge_aux_co_strength(
    df_co: pd.DataFrame,
    df_aux: pd.DataFrame
):
    """
    Merges auxiliary metadata (type and year) into the co-citation DataFrame.
    """
    all_ids = pd.concat([df_co["source1"], df_co["source2"]])
    all_ids = all_ids.unique()
    df_aux = df_aux[df_aux["node_id"].isin(all_ids)]
    meta = df_aux[["node_id","type","year"]]

    df1 = (
        df_co
        .merge(meta, how="left", left_on="source1", right_on="node_id")
        .rename(columns={"type":"type1", "year":"year1"})
        .drop(columns=["node_id"])
    )

    final = (
        df1
        .merge(meta, how="left", left_on="source2", right_on="node_id")
        .rename(columns={"type":"type2", "year":"year2"})
        .drop(columns=["node_id"])
    )
    return final

@click.command(help='Compute co-citation strength from an edgelist CSV file.')
@click.option(
    '--edgelist_path',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help='Path to the edgelist CSV file.'
)
@click.option(
    '--output_aux_path',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to the output.aux, if your edgelist has an output.aux file (optional).'
)
@click.option(
    '--csv_path',
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help='Path to save the co-citation CSV file.'
)
@click.option(
    '--in_degree_threshold',
    type=int,
    default=0,
    help='Minimum in-degree threshold for target nodes.'
)
@click.option(
    '--n_nodes',
    type=int,
    default=None,
    help='Number of top target nodes to consider.'
)
@click.option(
    '--threshold-co-strength',
    'threshold_co_strength',
    type=int, default=0,
    help='Minimum co-citation strength to consider. Co-strength must be GREATER than this value.'
)
def click_co_citation(edgelist_path: Path, output_aux_path: Path, csv_path: Path, in_degree_threshold: int, n_nodes: int, threshold_co_strength: int):
    top_ids, sources_by_target = preprocess_edgelist(edgelist_path, in_degree_threshold, n_nodes)
    total_pairs = len(top_ids) * (len(top_ids) - 1) // 2
    start = time.perf_counter()
    #avoid storing all results in RAM memory
    with open(csv_path, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                            quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['source1', 'source2', 'Co-strength'])
        with click.progressbar(length=total_pairs, label='Computing co-citation') as bar:
            for i, (a, b) in enumerate(combinations(top_ids, 2), start=1):
                elapsed = time.perf_counter() - start
                bar.label = (
                    f'Co-citation {i}/{total_pairs} | '
                    f'Elapsed: {elapsed:.1f}s'
                )

                s1 = sources_by_target.get(a, set())
                s2 = sources_by_target.get(b, set())
                strength = len(s1 & s2) if len(s1) < len(s2) else len(s2 & s1)
                bar.update(1)
                if strength > threshold_co_strength:
                    # Write to CSV
                    spamwriter.writerow([a, b, strength])

    if output_aux_path is not None:
        df = pd.read_csv(csv_path)
        df_out = pd.read_csv(output_aux_path)
        df = merge_aux_co_strength(df, df_out)
        df.to_csv(csv_path, index=False)

    print(f"Co-citation data saved to {csv_path}")

if __name__ == "__main__":
    click_co_citation()