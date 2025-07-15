from typing import Dict, Sequence, Set, Tuple
import pandas as pd
from pathlib import Path
from itertools import combinations
import click
import csv
import time

def preprocess_edgelist(
    edgelist_path: Path,
    # output_aux_path: Path = None,
    out_degree_threshold: int = 0,
    n_nodes: int = None
) -> Tuple[Sequence[int], Dict[int, Set[int]]]:
    """
    Returns:
      top_ids:     list of the top-n citing IDs with out_degree > threshold
      targets_by_source: mapping source ID → set of target IDs
    """
    df = pd.read_csv(edgelist_path)

    first_col = df.columns[0]
    second_col = df.columns[1]

    df = df.rename(columns={first_col: '#source', second_col: 'target'})

    # if output_aux_path is not None:
    #     df_aux = pd.read_csv(output_aux_path)
    #     df_aux = df_aux[df_aux["type"] == "seed"]
    #     ids_seeds = df_aux["node_id"].value_counts().index.tolist()
    #     df = df[
    #         df['#source'].isin(ids_seeds) &
    #         df['target'].isin(ids_seeds)
    #     ]

    vc = df['#source'].value_counts()
    filtered = vc[vc > out_degree_threshold]

    if n_nodes is None:
        top_ids = filtered.index.tolist()
    else:
        top_ids = filtered.head(n_nodes).index.tolist()

    targets_by_source = (
        df[df['#source'].isin(top_ids)]
          .groupby('#source')['target']
          .agg(set)
          .to_dict()
    )

    return top_ids, targets_by_source

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

@click.command(help='Compute bibliographic coupling from an edgelist CSV file.')
@click.option(
    '--edgelist_path',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help='Path to the edgelist CSV file.'
)
@click.option(
    '--output_aux_path',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to the auxiliary CSV file with seed nodes (optional).'
)
@click.option(
    '--csv_path',
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help='Path to save the bibliographic coupling CSV file.'
)
@click.option(
    '--out_degree_threshold',
    type=int,
    default=0,
    help='Minimum out-degree threshold for source nodes.'
)
@click.option(
    '--n_nodes',
    type=int,
    default=None,
    help='Number of top source nodes to consider.'
)
@click.option(
    '--threshold-biblio-strength',
    'threshold_bi_strength',
    type=int, default=0,
    help='Minimum biblio-strength to consider. Biblio-strength must be GREATER than this value.'
)
def click_bibliographic_coupling(
    edgelist_path: Path,
    output_aux_path: Path,
    csv_path: Path,
    out_degree_threshold: int,
    n_nodes: int,
    threshold_bi_strength: int = 0
):
    top_ids, targets_by_source = preprocess_edgelist(
        edgelist_path,
        out_degree_threshold,
        n_nodes
    )
    total_pairs = len(top_ids) * (len(top_ids) - 1) // 2
    start = time.perf_counter()
    with open(csv_path, 'w', newline='') as csvfile:
        spawriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        spawriter.writerow(['source1', 'source2', 'Biblio-strength'])
        with click.progressbar(length=total_pairs, label='Computing bibliographic coupling') as bar:
            for i, (a, b) in enumerate(combinations(top_ids, 2), start=1):
                elapsed = time.perf_counter() - start
                bar.label = (
                    f'Biblio {i}/{total_pairs} |'
                    f'Elapsed: {elapsed:.2f}s'
                )

                c1 = targets_by_source.get(a, set())
                c2 = targets_by_source.get(b, set())
                strength = len(c1 & c2) if len(c1) < len(c2) else len(c2 & c1)
                bar.update(1)
                if strength > threshold_bi_strength:
                    spawriter.writerow([a, b, strength])

    if output_aux_path is not None:
        df = pd.read_csv(csv_path)
        df_out = pd.read_csv(output_aux_path)
        df = merge_aux_co_strength(df, df_out)
        df.to_csv(csv_path, index=False)
    print(f'Bibliographic coupling data saved to {csv_path}')

if __name__ == '__main__':
    click_bibliographic_coupling()