import matplotlib.pyplot as plt
import pandas as pd
import click
from pathlib import Path
import numpy as np

def save_table(df: pd.DataFrame,
               filename: str = "table.png",
               scale: tuple[float, float] = (1, 1.5),
               font_size: int = 10,
               dpi: int = 300):

    width = max(4, len(df.columns) * 1.2)
    height = max(2, df.shape[0] * 0.5 + 1)
    fig, ax = plt.subplots(figsize=(width, height))

    ax.axis('off')

    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(*scale)

    plt.savefig(filename, bbox_inches='tight', dpi=dpi)
    plt.close(fig)

def get_metrics_by_id(exp_path, cluster_id, gamma, ss_ids):
    exp_path = str(exp_path)
    gamma_no_point = str(gamma).replace(".", "")
    exp = "ss_"
    exp_name = exp_path.split("/")[-1]
    if "ra" in exp_name:
        exp += "ra_"
    else:
        exp += "sa_"

    if "er" in exp_name:
        exp += "er"
    else:
        exp += "sj"


    clustering_metrics_file = "clustering_metrics_" + str(gamma_no_point) + ".csv"
    leiden_clustering_file = "leidenClustering_" + str(gamma) + ".csv"
    output_aux_file = "output.aux"

    df_metrics = pd.read_csv(f"{exp_path}/output/{clustering_metrics_file}")
    df_leiden = pd.read_csv(f"{exp_path}/output/{leiden_clustering_file}")
    df_aux = pd.read_csv(f"{exp_path}/output/{output_aux_file}", low_memory=False)

    df_leiden = df_leiden[df_leiden["cluster_id"] == cluster_id]
    df_leiden = df_leiden[df_leiden["node_id"].isin(ss_ids)]
    ids_leiden = df_leiden["node_id"].to_list()
    df_aux = df_aux[df_aux["node_id"].isin(ids_leiden)]
    fit_list = df_aux["fit_peak_value"].to_list()
    string = ""
    for fit in fit_list:
        if fit == 10000:
            string += "10k "
        elif fit == 100000:
            string += "100k "
        elif fit == 1000000:
            string += "1M "

    num_ss = df_leiden.shape[0]

    df_metrics = df_metrics[df_metrics["cluster_id"] == cluster_id]
    df_metrics.drop(columns=["intra_edges", "boundary_edges", "normalized_density"], inplace=True)
    df_metrics["experiment"] = exp
    df_metrics["ss_count"] = num_ss
    df_metrics["ss_fit"] = string
    

    dict_metrics = df_metrics.to_dict(orient="records")
    return dict_metrics

def get_pheno_network(edgelist_path, aux_path):

    df_edge_list = pd.read_csv(edgelist_path)
    df_aux = pd.read_csv(aux_path, low_memory=False)
    weights = {}

    ss_ids = df_aux[df_aux["fit_peak_value"] > 1000]["node_id"].to_list()
    for id_ in ss_ids:

        citings_ids = df_edge_list[df_edge_list["target"] == id_]["#source"].tolist()

        filtered_aux = df_aux[df_aux["node_id"].isin(citings_ids)]
        filtered_aux = filtered_aux[filtered_aux["type"] == "agent"]

        average_pa_weight_net = filtered_aux.groupby("year")["pa_weight"].mean().reset_index()
        average_rec_weight_net = filtered_aux.groupby("year")["rec_weight"].mean().reset_index()
        average_fit_weight_net = filtered_aux.groupby("year")["fit_weight"].mean().reset_index()

        ss_fit = df_aux[df_aux["node_id"] == id_]["fit_peak_value"].values[0]
        weights[str(ss_fit)] = {
            "pa_weight": average_pa_weight_net,
            "rec_weight": average_rec_weight_net,
            "fit_weight": average_fit_weight_net
        }
    return weights

def get_pheno_cluster(leidenClustering_path, aux_path):

    df_leiden = pd.read_csv(leidenClustering_path)
    df_aux = pd.read_csv(aux_path, low_memory=False)
    weights = {}

    ss_ids = df_aux[df_aux["fit_peak_value"] > 1000]["node_id"].to_list()
    for id in ss_ids:
        cluster_id = df_leiden[df_leiden["node_id"] == id]["cluster_id"].values[0]
        all_cluster_ids = df_leiden[df_leiden['cluster_id'] == cluster_id]["node_id"].to_list()

        filtered_aux = df_aux[df_aux['node_id'].isin(all_cluster_ids)]
        filtered_aux = filtered_aux[filtered_aux["type"] == "agent"]

        average_pa_weight = filtered_aux.groupby("year")["pa_weight"].mean().reset_index()
        average_rec_weight = filtered_aux.groupby("year")["rec_weight"].mean().reset_index()
        average_fit_weight = filtered_aux.groupby("year")["fit_weight"].mean().reset_index()
        
        if cluster_id not in weights:
            weights[str(cluster_id)] = {
                "pa_weight": average_pa_weight,
                "rec_weight": average_rec_weight,
                "fit_weight": average_fit_weight
            }

    return weights, ss_ids


@click.command()
@click.option(
    "--exp-dirs", "-e",
    multiple=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="One or more experiment directory paths"
)
@click.option(
    "--gamma", "-g",
    type=float,
    required=True,
    help="Gamma value for the leidenClustering files"
)
def cli(exp_dirs, gamma):

    for exp_path in exp_dirs:

        leiden_clustering_file = exp_path / "output" / f"leidenClustering_{str(gamma)}.csv"
        edgelist_path = exp_path / "output" / "output.edgelist"
        aux_path = exp_path / "output" / "output.aux"

        net_weights = get_pheno_network(edgelist_path, aux_path)
        cluster_weights, ss_ids = get_pheno_cluster(leiden_clustering_file, aux_path)

        n = len(net_weights)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=True)
        axes = np.atleast_1d(axes)

        for ax, (key, value) in zip(axes, net_weights.items()):
            pa = value["pa_weight"].sort_values("year")
            rec = value["rec_weight"].sort_values("year")
            fit = value["fit_weight"].sort_values("year")

            ax.plot(pa["year"], pa["pa_weight"],      label="Avg PA weight")
            ax.plot(rec["year"], rec["rec_weight"],   label="Avg Recency weight")
            ax.plot(fit["year"], fit["fit_weight"],   label="Avg Fitness weight")

            ax.set_title(f"ss_fit = {key}")
            ax.set_xlabel("Year")
            ax.grid(True)

        axes[0].set_ylabel("Average weight")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.92))
        fig.suptitle("AVG pheno weight for citing ss agents in the network", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig("all_weights_net.png", dpi=300)
        plt.close(fig)

        n = len(cluster_weights)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=True)
        axes = np.atleast_1d(axes)

        for ax, (key, value) in zip(axes, cluster_weights.items()):
            pa = value["pa_weight"].sort_values("year")
            rec = value["rec_weight"].sort_values("year")
            fit = value["fit_weight"].sort_values("year")

            ax.plot(pa["year"], pa["pa_weight"],      label="Avg PA weight")
            ax.plot(rec["year"], rec["rec_weight"],   label="Avg Recency weight")
            ax.plot(fit["year"], fit["fit_weight"],   label="Avg Fitness weight")

            ax.set_title(f"cluster_id = {key}")
            ax.set_xlabel("Year")
            ax.grid(True)

        axes[0].set_ylabel("Average weight")
        handles, labels = axes[0].get_legend_handles_labels()

        fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.92))
        fig.suptitle("AVG pheno weight for citing agents in ss clusters", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig("all_weights_cluster.png", dpi=300)
        plt.close(fig)

            
        used_cluster_ids = []
        df = pd.DataFrame()
        for key, value in cluster_weights.items():
            if key in used_cluster_ids:
                continue
            used_cluster_ids.append(key)
            dic = get_metrics_by_id(exp_path, int(key), gamma, ss_ids)
            df = pd.concat([df, pd.DataFrame(dic)], ignore_index=True)
        save_table(df)

if __name__ == "__main__":
    cli()