import pandas as pd

def get_metrics(outaux, clustering_metrics, leiden_clustering) -> dict:
    df_outaux = pd.read_csv(outaux)
    df_outaux = df_outaux.set_index('node_id') 

    df_clustering_metrics = pd.read_csv(clustering_metrics)
    df_leiden_clustering = pd.read_csv(leiden_clustering)

    metrics = {}

    cluster_ids = df_clustering_metrics['cluster_id'].unique()
    for cluster_id in cluster_ids:
        recency_w = 0
        pa_w= 0
        fitness_w = 0
        fitness = 0
        seeds = 0
        agents = 0

        node_ids = df_leiden_clustering[df_leiden_clustering['cluster_id'] == cluster_id]['node_id'].tolist()

        for node_id in node_ids:
            if node_id in df_outaux.index:
                row = df_outaux.loc[node_id]
                if row['type'] == 'seed':
                    seeds += 1
                else:
                    agents += 1
                    recency_w += row['rec_weight']
                    pa_w += row['pa_weight']
                    fitness_w += row['fit_weight']
                fitness += row['fit_peak_value']
        
        if agents == 0: #cluster with only seeds
            avg_recency_w = None
            avg_pa_w = None
            avg_fitness_w = None
        else:
            size = len(node_ids)
            avg_recency_w = f"{float(recency_w / agents):.4f}"
            avg_pa_w = f"{float(pa_w / agents):.4f}"
            avg_fitness_w = f"{float(fitness_w / agents):.4f}"

        avg_fitness = f"{float(fitness / size):.4f}"
            

        metrics[cluster_id] = {
            'avg_recency_w': avg_recency_w,
            'avg_pa_w': avg_pa_w,
            'avg_fitness_w': avg_fitness_w,
            'avg_fitness': avg_fitness,
            'seeds': int(seeds),
            'agents': int(agents)
        }



    metrics_df = pd.DataFrame.from_dict(metrics, orient='index')
    metrics_df.index.name = 'cluster_id'  
    df_clustering_metrics = df_clustering_metrics.merge(metrics_df, on='cluster_id', how='left')
    return df_clustering_metrics



if __name__ == "__main__":
    outaux = 'test_data/7_6_fmf_standard_no_ss_ra_new_pubmed/output/output.aux'
    clustering_metrics = 'test_data/7_6_fmf_standard_no_ss_ra_new_pubmed/output/clustering_metrics_01.csv'
    leiden_clustering = 'test_data/7_6_fmf_standard_no_ss_ra_new_pubmed/output/leidenClustering_0.1.csv'

    metrics_df = get_metrics(outaux, clustering_metrics, leiden_clustering)
    metrics_df.to_csv('test_data/7_6_fmf_standard_no_ss_ra_new_pubmed/output/clustering_metrics_with_weights.csv', index=False)

