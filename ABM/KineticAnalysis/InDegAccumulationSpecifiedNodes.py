import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import click


def get_in_degree_cummulative(node_df, edge_df, year):
    to_node_counts = edge_df['target'].value_counts()

    node_df[f'in_degree_until_{year}'] = node_df['node_id'].map(to_node_counts)



def extract_years_from_filenames(directory):
    years = set()
    year_pattern = re.compile(r'output\.edgelist_(\d{4})')

    for filename in os.listdir(directory):
        match = year_pattern.match(filename)
        if match:
            years.add(int(match.group(1)))
            
    years = list(years)

    return sorted(years)

@click.command()
@click.option('--input-dir', required= True, default= "", type=click.Path(exists=True), help='Name of the input directory')
@click.option('--nodes-file', required= True, default= "", type=click.Path(exists=True), help='Name of the txt file with node IDs')
@click.option('--output-dir', required= False, default= "", type=click.Path(), help='Name of the output directory')
@click.option('--exec-label', required= False, default= "", type=str, help='Name of the execution folder')
def plot_in_deg(input_dir, nodes_file, output_dir, exec_label):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, exec_label), exist_ok=True)
    
    with open(nodes_file, 'r', encoding='utf-8') as f:
        nodes = [int(line.strip()) for line in f.readlines()]

    output_dir = os.path.join(output_dir, exec_label)
    
    node_df = pd.read_csv(os.path.join(input_dir, "output.aux"))
    node_df = node_df[node_df["node_id"].isin(nodes)]
    
    years = extract_years_from_filenames(input_dir)
    
    for year in years:
        
        edgelist_filename = f"output.edgelist_{year}"
        edge_df = pd.read_csv(os.path.join(input_dir, edgelist_filename))
        edge_df = edge_df[edge_df["target"].isin(node_df["node_id"])]
        get_in_degree_cummulative(node_df=node_df, edge_df=edge_df, year=year)
        
    edgelist_filename = f"output.edgelist"
    edge_df = pd.read_csv(os.path.join(input_dir, edgelist_filename))
    edge_df = edge_df[edge_df["target"].isin(node_df["node_id"])]
    get_in_degree_cummulative(node_df=node_df, edge_df=edge_df, year=max(years)+1)

        
    node_df.to_csv(os.path.join(output_dir, f"in_degree_accumulation_specified_nodes_{exec_label}.csv"))

    all_dfs = {}
    
    for nid in nodes:
        all_dfs[f"Node {nid}"] = node_df[node_df["node_id"] == nid]

    years.append(max(years)+1)

    plt.figure(figsize=(12, 6))

    for t in all_dfs:
        in_degrees = [all_dfs[t][f"in_degree_until_{year}"] for year in years]
        plt.plot(years, in_degrees, marker='o', label=t)

    plt.title(f"In-degree over time for Nodes - {exec_label}")
    plt.xlabel("Year")
    plt.ylabel("In-degree")
    plt.yscale("log")
    plt.ylim(1, 10000) 
    plt.yticks([1, 10, 100, 1000, 10000, 100000, 1000000])
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir,f"in_degree_accumulation_specified_{exec_label}.png"))
    plt.show()

if __name__ == "__main__":
    plot_in_deg()