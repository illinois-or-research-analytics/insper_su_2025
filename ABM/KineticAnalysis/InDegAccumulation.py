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
@click.option('--output-dir', required= False, default= "", type=click.Path(), help='Name of the output directory')
@click.option('--exec-label', required= False, default= "", type=str, help='Name of the execution folder')
def plot_in_deg(input_dir, output_dir, exec_label):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, exec_label), exist_ok=True)

    output_dir = os.path.join(output_dir, exec_label)

    node_df = pd.read_csv(os.path.join(input_dir, "output.aux"))
    node_df = node_df[node_df["type"] == "agent"]
        
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

        
    node_df.to_csv(os.path.join(output_dir, f"in_degree_accumulation_{exec_label}.csv"))

    med_d = node_df[f"in_degree_until_{max(years)+1}"].median()
    d_75 = node_df[f"in_degree_until_{max(years)+1}"].quantile(0.75)
    d_25 = node_df[f"in_degree_until_{max(years)+1}"].quantile(0.25)

    # print(f"Median in-degree - {med_d}")
    # print(f"75th percentile - {d_75}")
    # print(f"25th percentile - {d_25}")

    med_node_id = node_df[node_df[f"in_degree_until_{max(years)+1}"] == med_d].iloc[0]["node_id"]
    d75_node_id = node_df[node_df[f"in_degree_until_{max(years)+1}"] == d_75].iloc[0]["node_id"]
    d25_node_id = node_df[node_df[f"in_degree_until_{max(years)+1}"] == d_25].iloc[0]["node_id"]

    # print(f"Node with median in-degree - {med_node_id}")
    # print(f"Node at 75th percentile - {d75_node_id}")
    # print(f"Node at 25th percentile - {d25_node_id}")


    med_df = node_df[node_df["node_id"] == med_node_id]
    d75_df = node_df[node_df["node_id"] == d75_node_id]
    d25_df = node_df[node_df["node_id"] == d25_node_id]

    all_dfs = {
        "percentile_75": d75_df,
        "median": med_df,
        "percentile_25": d25_df
    }
    
    if len(node_df[node_df["planted_nodes_line_number"]>0]) == 0:

        max_d = node_df[f"in_degree_until_{max(years)+1}"].max()
        # print(f"Max in-degree - {max_d}")
        max_node_id = node_df[node_df[f"in_degree_until_{max(years)+1}"] == max_d].iloc[0]["node_id"]
        # print(f"Node with max in-degree - {max_node_id}")
        max_df = node_df[node_df["node_id"] == max_node_id]
        all_dfs["max"] = max_df

    else:

        line_num = [n for n in node_df["planted_nodes_line_number"].dropna().unique().tolist() if n >= 0]
        for i in line_num:    

            node = node_df[node_df["planted_nodes_line_number"] == i]["node_id"].to_list()[0]
            
            print(node)
            
            local_node = node_df[node_df["node_id"] == node]
            
            print(local_node['fit_peak_value'].tolist()[0])
            
            all_dfs[f"Fitness {local_node['fit_peak_value'].tolist()[0]}"] = local_node
        
        


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
    plt.savefig(os.path.join(output_dir,f"in_degree_accumulation_{exec_label}.png"))
    plt.show()
    
    print("Done")

if __name__ == "__main__":
    plot_in_deg()