import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import click
import math


def get_in_degree_cummulative(node_df, edge_df, year):
    to_node_counts = edge_df['target'].value_counts()

    node_df[str(year)] = node_df['node_id'].map(to_node_counts)


def extract_years_from_filenames(directory):
    years = set()
    year_pattern = re.compile(r'output\.edgelist_(\d{4})')

    for filename in os.listdir(directory):
        match = year_pattern.match(filename)
        if match:
            years.add(int(match.group(1)))
            
    years = list(years)

    return sorted(years)


def plot_quartile_fill(df, label, color, years):

    min_vals = df[years].min(axis=0)
    max_vals = df[years].max(axis=0)

    plt.fill_between(years, min_vals, max_vals, alpha=0.2, color=color, label=label)



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
    
    # node_df.to_csv(os.path.join(output_dir, f"in_degree_accumulation_{exec_label}.csv"))
    
    years.append(max(years)+1)
    years = sorted(years)

    node_df = node_df.dropna(subset=[str(max(years))])

    node_df['quartile'] = pd.qcut(
        node_df[str(max(years))],
        q=4,
        labels=[1, 2, 3, 4],
        duplicates='drop'
    ).astype(int)

    node_df['quartile'] = node_df['quartile'].astype(int)

    node_df_q1 = node_df[node_df['quartile'] == 1]
    node_df_q2 = node_df[node_df['quartile'] == 2]
    node_df_q3 = node_df[node_df['quartile'] == 3]
    node_df_q4 = node_df[node_df['quartile'] == 4]


    plt.figure(figsize=(12, 6))
    
    year_cols = [str(year) for year in years]

    plot_quartile_fill(node_df_q1, label='Q1 (lowest)', color='red', years=year_cols)
    plot_quartile_fill(node_df_q2, label='Q2', color='blue', years=year_cols)
    plot_quartile_fill(node_df_q3, label='Q3', color='green', years=year_cols)
    plot_quartile_fill(node_df_q4, label='Q4 (highest)', color='orange', years=year_cols)
    
    med_vals = node_df[year_cols].median(axis=0)
    
    plt.plot(year_cols, med_vals, label="Median In-Degree Value (Overall)", color="black", linewidth=2)
    
    n = math.ceil(math.log10(node_df[str(max(years))].max()))

    plt.title(f"In-degree over time for Nodes Quartiles - {exec_label}")
    plt.xlabel("Year")
    plt.ylabel("In-degree")
    plt.yscale("log")
    plt.xticks(rotation=90)
    plt.ylim(1, 10000) 
    plt.yticks([10**i for i in range(1, n + 1)])
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir,f"in_degree_accumulation_by_quartiles_{exec_label}.png"))
    plt.show()
    
    print("Done")

if __name__ == "__main__":
    plot_in_deg()