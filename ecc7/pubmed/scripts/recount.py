import os
import time
import click
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from Bio import Entrez

# CONFIG
Entrez.email = "estherccr@al.insper.edu.br"
Entrez.api_key = "28129c1de2ddd5ea6d4ae116eaf755977609"  

def get_references(xml_data, id_types):
    """
    Counts all references in a PubMed article XML and returns the count and ID types.
    """
    count = 0
    root = ET.fromstring(xml_data)

    for ref in root.findall(".//Reference"):
        for article_id in ref.findall(".//ArticleId"):
            id_type = article_id.attrib.get("IdType", "N/A")
            if article_id is not None:
                count += 1
            id_types[id_type] = id_types.get(id_type, 0) + 1

    return count, id_types

def fetch_or_load_xml(pmid, xml_cache_dir):
    """
    Reads an XML file from cache or fetches it from PubMed if not cached.
    """
    os.makedirs(xml_cache_dir, exist_ok=True)
    xml_path = os.path.join(xml_cache_dir, f"{pmid}.xml")

    if os.path.exists(xml_path):
        with open(xml_path, 'r') as f:
            return f.read()

    try:
        handle = Entrez.efetch(db="pubmed", id=str(pmid), retmode="xml")
        xml_data = handle.read().decode('utf-8')
        handle.close()

        with open(xml_path, 'w') as f:
            f.write(xml_data)

        time.sleep(0.3) 
        return xml_data

    except Exception as e:
        print(f"Error fetching PMID {pmid}: {e}")
        return None

@click.command(help='Recounts references using cached XMLs or fetching them, and outputs a clean new CSV.')
@click.option('--input_csv', type=str, required=True, help='Path to the input CSV with "#node_id"')
@click.option('--output_dir', type=str, default='../outputs', help='Directory to save the new CSV')
@click.option('--filename', type=str, default='recounted_clean.csv', help='Name of the output CSV')
@click.option('--xml_cache_dir', type=str, default='../outputs/xml_cache', help='Directory to store XML files')
def recount_and_save(input_csv, output_dir, filename, xml_cache_dir):
    df = pd.read_csv(input_csv)
    updated_pmids = []
    updated_refs = []
    id_types = {}

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    id_type_file = os.path.join(output_dir, f'id_types_{filename}')

    for pmid in tqdm(df['#node_id']):
        xml_data = fetch_or_load_xml(pmid, xml_cache_dir)
        if xml_data is None:
            continue

        try:
            ref_count, id_types = get_references(xml_data, id_types)
            updated_pmids.append(pmid)
            updated_refs.append(ref_count)
        except Exception as e:
            print(f"Error parsing PMID {pmid}: {e}")

    new_df = pd.DataFrame({
        '#node_id': updated_pmids,
        'out_degree': updated_refs
    })
    new_df.to_csv(output_path, index=False)
    print(f"CSV saved to: {output_path}")


    df_id_types = pd.DataFrame.from_dict(id_types, orient='index', columns=['count'])
    df_id_types.index.name = 'type'
    df_id_types.reset_index(inplace=True)
    df_id_types.to_csv(id_type_file, index=False)
    print(f"ID types saved to: {id_type_file}")

if __name__ == "__main__":
    recount_and_save()
