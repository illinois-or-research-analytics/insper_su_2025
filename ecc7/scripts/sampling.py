import os
import pandas as pd
import numpy as np
import time
import seaborn as sns
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from tqdm import tqdm
from Bio import Entrez

Entrez.email = "ecc7@illinois.edu"

def is_journal_article(xml_data):
    """
    Checks whether a PubMed article (in XML format) is of type 'Journal Article'

    Parameters:
        xml_data (str): Raw XML string of the article

    Returns:
        bool: True if the article is a journal article, False otherwise
    """
    root = ET.fromstring(xml_data)
    for article in root.findall(".//PubmedArticle"):
        for article_type in article.findall(".//PublicationType"):
            if article_type.text == "Journal Article":
                return True
            
    return False

def get_references(xml_data):
    """
    Fetches the list of referenced PubMed IDs for a given article

    Parameters:
        xml_data (str): Raw XML string of the article

    Returns:
        count (int): Number of references found
    """
    count = 0
    root = ET.fromstring(xml_data)
    for ref in root.findall(".//Reference"):
        article_id = ref.find(".//ArticleId[@IdType='pubmed']")
        if article_id is not None:
            count += 1
    return count

def build_sample_from_range(id_2020, id_2025, sample_size=10, output_dir="../outputs", filename="pubmed_references.csv"):
    """
    Generates a sample of PubMed journal articles with at least 5 references,
    selected randomly between a given range of PMIDs

    Parameters:
        id_2020 (str): PMID from year 2020
        id_2025 (str): PMID from year 2025
        sample_size (int): Number of valid samples to collect
        output_dir (str): Directory to store the CSV output
        filename (str): Name of the output CSV file

    Returns:
        pd.DataFrame: DataFrame containing sampled PMIDs and their out-degree (reference count)
    """
    valid_pmids = []
    refs = []

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)

    # Sample articles until the desired number of valid articles is collected
    for _ in tqdm(range(sample_size)):
        pmid = str(np.random.choice(range(int(id_2020), int(id_2025) + 1)))

        while True:
            try:
                handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
                xml_data = handle.read()
                handle.close()

                reference_pmids = get_references(xml_data)

                # Accept only if the article has at least 5 references
                if reference_pmids >= 5 and is_journal_article(xml_data):
                    valid_pmids.append(pmid)
                    refs.append(reference_pmids)
                    break
                else:
                    # Try another random PMID
                    pmid = str(np.random.choice(range(int(id_2020), int(id_2025)+ 1)))
            except:
                # If Entrez fetch fails, wait and retry
                time.sleep(5)

    # Save results to CSV in format expected by the ABM
    df = pd.DataFrame({
        'node_id': valid_pmids,
        'out_degree': refs
    })
    df.to_csv(output_path, index=False)

    return df

# Example run: generate a sample of 10,000 articles
df = build_sample_from_range(
    id_2020='31901868',
    id_2025='40486797',
    sample_size=10000,
    output_dir="../outputs",
    filename="pubmed_references_10k.csv"
)