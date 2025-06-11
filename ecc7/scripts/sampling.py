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
    root = ET.fromstring(xml_data)
    for article in root.findall(".//PubmedArticle"):
        for article_type in article.findall(".//PublicationType"):
            if article_type.text == "Journal Article":
                return True
    return False

def get_references(xml_data):
    reference_pmids = []
    root = ET.fromstring(xml_data)
    for ref in root.findall(".//Reference"):
        article_id = ref.find(".//ArticleId[@IdType='pubmed']")
        if article_id is not None:
            reference_pmids.append(article_id.text)

    return reference_pmids

def build_sample_from_range(id_2020, id_2025, sample_size=10, output_dir="outputs", filename="pubmed_references.csv"):
    valid_pmids = []
    refs = []

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)

    for _ in tqdm(range(sample_size)):
        pmid = str(np.random.choice(range(int(id_2020), int(id_2025))))

        while True:
            try:
                handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
                xml_data = handle.read()
                handle.close()
                reference_pmids = get_references(xml_data)
                if len(reference_pmids) >= 5 and is_journal_article(xml_data):
                    valid_pmids.append(pmid)
                    refs.append(len(reference_pmids))
                    break
                else:
                    pmid = str(np.random.choice(range(int(id_2020), int(id_2025))))
            except:
                time.sleep(10)

    df = pd.DataFrame({
        'node_id': valid_pmids,
        'out_degree': refs
    })
    df.to_csv(output_path, index=False)
    return df

# Example usage:
df = build_sample_from_range(id_2020='31901868', id_2025='40486797', sample_size=10000, output_dir="outputs", filename="pubmed_references_10k.csv")

