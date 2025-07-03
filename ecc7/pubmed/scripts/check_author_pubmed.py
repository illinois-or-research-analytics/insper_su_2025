import os
import re
import time
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from Bio import Entrez

Entrez.email = "ecc7@illinois.edu" 

def extract_country(affiliation_text):
    """
    Extracts the last word or segment after last comma as country from affiliation.
    """
    if affiliation_text:
        last = affiliation_text.strip().split(",")[-1].strip()
        return re.sub(r'[^\w\s]', '', last)
    return "Unknown"

def get_affiliations_and_year(xml_data):
    """
    Parses XML data to extract affiliations and publication year.
    """
    affiliations = []
    year = None

    root = ET.fromstring(xml_data)

    # Affiliations
    for aff in root.findall(".//AffiliationInfo/Affiliation"):
        if aff.text:
            affiliations.append(aff.text)

    year = root.findtext(".//PubDate/Year")

    return affiliations, year

def main():
    output_dir = "../outputs"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "exosome_affiliations.csv")

    print("Fetching PMIDs for exosome-related articles from PubMed...")
    handle = Entrez.esearch(
        db="pubmed",
        term='("exosomes"[MeSH Terms]',
        retmax=100000,
        mindate="1983", maxdate="2025"
    )
    record = Entrez.read(handle)
    pmids = record["IdList"]
    print(f"Found {len(pmids)} articles.")

    # Step 2: Iterate over PMIDs and extract data
    rows = []

    for pmid in tqdm(pmids, desc="Processing articles"):
        success = False
        while not success:
            try:
                fetch = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
                xml_data = fetch.read()
                fetch.close()

                affiliations, year = get_affiliations_and_year(xml_data)
                countries = [extract_country(aff) for aff in affiliations]

                # For each unique country mentioned in the paper
                for country in set(countries):
                    rows.append({
                        "PMID": pmid,
                        "Year": year,
                        "Country": country
                    })
                success = True
            except Exception as e:
                print(f"Error on PMID {pmid}, retrying in 5s...")
                time.sleep(5)

    # Step 3: Save results
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"Saved data to {output_path}")

if __name__ == "__main__":
    main()
