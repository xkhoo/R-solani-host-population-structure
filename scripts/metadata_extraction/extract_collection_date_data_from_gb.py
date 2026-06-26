from Bio import Entrez
import csv

# Read the accessions from a file
accessions_file = 'data/accessions.txt'
with open(accessions_file) as f:
    ids = f.read().split('\n')

# Fetch the entries from Entrez
Entrez.email = 'xxxxx@gmail.com'  # Insert your email here
handle = Entrez.efetch('nuccore', id=ids, retmode='xml')
response = Entrez.read(handle)

# Parse the entries to get the collection_date
def extract_countries(entry):
    sources = [feature for feature in entry['GBSeq_feature-table']
               if feature['GBFeature_key'] == 'source']

    for source in sources:
        qualifiers = [qual for qual in source['GBFeature_quals']
                      if qual['GBQualifier_name'] == 'collection_date']
        
        for qualifier in qualifiers:
            yield qualifier['GBQualifier_value']

for entry in response:
    accession = entry['GBSeq_primary-accession']
    for collection_date in extract_countries(entry):
        print(accession, collection_date, sep=',')

# Write output to CSV
with open('output/collection_date.csv', 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(['accession', 'collection_date'])

    for entry in response:
        accession = entry['GBSeq_primary-accession']
        for host in extract_countries(entry):
            writer.writerow([accession, host])
