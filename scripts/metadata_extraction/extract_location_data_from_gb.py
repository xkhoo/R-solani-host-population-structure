from Bio import Entrez
import csv

# Read the accessions from a file
accessions_file = 'data/accessions.txt'
with open(accessions_file) as f:
    ids = f.read().split('\n')

# Fetch the entries from Entrez
Entrez.email = 'xxxx@gmail.com'  # Insert your email here
handle = Entrez.efetch('nuccore', id=ids, retmode='xml')
response = Entrez.read(handle)

# Parse the entries to get the country
def extract_countries(entry):
    sources = [feature for feature in entry['GBSeq_feature-table']
               if feature['GBFeature_key'] == 'source']

    for source in sources:
        for qual in source['GBFeature_quals']:
            if qual['GBQualifier_name'] in ['country', 'geo_loc_name']:
                yield qual['GBQualifier_value']
        
for entry in response:
    accession = entry['GBSeq_primary-accession']
    locations = list(extract_countries(entry))

    if locations:
        for location in locations:
            print(accession, location, sep=',')
    else:
        print(accession, 'NA', sep=',')

# Write output to CSV
with open('output/location.csv', 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(['accession', 'location'])

    for entry in response:
        accession = entry['GBSeq_primary-accession']
        for host in extract_countries(entry):
            writer.writerow([accession, host])
