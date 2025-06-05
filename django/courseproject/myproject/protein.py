import requests

'''
This module describes the logic necessary to retrieve protein information from the UniProt API based on the user-entered UniProt ID.
'''

# Get repsonse from UniProt API
def uniprot_api(input_protein):
    url = f"https://rest.uniprot.org/uniprotkb/{input_protein}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        None

# Get UniProt source information
def get_uniprot_source_info(input_protein):
    data = uniprot_api(input_protein)
    if data is None:
        return 
    else:
        entry_audit = data.get("entryAudit", {})
        entry_type = data.get("entryType")
        return {
            "Entry Type": entry_type, # entry_id
            "First Public Date": entry_audit.get("firstPublicDate", "Unavailable"), # created_date
            "Last Annotation Update": entry_audit.get("lastAnnotationUpdateDate", "Unavailable"), # last_updated
            "Last Sequence Update": entry_audit.get("lastSequenceUpdateDate", "Unavailable"),
            "Entry Version": entry_audit.get("entryVersion", "Unavailable"), # entry_version
            "Sequence Version": entry_audit.get("sequenceVersion", "Unavailable"),
        }
    
# Get general protein information
def get_protein_data(input_protein):
    data = uniprot_api(input_protein)  
    if data is None:
        return None
    else:
        pdb_ids = [
            ref['id'] for ref in data.get('uniProtKBCrossReferences', [])
            if ref['database'] == 'PDB'
        ]
        gene = data.get('genes', [{}])[0].get('geneName', {}).get('value', 'Unavailable')
        organism = data.get('organism', {}).get('scientificName', 'Unavailable')
        protein_name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unavailable')
        sequence = data.get('sequence', {}).get('value', 'Unavailable')
        seq_length = data.get('sequence', {}).get('length', 'Unavailable')
        return {
            "PDB IDs": ', '.join(pdb_ids) if pdb_ids else "Unavailable",
            "Gene": gene,
            "Organism": organism,
            "Protein Name": protein_name,
            "Sequence": sequence,
            "Length": seq_length
        }