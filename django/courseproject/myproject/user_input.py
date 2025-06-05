import requests
from .protein import get_protein_data


def get_user_input(protein_id):
    protein_id = protein_id.strip()
    # If user enters a UniProt ID
    if len(protein_id) == 6 and protein_id.isalnum() and protein_id[0].isalpha(): 
        protein_id = protein_id.upper()
        protein_info = get_protein_data(protein_id)
        if protein_info:
            return "success", protein_id, protein_info  # Valid UniProt ID with data
        else:
            return "no_data", protein_id, "No data available for the given UniProt ID"  # Valid UniProt ID, but no data
    # If user enters a PDB ID
    elif len(protein_id) == 4 and protein_id.isalnum(): 
        protein_id = protein_id.lower()
        uniprot_ids = pdb_to_uniprot_map(protein_id)     
        if uniprot_ids: # Multiple UniProt IDs found for the PDB ID
            if len(uniprot_ids) > 1:
                return "multiple_ids", uniprot_ids, None  
            else:
                protein_info = get_protein_data(uniprot_ids[0])
                if protein_info:
                    return "success", uniprot_ids[0], protein_info  # Valid PDB ID, with single UniProt ID
                else:
                    return "no_data", uniprot_ids[0], "No data available for the UniProt ID mapped from PDB ID"
        else: # No UniProt ID associated with the PDB ID
            return "no_uniprot", None, f"No UniProt ID found for PDB ID {protein_id}." 
    else: # Invalid input
        return "invalid_id", None, "Invalid input. Please enter a valid UniProt ID or PDB ID"



# Function to map PDB ID to UniProt ID if user enters a PDB ID
def pdb_to_uniprot_map(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"    
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        if pdb_id in data:
            uniprot_dict = data[pdb_id].get('UniProt', {})
            uniprot_ids = list(uniprot_dict.keys())        
            return uniprot_ids
        else:
            return None
    except requests.exceptions.RequestException:
        return None