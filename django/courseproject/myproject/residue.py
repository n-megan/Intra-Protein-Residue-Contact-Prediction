import requests
from Bio.PDB import MMCIFParser
import numpy as np
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

'''
This module contains the Residue class and the amino acid mapping dictionary.
After obtaining the protein information from the user input, residue information (structural data) is obtained from the mmCIF file from the AlphaFold API.
'''

# Amino acid mapping for 3-letter code to full name
def amino_acid_mapping(residue):
    aa_full_names = {
        'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine', 'ASP': 'Aspartic Acid',
        'CYS': 'Cysteine', 'GLU': 'Glutamic Acid', 'GLN': 'Glutamine', 'GLY': 'Glycine',
        'HIS': 'Histidine', 'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
        'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline', 'SER': 'Serine',
        'THR': 'Threonine', 'TRP': 'Tryptophan', 'TYR': 'Tyrosine', 'VAL': 'Valine'
    }
    return aa_full_names.get(residue, None)  # Return the full name of the residue


class ResidueInfo:
    def __init__(self, input_protein):
        self.input_protein = input_protein

    # Get response from AlphaFold API 
    def alphafold_api(self):
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{self.input_protein}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            None

    # Get mmCIF file from AlphaFold API
    def get_mmcif_file(self):
        data = self.alphafold_api() 
        if data is None:
            return None
        mmcif_url = data[0]['cifUrl']
        mmcif_response = requests.get(mmcif_url)
        if mmcif_response.status_code == 200:
            mmcif_data = mmcif_response.text
            parser = MMCIFParser()
            mmcif_io = io.StringIO(mmcif_data)
            structure = parser.get_structure("mmcif_structure", mmcif_io)
            return structure  
        else:
            return None

    # Get AlphaFold source information
    def get_alphafold_source_info(self):
        data = self.alphafold_api()
        if data is None:
            return None       
        entry_id = data[0]['entryId']
        model_created_date = data[0]['modelCreatedDate']
        latest_version = data[0]['latestVersion']
        is_reviewed = data[0]['isReviewed']
        return {
            "Entry ID": entry_id, 
            "Created Date": model_created_date, # created_date
            "Latest Version": latest_version, # latest version
            "Is Reviewed": is_reviewed # is_reviewed
        }

    # Get all residues within a protein
    def get_residue_data(self, structure):
        total_residues = []
        if structure is not None:
            for model in structure:
                for chain in model:
                    for residue in chain.get_residues():
                        total_residues.append(residue)  
        return total_residues

    # This function returns detailed residue information, including chain ID and residue position
    def get_residue_details(self, structure):
        residue_info = []
        if structure is not None:
            for model in structure:
                for chain in model:
                    for residue in chain.get_residues():
                        residue_name = residue.get_resname()  # e.g., 'ALA', 'CYS'
                        residue_full_name = amino_acid_mapping(residue_name) # Full name of residue
                        residue_pos = residue.get_id()[1]  # Residue position (index) in chain
                        chain_id = chain.get_id()  # Chain ID (e.g., 'A', 'B', etc.)
                        
                        # Store the generic residue info for database or results
                        residue_info.append({
                            'residue_name': residue_name,
                            'residue_full_name': residue_full_name,
                            'residue_pos': residue_pos,
                            'chain_id': chain_id,
                            'total_residues': len(chain) 
                        })
        return residue_info

    # Get all atomic coordinates of specified atom type for each residue
    def get_coords(self, residue, atom_type):
        coords = []
        for atom in residue:
            if atom.get_name() == atom_type:
                coords.append(atom.get_coord())
        return coords if coords else None

    # Get unique residue pairs within the protein
    def get_residue_pairs(self, total_residues):
        residue_pairs = []
        for i, res1 in enumerate(total_residues):
            for j, res2 in enumerate(total_residues):
                if i < j: 
                    residue_pairs.append((res1, res2))
        return residue_pairs

    # Calculate the Euclidean distance between all CA or CB atoms for each residue pair
    def calc_euc_distance(self, res1, res2, atom_type="CA"):
        coords1 = self.get_coords(res1, atom_type)
        coords2 = self.get_coords(res2, atom_type)
        if not coords1 or not coords2:
            return []
        distances = [np.linalg.norm(coord1 - coord2) for coord1 in coords1 for coord2 in coords2]
        return distances

    # Determine residue contacts based on distances between all CA or CB atoms
    def calc_residue_contacts(self, residue_pairs, atom_type=None, threshold=5.0):
        residue_contacts = []
        for res1, res2 in residue_pairs:
            distances = self.calc_euc_distance(res1, res2, atom_type)
            for distance in distances:
                if distance <= threshold:
                    residue_contacts.append(((res1.id, res2.id), distance, res1, res2))
        return residue_contacts

    # Collect the CA and CB residue contacts with distance cutoff
    def get_residue_contacts(self, structure, threshold):
        total_residues = self.get_residue_data(structure)
        residue_pairs = self.get_residue_pairs(total_residues)
        ca_contacts = self.calc_residue_contacts(residue_pairs, atom_type="CA", threshold=threshold)
        cb_contacts = self.calc_residue_contacts(residue_pairs, atom_type="CB", threshold=threshold)
        self.contact_data = {'ca_contacts': ca_contacts, 'cb_contacts': cb_contacts}
        return ca_contacts, cb_contacts
    
    # Get CA residue contacts information
    def ca_residue_contacts(self, contacts, atom_type="CA"): 
        ca_contact_data = []
        if not contacts:
            print(f"No {atom_type} residue contacts found.")
            return ca_contact_data
        for contact in contacts:
            (res1_id, res2_id), distance, res1, res2 = contact
            coord1 = self.get_coords(res1, atom_type)
            coord2 = self.get_coords(res2, atom_type)
            if coord1 is not None and coord2 is not None:
                ca_contact_data.append({
                    "res1_full_name": res1.get_resname(),
                    "res2_full_name": res2.get_resname(),
                    "res1_pos": res1_id[1],
                    "res2_pos": res2_id[1],
                    "distance": round(distance, 2),
                    "res1_coord1": coord1,
                    "res2_coord2": coord2
                })
            else:
                print(f"Coordinates for {atom_type} atoms in residue pair {res1.get_resname()} {res1_id[1]} - {res2.get_resname()} {res2_id[1]} could not be retrieved.")
        return ca_contact_data

    # Get CB residue contacts information
    def cb_residue_contacts(self, contacts, atom_type="CB"): 
        cb_contact_data = []
        if not contacts:
            print(f"No {atom_type} residue contacts found.")
            return cb_contact_data
        for contact in contacts:
            (res1_id, res2_id), distance, res1, res2 = contact
            coord1 = self.get_coords(res1, atom_type)
            coord2 = self.get_coords(res2, atom_type)
            if coord1 is not None and coord2 is not None:
                cb_contact_data.append({
                    "res1_full_name": res1.get_resname(),
                    "res2_full_name": res2.get_resname(),
                    "res1_pos": res1_id[1],
                    "res2_pos": res2_id[1],
                    "distance": round(distance, 2),
                    "res1_coord1": coord1,
                    "res2_coord2": coord2
                })
            else:
                print(f"Coordinates for {atom_type} atoms in residue pair {res1.get_resname()} {res1_id[1]} - {res2.get_resname()} {res2_id[1]} could not be retrieved.")
        return cb_contact_data

    # Create contact map from residue pairs
    def create_contact_map(self, contact_data, num_residues):
        contact_map = np.zeros((num_residues, num_residues)) # Initialize square matrix
        for contact in contact_data: # Unpack the list of tuples
            (residue1, residue2), distance, res1_obj, res2_obj = contact
            res1_pos = residue1[1] - 1  # Adjust for zero-indexing
            res2_pos = residue2[1] - 1  # Adjust for zero-indexing
            contact_map[res1_pos][res2_pos] = distance  # Symmetric contact map
            contact_map[res2_pos][res1_pos] = distance  # Symmetric contact map
        return contact_map

    # Plot the contact map to display on GUI
    def plot_contact_map(self, contact_map, title, save_path):
        plt.figure(figsize=(6, 4))
        plt.imshow(contact_map, cmap='bone', interpolation='nearest')
        plt.title(title)
        plt.colorbar()
        plt.savefig(save_path, bbox_inches='tight')  # Save the plot with tight bounding box to avoid clipping
        plt.clf()
