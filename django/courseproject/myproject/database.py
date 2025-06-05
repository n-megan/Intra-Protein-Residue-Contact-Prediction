# views.py (or elsewhere)

from .models import ProteinIdentifier, Protein, Residue, ResidueInteraction
from protein import get_protein_data  # Assuming this function returns protein data

def insert_protein_data(protein_info):
    # Step 1: Insert into Protein_Identifier
    protein_identifier = ProteinIdentifier.objects.create(
        uniprot_id=protein_info['uniprot_id'],
        pdb_id=protein_info['pdb_id']
    )

    # Step 2: Insert into Protein, referencing Protein_Identifier's id_number
    protein = Protein.objects.create(
        protein_name=protein_info['Protein Name'],
        gene=protein_info['Gene'],
        organism=protein_info['Organism'],
        sequence=protein_info['Sequence'],
        protein_identifier=protein_identifier
    )

    # Step 3: Insert Residue data, referencing Protein's protein_id
    for residue in protein_info['residues']:
        Residue.objects.create(
            total_residues=protein_info['total_residues'],
            residue_name=residue['residue_name'],
            residue_full_name=residue['residue_full_name'],
            residue_pos=residue['residue_pos'],
            chain_id=residue['chain_id'],
            protein=protein
        )

    # Step 4: Insert Residue Interactions, referencing Protein's protein_id
    for interaction in protein_info['interactions']:
        ResidueInteraction.objects.create(
            contact_type=interaction['contact_type'],
            res1_name=interaction['res1_name'],
            res1_pos=interaction['res1_pos'],
            res1_coords=', '.join(map(str, interaction['res1_coords'])),
            res2_name=interaction['res2_name'],
            res2_pos=interaction['res2_pos'],
            res2_coords=', '.join(map(str, interaction['res2_coords'])),
            distance=interaction['distance'],
            protein=protein
        )

    return protein.id  # Return the protein ID
