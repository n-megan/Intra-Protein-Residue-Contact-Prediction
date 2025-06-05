from django.shortcuts import render, redirect
from django.contrib import messages
from django.http import JsonResponse
from django.conf import settings
from django.db import transaction
import os
from .user_input import *
from .protein import *
from .residue import *
from .models import ProteinIdentifier, PDB_IDS, Protein, Residue, Residue_Interactions, Source

'''
Django's view module designed to handle protein data retrieval, processing, and database storage.
'''

# AJAX request to fetch protein data
def get_protein_info_ajax(request, uniprot_id):
    protein_info = get_protein_data(uniprot_id)
    if protein_info:
        return JsonResponse(protein_info)
    else:
        messages.error(request, 'Protein structure unavailable for the provided input or invalid input.')
        return redirect('view_user_input')


# Helper function to process protein data
def process_protein_data(protein_id, protein_info):
    residue_info = ResidueInfo(protein_id)
    structure = residue_info.get_mmcif_file()
    if structure:
        uniprot_source_info = get_uniprot_source_info(protein_id)
        alphafold_source_info = residue_info.get_alphafold_source_info()
        residue_data = residue_info.get_residue_details(structure)
        ca_contacts, cb_contacts = residue_info.get_residue_contacts(structure, threshold=5.0)
        ca_contact_data = residue_info.ca_residue_contacts(ca_contacts)
        cb_contact_data = residue_info.cb_residue_contacts(cb_contacts)     

        # Create contact map
        num_residues = len(residue_info.get_residue_data(structure))
        static_path = os.path.join(settings.BASE_DIR, 'myproject', 'static', 'images')   
        if not os.path.exists(static_path):
            os.makedirs(static_path)    
        ca_contact_map_path = os.path.join(static_path, f'ca_contact_map_{protein_id}.png')
        cb_contact_map_path = os.path.join(static_path, f'cb_contact_map_{protein_id}.png')   
        residue_info.plot_contact_map(residue_info.create_contact_map(ca_contacts, num_residues), "Contact Map", ca_contact_map_path)
        residue_info.plot_contact_map(residue_info.create_contact_map(cb_contacts, num_residues), "Contact Map", cb_contact_map_path)
 
        chain_ids = [residue['chain_id'] for residue in residue_data]      
        chain_id = chain_ids[0] if chain_ids else None
        return {
            'protein_id': protein_id,
            'protein_info': protein_info,
            'chain_id': chain_id,
            'ca_contact_data': ca_contact_data,
            'cb_contact_data': cb_contact_data,
            'ca_contact_map_path': f'images/ca_contact_map_{protein_id}.png',
            'cb_contact_map_path': f'images/cb_contact_map_{protein_id}.png',
            'residue_data': residue_data, 
            'uniprot_source_info': uniprot_source_info,
            'alphafold_source_info': alphafold_source_info
        }
    else:  # If structure is not available, return protein info only
        error_message = "No structural data available to determine residue contacts"
        return {
            'protein_id': protein_id,
            'protein_info': protein_info,
            'error': error_message
        }


'''
This class processes protein data (UniProt ID, residue details, protein structure) and saves it to the database. 
It manages the creation of related instances (Protein_Identifier, PDB_IDS, Protein, Residue, Residue_Interactions and Source), handles exceptions and ensures data integrity using transactions.
'''
class DatabaseData:
    def __init__(self, result):
        self.result = result
        self.protein_info = result.get("protein_info", {})
        self.protein_id = result.get("protein_id")
        self.uniprot_instance = self.create_uniprot_instance()
        self.protein_instance = self.create_protein_instance()
        self.uniprot_source_info = result.get("uniprot_source_info", [])
        self.alphafold_source_info = result.get("alphafold_source_info", [])

    def save_data(self):
        try:
            with transaction.atomic():
                self.process_pdb_ids()
                self.create_residue_instances()
                self.create_residue_interactions(self.result, self.protein_instance)
                self.save_source_info(self.uniprot_source_info, self.alphafold_source_info, self.protein_instance)
                print("Data successfully saved.")
        except Exception as e:
            print(f"Error saving data: {e}")
            raise

    def create_uniprot_instance(self):
        try:
            uniprot_instance, created = ProteinIdentifier.objects.get_or_create(uniprot_id=self.protein_id)
            return uniprot_instance
        except Exception as e:
            print(f"Error creating UniProt instance: {e}")
            return None

    def create_protein_instance(self):
        try:
            protein_instance, created = Protein.objects.get_or_create(
                id_number=self.uniprot_instance,
                protein_name=self.protein_info.get("Protein Name", ""),
                gene=self.protein_info.get("Gene", ""),
                organism=self.protein_info.get("Organism", ""),
                total_residues=len(self.protein_info.get("residue_data", [])),
                sequence=self.protein_info.get("Sequence", "")
            )
            print(f"Protein {self.protein_info.get('Protein Name', '')} {'created' if created else 'already exists'}.")
            return protein_instance
        except Exception as e:
            print(f"Error creating Protein instance: {e}")
            return None

    def process_pdb_ids(self):
        pdb_id = self.protein_info.get("PDB IDs", [])
        if isinstance(pdb_id, str):
            pdb_id = pdb_id.split(",")        
        for pdb in pdb_id:
            pdb = pdb.strip()
            try:
                existing_pdb = PDB_IDS.objects.filter(pdb_id=pdb).first()
                if existing_pdb:
                    existing_pdbids = PDB_IDS.objects.filter(pdb_id=pdb, uniprot_id=self.uniprot_instance).first()
                    if not existing_pdbids:
                        PDB_IDS.objects.create(uniprot_id=self.uniprot_instance, pdb_id=pdb)
                else:
                    PDB_IDS.objects.create(uniprot_id=self.uniprot_instance, pdb_id=pdb)
            except Exception as e:
                print(f"Error processing PDB ID {pdb}: {e}")

    def create_residue_instances(self):
        residue_data = self.result.get("residue_data", [])
        for residue in residue_data:
            try:
                Residue.objects.create(
                    protein_id=self.protein_instance,
                    residue_name=residue.get("residue_name"),
                    residue_full_name=residue.get("residue_full_name"),
                    residue_pos=residue.get("residue_pos"),
                    chain_id=residue.get("chain_id"),
                )
            except Exception as e:
                print(f"Error creating Residue instance: {e}")

    def create_residue_interactions(self, result, protein_instance):
        def create_contact(contact_data, contact_type):
            for contact in contact_data:
                Residue_Interactions.objects.create(
                    contact_type=contact_type,
                    res1_name=contact["res1_full_name"],
                    res1_pos=contact["res1_pos"],
                    res1_coords=",".join(map(str, contact["res1_coord1"][:3])),
                    res2_name=contact["res2_full_name"],
                    res2_pos=contact["res2_pos"],
                    res2_coords=",".join(map(str, contact["res2_coord2"][:3])),
                    distance=contact["distance"],
                    protein_id=protein_instance
                )
        create_contact(result.get("ca_contact_data", []), "alpha-carbon")
        create_contact(result.get("cb_contact_data", []), "beta-carbon")

    def save_source_info(self, uniprot_source_info, alphafold_source_info, protein_instance):
        # print(f"UniProt Source Info: {uniprot_source_info}")
        # print(f"AlphaFold Source Info: {alphafold_source_info}")
        try:
            Source.objects.create(
                source_name="UniProt",
                source_version=uniprot_source_info.get("Entry Version", ""),
                created_date=uniprot_source_info.get("First Public Date", ""),
                protein_id=protein_instance,
            )
            Source.objects.create(
                source_name="AlphaFold",
                source_version=alphafold_source_info.get("Latest Version", ""),
                created_date=alphafold_source_info.get("Created Date", ""),
                protein_id=protein_instance,
            )
        except Exception as e:
            print(f"Error saving source info: {e}")


# Function to view results for multiple UniProt IDs
def view_user_input(request):
    if request.method == 'POST':
        protein_input = request.POST.get('protein_id')
        status, protein_id, protein_info = get_user_input(protein_input)          
        if status == "invalid_id":
            messages.error(request, protein_info)  
            return redirect('view_user_input')       
        if status == "no_data":
            messages.error(request, protein_info) 
            return redirect('view_user_input')       
        if status == "no_uniprot":
            messages.error(request, protein_info)
            return redirect('view_user_input')       
        if status == "multiple_ids":
            messages.error(request, f"Multiple UniProt IDs found for {protein_input}")
            return render(request, 'select_id.html', {
                'pdb_id': protein_input,
                'uniprot_ids': protein_id,                 
            })
        if status == "success":
            result = process_protein_data(protein_id, protein_info)
            database_data = DatabaseData(result)
            database_data.save_data()
            if 'error' in result:
                messages.error(request, result['error'])
            return render(request, 'results.html', result)   
    return render(request, 'user_input.html')


# Function to view results directly from the user input form 
def view_results(request, uniprot_id):
    status, protein_id, protein_info = get_user_input(uniprot_id)
    if status == "success":
        result = process_protein_data(protein_id, protein_info)
        database_data = DatabaseData(result)
        database_data.save_data()
        if "error" in result:
            messages.error(request, result["error"])
            return redirect("view_user_input") 
        return render(request, "results.html", result)
    else:
        messages.error(request, "Protein information not available or invalid input")
        return redirect("view_user_input") 
    