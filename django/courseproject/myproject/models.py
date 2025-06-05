from django.db import models

'''
Django's models.py is used to define the data models (tables) for the database backend supported by Django.
'''

class ProteinIdentifier(models.Model):
    id_number = models.AutoField(primary_key=True)
    uniprot_id = models.CharField(unique=True) 
    class Meta:
        db_table = 'protein_identifier'
        managed = False 
    def __str__(self):
        return self.uniprot_id


class PDB_IDS(models.Model):
    pdbids_id = models.AutoField(primary_key=True)
    uniprot_id = models.ForeignKey(
        ProteinIdentifier,
        related_name='pdb_ids',
        on_delete=models.CASCADE,
        db_column='uniprot_id', 
        to_field='uniprot_id',   
    )
    pdb_id = models.CharField(max_length=10)  
    class Meta:
        db_table = 'pdb_ids'
        managed = False
    def __str__(self):
        return f"{self.uniprot_id}, {self.pdb_id}"


class Protein(models.Model):
    protein_id = models.AutoField(primary_key=True)
    id_number = models.ForeignKey(
        ProteinIdentifier,
        related_name='proteins', 
        on_delete=models.CASCADE,
        db_column='id_number'
    )
    gene = models.CharField(max_length=255, blank=True, null=True)
    organism = models.CharField(max_length=255, blank=True, null=True)
    protein_name = models.CharField(max_length=255, blank=True, null=True)
    total_residues = models.IntegerField(default=0)  
    sequence = models.TextField(blank=True, null=True)
    class Meta:
        db_table = 'protein'
        managed = False
    def __str__(self):
        return f"Protein: {self.protein_name}, Gene: {self.gene}, Organism: {self.organism}"


class Residue(models.Model):
    residue_id = models.AutoField(primary_key=True)
    protein_id = models.ForeignKey(
        Protein,
        related_name='residues', 
        on_delete=models.CASCADE,
        db_column='protein_id'
    )
    chain_id = models.CharField(max_length=50)
    residue_name = models.CharField(max_length=50)
    residue_full_name = models.CharField(max_length=255)
    residue_pos = models.IntegerField()
    class Meta:
        db_table = 'residue'  
    def __str__(self):
        return f"Residue {self.residue_name} at position {self.residue_pos} in chain {self.chain_id}"


class Residue_Interactions(models.Model):
    contact_id = models.AutoField(primary_key=True)
    contact_type = models.CharField(max_length=100)
    res1_name = models.CharField(max_length=50)
    res1_pos = models.IntegerField()
    res1_coords = models.CharField(max_length=255)  # Store as a string (e.g., comma-separated values)
    res2_name = models.CharField(max_length=50)
    res2_pos = models.IntegerField()
    res2_coords = models.CharField(max_length=255)  # Store as a string (e.g., comma-separated values)
    distance = models.FloatField()
    protein_id = models.ForeignKey(
        Protein, 
        related_name='interactions',
        on_delete=models.CASCADE,
        db_column='protein_id')
    class Meta:
        db_table = 'residue_interactions'  
        managed = False
    def __str__(self):
        return f"Interaction between {self.res1_name} and {self.res2_name}"


class Source(models.Model):
    entry_id = models.AutoField(primary_key=True)  
    source_name = models.CharField(max_length=255)  # Name of the data source
    source_version = models.CharField(max_length=50)  # Version of the source or data entry
    created_date = models.DateField()  # Date the entry was created
    protein_id = models.ForeignKey(
        Protein, 
        related_name='source',
        on_delete=models.CASCADE,
        db_column='protein_id')
    class Meta:
        db_table = 'source'
        managed = False  
    def __str__(self):
        return f"{self.source_name}, {self.entry_version}"
