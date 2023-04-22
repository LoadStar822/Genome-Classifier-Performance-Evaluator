# coding:utf-8
"""
Author  : Tian
Time    : 2023-04-22 13:29
Desc    :This script downloads the FASTA sequences of a list of organisms.
"""
from Bio import Entrez

# Set your email address to be used in the Entrez API
Entrez.email = "your.email@example.com"

# Define a dictionary of organisms and their associated taxonomy group (bacteria or archaea)
organisms = {
    "Escherichia coli": "bacteria",
    "Staphylococcus aureus": "bacteria",
    "Streptococcus pyogenes": "bacteria",
    "Bacillus subtilis": "bacteria",
    "Pseudomonas aeruginosa": "bacteria",
    "Mycobacterium tuberculosis": "bacteria",
    "Lactobacillus acidophilus": "bacteria",
    "Clostridium difficile": "bacteria",
    "Helicobacter pylori": "bacteria",
    "Salmonella enterica": "bacteria",
    "Vibrio cholerae": "bacteria",
    "Shigella flexneri": "bacteria",
    "Campylobacter jejuni": "bacteria",
    "Listeria monocytogenes": "bacteria",
    "Neisseria meningitidis": "bacteria",
    "Haemophilus influenzae": "bacteria",
    "Corynebacterium diphtheriae": "bacteria",
    "Legionella pneumophila": "bacteria",
    "Enterococcus faecalis": "bacteria",
    "Chlamydia trachomatis": "bacteria",
    "Methanocaldococcus jannaschii": "archaea",
    "Methanosarcina mazei": "archaea",
    "Pyrococcus furiosus": "archaea",
    "Sulfolobus solfataricus": "archaea",
    "Halobacterium salinarum": "archaea",
    "Thermococcus kodakarensis": "archaea",
    "Methanobrevibacter smithii": "archaea",
    "Nanoarchaeum equitans": "archaea",
    "Aeropyrum pernix": "archaea",
    "Methanococcus maripaludis": "archaea"
}


# Define a function to search for an organism in the given database
def search_organism(db, organism):
    term = f"{organism}[Organism]"
    handle = Entrez.esearch(db=db, term=term, retmax=1)
    return Entrez.read(handle)


# Define a function to fetch a sequence from the database using its ID
def fetch_sequence(db, seq_id):
    handle = Entrez.efetch(db=db, id=seq_id, rettype="fasta", retmode="text")
    return handle.read()


# Define a function to save a FASTA sequence to a file
def save_fasta_file(organism, fasta_seq):
    with open(f"{organism}.fasta", "w") as f:
        f.write(fasta_seq)


# Define a function to download sequences for a list of organisms
def download_organism_sequences(organisms):
    for organism, db in organisms.items():
        record = search_organism(db, organism)
        if record["IdList"]:
            seq_id = record["IdList"][0]
            fasta_seq = fetch_sequence(db, seq_id)
            save_fasta_file(organism, fasta_seq)
            print(f"Downloaded {organism}")
        else:
            print(f"No record found for {organism}")


download_organism_sequences(organisms)
