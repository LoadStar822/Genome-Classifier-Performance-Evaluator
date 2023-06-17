# coding:utf-8
"""
Author  : Tian
Time    : 2023-04-22 13:29
Desc    :This script downloads the FASTA sequences of a list of organisms.
"""
import os
from Bio import Entrez, SeqIO
from io import StringIO


Entrez.email = ''
Entrez.api_key = ''

species_list = [
    'Escherichia coli', 'Staphylococcus aureus', 'Streptococcus pyogenes',
    'Bacillus subtilis', 'Pseudomonas aeruginosa', 'Mycobacterium tuberculosis',
    'Lactobacillus acidophilus', 'Clostridium difficile', 'Helicobacter pylori',
    'Salmonella enterica', 'Vibrio cholerae', 'Shigella flexneri',
    'Campylobacter jejuni', 'Listeria monocytogenes', 'Neisseria meningitidis',
    'Haemophilus influenzae', 'Corynebacterium diphtheriae', 'Legionella pneumophila',
    'Enterococcus faecalis', 'Chlamydia trachomatis', 'Methanocaldococcus jannaschii',
    'Methanosarcina mazei', 'Pyrococcus furiosus', 'Sulfolobus solfataricus',
    'Halobacterium salinarum', 'Thermococcus kodakarensis', 'Methanobrevibacter smithii',
    'Nanoarchaeum equitans', 'Aeropyrum pernix', 'Methanococcus maripaludis'
]

output_folder = 'fasta_sequences'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for species in species_list:
    print(f'Downloading {species}...')
    search_query = f'"{species}"[Organism] AND "complete genome"[title]'
    search_results = Entrez.read(Entrez.esearch(db='nucleotide', term=search_query, retmax=1))

    if search_results['Count'] != '0':
        sequence_id = search_results['IdList'][0]
        record = Entrez.efetch(db='nucleotide', id=sequence_id, rettype='fasta', retmode='text')
        fasta = record.read()

        # Extract strain or variant information
        fasta_parsed = SeqIO.read(StringIO(fasta), "fasta")
        description_parts = fasta_parsed.description.split(' ')
        strain_or_variant = ' '.join(description_parts[1:]) if len(description_parts) > 1 else 'Unknown'
        print(f'Downloaded strain or variant: {strain_or_variant}')

        output_file = os.path.join(output_folder, f'{species.replace(" ", "_")}.fasta')
        with open(output_file, 'w') as f:
            f.write(fasta)
        print(f'Successfully downloaded {species}.')
    else:
        print(f'No complete genome found for {species}.')