"""
ensemblJSONtoFASTA.py
Marcus Viscardi     Aug 5, 2020

parse the output of REST GET family/id/:id to a fasta format (mostly to feed to UGENE)
"""

import json

json_file = "200805_Ensembl_PTHR15696_SF1.json"
fasta_file = "200805_Ensembl_PTHR15696_SF1.fasta"

if __name__ == '__main__':
    with open(json_file, 'r') as json_file:
        prot_family_dict = json.load(json_file)
    
    with open(fasta_file, 'w', encoding='utf-8') as fasta_file:
        for protein in prot_family_dict['members']:
            # print(protein)
            fasta_file.write(f">{protein['gene_stable_id']}(proteinID:{protein['protein_stable_id']})"
                             f"—{protein['genome']}\n")
            fasta_file.write(f"{protein['protein_seq']}\n")
