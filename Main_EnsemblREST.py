"""
Main_EnsemblREST.py
Marcus Viscardi     Aug 5, 2020

Trying to pull everything together into one script to make later reuse easier
"""

import requests, sys, json
from pprint import pprint
from time import strftime

SERVER = "https://rest.ensembl.org"


def fetch_endpoint(server, request, content_type='application/json', params={}):
    """
    Fetch an endpoint from the server, allow overriding of default content-type
    """
    url = server+request
    headers={
        "Accept": content_type
    }
    r = requests.get(url, headers=headers, params=params)

    if not r.ok:
        r.raise_for_status()
        raise Exception("Endpoint failed.")

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text


def fetch_endpoint_POST(server, request, data, content_type='application/json', params={}):
    url = server + request
    headers = {
        "Content-Type": content_type
    }
    r = requests.post(url, headers=headers, json=data, params=params)

    if not r.ok:
        r.raise_for_status()
        raise Exception("Posting to endpoint failed.")

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text


def proteinFamilyParentSequences(protein_family: str, parent_seq_type: str = 'cDNA') -> str:
    if parent_seq_type.lower() == "cdna":
        parent_seq_type = "cDNA"
        seq_type = "cdna"
    elif parent_seq_type.lower() == "cds":
        parent_seq_type = "CDS"
        seq_type = "cds"
    else:
        raise TypeError
    
    # First pull all members of the protein family
    fam_params = {
        'sequence': 'none',
        'aligned': '0',
        'member_source': 'ensembl'
    }
    protein_family_dict = fetch_endpoint(SERVER, f"/family/id/{protein_family}", params=fam_params)
    
    # Isolate just the gene IDs from here
    members_geneIDs = set(protein['gene_stable_id'] for protein in protein_family_dict['members'])
    
    # Convert back and forth to a set to remove redundant genes, these will show up again when looking for transcripts
    #   (I think...)
    members_geneIDs = list(members_geneIDs)
    
    # Set output file name and sequence call parameters
    fasta_output = f"{strftime('%y%m%d')}_Ensembl_{protein_family}_{parent_seq_type}.fasta"
    seq_params = {
        'type': seq_type,
    }

    # The maximum post size seems to be 50 for Ensembl, iterate through the IDs in sets of 49!
    for i in range(0, len(members_geneIDs), 50):
        try:
            data = {"ids": members_geneIDs[i:i + 49]}
            print(f"Running geneIDs {i} - {i + 49}")
        except IndexError:
            data = {"ids": members_geneIDs[i:]}
            print(f"Running geneIDs {i} - {len(members_geneIDs)}")
        
        family_seqs_list = fetch_endpoint_POST(SERVER, "/sequence/id", data, params=seq_params)

        print(family_seqs_list)
        with open(fasta_output, 'a', encoding='utf-8') as f:  # Append to the file with each loop
            for hit in family_seqs_list:
                f.write(f">{hit['query']}({hit['id']})->{parent_seq_type}\n")
                f.write(f"{hit['seq']}\n")
    return fasta_output


if __name__ == '__main__':
    # proteinFamilyParentSequences("PTHR15696_SF1", cds_or_cdna='CDS')
    
    protein_fam_params = {
        'sequence': 'protein',
        'aligned': '1',
        # 'member_source': 'ensembl',
    }
    protein_family = "PTHR15696_SF1"
    
    protein_family_dictionary = fetch_endpoint(SERVER, f"/family/id/{protein_family}", params=protein_fam_params)
    pprint(protein_family_dictionary)
    with open(f"./{strftime('%y%m%d')}_{protein_family}_proteins.json", 'w') as f:
        json.dump(protein_family_dictionary, f)
    fasta_style_dict = {f"{hit['protein_stable_id']}({hit['source_name']})->Protein": hit['protein_alignment']
                        for hit in protein_family_dictionary['members']}
    with open(f"./{strftime('%y%m%d')}_{protein_family}_proteins.fasta", 'w') as fasta_out:
        for name, seq in fasta_style_dict.items():
            print(name)
            fasta_out.write(f">{name}\n{seq}\n")
