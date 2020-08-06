"""
Main_EnsemblREST.py
Marcus Viscardi     Aug 5, 2020

Trying to pull everything together into one script to make later reuse easier
"""

import requests, sys, json
from pprint import pprint

SERVER = "https://rest.ensembl.org"


def fetch_endpoint(server, request, content_type):
    """
    Fetch an endpoint from the server, allow overriding of default content-type
    """
    r = requests.get(server + request,
                     headers={"Accept": content_type})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text


def fetch_endpoint_POST(server, request, data, content_type='application/json'):
    r = requests.post(server + request,
                      headers={"Content-Type": content_type},
                      data=data)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text


if __name__ == '__main__':
    # First pull all members of the protein family
    prot_family_dict = fetch_endpoint(SERVER,
                                 "/family/id/PTHR15696_SF1?sequence=none&aligned=0&member_source=ensembl",
                                 'application/json')
    
    # Isolate just the gene IDs from here
    members_geneIDs = []
    for protein in prot_family_dict['members']:
        members_geneIDs.append(protein['gene_stable_id'])
    
    # Convert back and forth to a set to remove redundant genes, these will show up again when looking for transcripts
    #   (I think...)
    members_geneIDs = list(set(members_geneIDs))
    
    # The maximum post size seems to be 50 for Ensembl, iterate through the IDs in sets of 49!
    for i in range(0, len(members_geneIDs), 50):
        try:
            data = f"""{{ "ids" : {members_geneIDs[i:i+49]} }}"""
            print(f"Running geneIDs {i} - {i+49}")
        except IndexError:  # Does the slice function cause an index error?
            data = f"""{{ "ids" : {members_geneIDs[i:]} }}"""
            print(f"Running geneIDs {i} - {len(members_geneIDs)}")
        data = data.replace("'", '"')
        # print(data)
        family_seqs_list = fetch_endpoint_POST(SERVER, "/sequence/id?type=cds", data)
        # print(family_seqs_list)
        fasta_output = "200805_Ensembl_PTHR15696_SF1_CDS.fasta"
        with open(fasta_output, 'a', encoding='utf-8') as f:  # Append to the file with each loop
            for hit in family_seqs_list:
                f.write(f">{hit['query']}({hit['id']})->CDS\n")
                f.write(f"{hit['seq']}\n")
