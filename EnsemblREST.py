"""
EnsemblREST.py
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
                      json={"type": "protein"},
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
    
    # The maximum post size is 50, so I'll start with just grabbing the first 50, eventually I'll do it all in pieces
    data = f"""{{ "ids" : {members_geneIDs[:10]} }}"""
    data = data.replace("'", '"')
    print(data)
    family_seqs_dict = fetch_endpoint_POST(SERVER, "/sequence/id", data)
    print(family_seqs_dict)
