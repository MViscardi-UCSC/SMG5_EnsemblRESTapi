"""
pullFamilyMembers.py
Marcus Viscardi     Aug 5, 2020

Pulling members of a protein family with GET family/id/:id
"""

import requests, sys, json
from pprint import pprint

protein_family = "PTHR15696_SF1"

server = "https://rest.ensembl.org"
ext = f"/family/id/{protein_family}?sequence=cdna&aligned=0&member_source=ensembl"

if __name__ == '__main__':
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    pprint(repr(decoded))
    with open(f"200805_Ensembl_{protein_family}.json", 'w') as write_to:
        json.dump(r.json(), write_to)
