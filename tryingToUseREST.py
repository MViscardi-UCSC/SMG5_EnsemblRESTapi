"""
tryingToUseREST.py
Marcus Viscardi     Aug 5, 2020

I am just going to try and use a RESTapi to do SOMETHING/ANYTHING
"""

import requests, sys

server = "https://rest.ensembl.org"
ext = "/sequence/id/ENSG00000198952?"

if __name__ == '__main__':
    r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    print(r.text)
