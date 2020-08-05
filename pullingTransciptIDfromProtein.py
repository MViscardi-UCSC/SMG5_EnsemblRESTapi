"""
pullingTranscriptIDfromProtein.py
Marcus Viscardi     Aug 8, 2020

Goal is to retrieve the Transcript or Gene ID by using a given protein ID and the "Parent" param
"""

import requests, sys
from pprint import pprint

if __name__ == '__main__':
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.post(server + ext, headers=headers, data='{ "ids" : ["PTHR15696_SF1", "ENSHCOP00000022461" ] }')
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    print(repr(decoded))
