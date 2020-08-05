import requests, sys

server = "https://rest.ensembl.org"
ext = "/sequence/id?type=cds"
headers = {"Content-Type": "application/json", "Accept": "application/json"}
r = requests.post(server + ext, headers=headers, data="""{ "ids" : ["ENSG00000198952" ] }""")

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
print(repr(decoded))