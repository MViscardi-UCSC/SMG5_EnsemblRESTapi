import requests, sys

server = "https://rest.ensembl.org"
ext = "/sequence/id"
headers = {"Content-Type": "application/json", "Accept": "application/json"}
r = requests.post(server + ext, headers=headers, data="""{ "ids" : ["ENSDNOG00000039862", "ENSMOCG00000016590" ] }""")

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
print(repr(decoded))