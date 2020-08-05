"""
EnsemblREST.py
Marcus Viscardi     Aug 5, 2020

Trying to pull everything together into one script to make later reuse easier
"""

import requests, sys, json
from pprint import pprint

SERVER = "https://rest.ensembl.org"


def fetch_endpoint(server, request, content_type, options: dict = None):
    """
    Fetch an endpoint from the server, allow overriding of default content-type
    """
    if options:
        r = requests.get(server+request,
                         headers={"Accept": content_type},
                         json=options)
    else:
        r = requests.get(server+request,
                         headers={"Accept": content_type})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text


def fetch_endpoint_POST(server, request, data, content_type='application/json', options: dict = None):

    if options:
        r = requests.post(server+request,
                          headers={"Content-Type": content_type},
                          data=data, json=options)
    else:
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
    