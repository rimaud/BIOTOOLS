import requests

def fetch_sequence(uniprot_id):
    response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta")
    return response.text if response.ok else None