import requests

# Função para obter as vias KEGG
def get_kegg_pathways(gene_id):
    try:
        url = f"https://rest.kegg.jp/link/pathway/{gene_id}"
        response = requests.get(url)

        if response.status_code == 200:
            pathways = response.text.splitlines()
            return pathways
        else:
            return f"Error fetching KEGG pathways: {response.status_code}, {response.text}"
    except Exception as e:
        return f"Error fetching KEGG pathways: {str(e)}"


# Função para obter as vias Reactome
def get_reactome_pathways(gene_id):
    try:
        url = f"https://reactome.org/ContentService/data/pathways/{gene_id}"
        response = requests.get(url)

        if response.status_code == 200:
            pathways = response.json()
            return [pathway['displayName'] for pathway in pathways]
        else:
            return f"Error fetching Reactome pathways: {response.status_code}, {response.text}"
    except Exception as e:
        return f"Error fetching Reactome pathways: {str(e)}"


# Função para obter os detalhes de uma via KEGG
def get_kegg_pathway_details(pathway_id):
    try:
        url = f"https://rest.kegg.jp/get/{pathway_id}"
        response = requests.get(url)

        if response.status_code == 200:
            return response.text
        else:
            return f"Error fetching KEGG pathway details: {response.status_code}, {response.text}"
    except Exception as e:
        return f"Error fetching KEGG pathway details: {str(e)}"


# Função para obter os detalhes de uma via Reactome
def get_reactome_pathway_details(pathway_name):
    try:
        url = f"https://reactome.org/ContentService/data/pathways/{pathway_name}"
        response = requests.get(url)

        if response.status_code == 200:
            return response.json()
        else:
            return f"Error fetching Reactome pathway details: {response.status_code}, {response.text}"
    except Exception as e:
        return f"Error fetching Reactome pathway details: {str(e)}"
