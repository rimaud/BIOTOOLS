import streamlit as st
import streamlit.components.v1 as components
# Standard library
import requests
import pandas as pd
from io import StringIO 

# Third-party
from Bio.Seq import Seq
from Bio import SeqIO
import py3Dmol
from Bio import Align
aligner = Align.PairwiseAligner()

# Local modules
from utils.uniprot import fetch_sequence as fetch_uniprot
from utils.pathways import get_kegg_pathways, get_reactome_pathways, get_kegg_pathway_details, get_reactome_pathway_details
from utils.alignment import plot_alignment
from utils.translate import translate_dna   

from Bio import Entrez
import streamlit as st

# Set Entrez email
from Bio import Entrez
Entrez.email = "seu.email@exemplo.com"


# Config
st.set_page_config(layout="wide", page_title="BioFacíl", page_icon="🧬")

# CSS
with open("assets/style.css") as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

# Sidebar com seletor único
st.sidebar.header("🧬 BioFacíl")
tool = st.sidebar.selectbox(
    "Tools",
    ["Sequence Fetch", "DNA Translation", "Sequence Alignment", "Pathway Analysis", "3D Structure"],
    key="tool_selector"
)

# Main App
if tool == "Sequence Fetch":
    st.header("🔍 Fetch Biological Sequences")
    col1, col2 = st.columns(2)

    with col1:
        db_type = st.radio("Database", ["NCBI Nucleotide", "NCBI Protein", "UniProt"])
        search_input = st.text_input("Enter ID or Gene Name (e.g. NM_001301717, TP53, P12345):")
        species = st.selectbox("Species", ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Escherichia coli", "Saccharomyces cerevisiae"])

    if st.button("Fetch"):
        try:
            if db_type.startswith("NCBI"):
                db = "nucleotide" if "Nucleotide" in db_type else "protein"

                # Primeiro: buscar via esearch
                query = f"{search_input}[Gene] AND {species}[Organism]"
                search_handle = Entrez.esearch(db=db, term=query, retmax=1)
                search_results = Entrez.read(search_handle)

                if not search_results["IdList"]:
                    st.error("No results found.")
                else:
                    seq_id = search_results["IdList"][0]

                    fetch_handle = Entrez.efetch(db=db, id=seq_id, rettype="fasta", retmode="text")
                    fasta_data = fetch_handle.read()
                    record = SeqIO.read(StringIO(fasta_data), "fasta")  # Aqui é necessário StringIO

                    st.success(f"Found: {record.description}")
                    st.code(str(record.seq))
                    st.download_button("Download FASTA", f">{record.id}\n{record.seq}", f"{search_input}.fasta")

            else:
                # UniProt
                seq = fetch_uniprot(search_input)
                st.code(seq)

        except Exception as e:
            st.error(f"Failed to fetch sequence: {str(e)}")

elif tool == "DNA Translation":
    st.header("🧬 DNA → Protein Translator")
    
    # Corrigindo a indentação
    dna_seq_input = st.text_area("Paste DNA sequence:")

    if dna_seq_input:  # Checa se a entrada não está vazia
        translated_protein = translate_dna(dna_seq_input)  # Usando a função de tradução
        if "Erro" not in translated_protein:  # Se não for um erro
            st.success(f"Protein sequence: {translated_protein}")
        else:
            st.error(translated_protein)  # Exibe o erro, se houver

elif tool == "Sequence Alignment":
    st.header("🔄 Sequence Alignment")
    seq1 = st.text_area("Sequence 1:", "ATGCGATACGT")
    seq2 = st.text_area("Sequence 2:", "ATGCGTTACGT")
    
    if st.button("Align"):
        with st.spinner("Aligning sequences..."):
            try:
                fig = plot_alignment(seq1, seq2)
                st.pyplot(fig)
            except Exception as e:
                st.error(f"Alignment failed: {str(e)}")

if tool == "Pathway Analysis":
    st.header("🔍 Pathway Analysis")
    
    # Entrada para ID ou nome da via
    pathway_input = st.text_input("Enter Pathway ID or Name (e.g. hsa04010 or MAPK signaling pathway):")
    gene_input = st.text_input("Enter Gene ID or Symbol (e.g. TP53, BRCA1):")

    # Se o usuário inserir o nome ou ID da via
    if pathway_input:
        st.subheader(f"KEGG Pathway Details for {pathway_input}")
        kegg_details = get_kegg_pathway_details(pathway_input)
        st.text(kegg_details)  # Exibe os detalhes da via no KEGG

        st.subheader(f"Reactome Pathway Details for {pathway_input}")
        reactome_details = get_reactome_pathway_details(pathway_input)
        st.text(reactome_details)  # Exibe os detalhes da via no Reactome

    # Se o usuário inserir um gene
    if gene_input:
        st.subheader(f"KEGG Pathways for Gene {gene_input}")
        kegg_pathways = get_kegg_pathways(gene_input)
        if isinstance(kegg_pathways, list) and kegg_pathways:
            for pathway in kegg_pathways:
                st.write(pathway)
        else:
            st.error(kegg_pathways)  # Caso não tenha encontrado pathways KEGG

        st.subheader(f"Reactome Pathways for Gene {gene_input}")
        reactome_pathways = get_reactome_pathways(gene_input)
        if isinstance(reactome_pathways, list) and reactome_pathways:
            for pathway in reactome_pathways:
                st.write(pathway)
        else:
            st.error(reactome_pathways)  # Caso não tenha encontrado pathways Reactome
elif tool == "3D Structure":
    st.header("🧊 Protein 3D Structure")
    pdb_id = st.text_input("PDB ID (e.g. 1A2C):")
    if pdb_id:
        try:
            view = py3Dmol.view(query=f"pdb:{pdb_id}")
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo()
            components.html(view._make_html(), height=500)
        except Exception as e:
            st.error(f"Failed to load 3D structure: {str(e)}")
