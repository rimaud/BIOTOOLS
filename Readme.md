

# BioFacíl

**BioFacíl** is a user-friendly Streamlit web application designed to simplify access to biological databases such as NCBI, UniProt, and KEGG/Reactome. It allows users to retrieve and visualize biological data like gene/protein sequences and pathway information without writing code.

## 🚀 Features

- 🔍 Search and fetch sequences from **NCBI**
- 🧬 Query proteins and metadata from **UniProt**
- 🧪 Explore metabolic pathways via **KEGG** and **Reactome**
- 📁 Upload and preview FASTA files
- 📊 Visualize pathway data in interactive tables and charts

## 📦 Installation

1. Clone the repository

2. Create a virtual environment (optional but recommended)

3. Install the dependencies:
pip install -r requirements.txt

4. Run the app
streamlit run app.py


# *Project Structure*

``` bioFacíl/                  # Root folder
      │
      ├── app.py                     # Main Streamlit application
      ├── requirements.txt           # Python dependencies
      ├── secrets.toml               # API keys (optional, for deployment)
      │
      ├── assets/                    # Static files (images, CSS)
      │   ├── logo.png               # App logo
      │   └── style.css              # Custom styles (optional)
      │
      ├── data/                      # Preprocessed data (examples)
      │   ├── example.fasta          # Example FASTA file
      │   └── pathways.csv           # Metabolic pathway data
      │
      ├── utils/                     # Helper modules
      │   ├── ncbi.py                # Functions for NCBI search
      │   ├── uniprot.py             # Functions for UniProt
      │   └── pathways.py            # Integration with KEGG/Reactome
      │
      └── README.md                  # Project documentation
```