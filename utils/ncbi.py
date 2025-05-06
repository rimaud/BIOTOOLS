from Bio import Entrez, SeqIO

def fetch_sequence(accession_id, db_type):
    handle = Entrez.efetch(db=db_type, id=accession_id, rettype="gb")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record