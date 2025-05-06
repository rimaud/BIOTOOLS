def translate_dna(dna_seq):
    try:
        # Garantir que a sequência contém apenas A, T, C, G
        dna_seq = dna_seq.replace(" ", "").replace("\n", "")  # Remover espaços e quebras de linha
        valid_bases = {"A", "T", "C", "G"}
        if not all(base in valid_bases for base in dna_seq):
            raise ValueError("Sequência de DNA inválida!")
        
        # Adicionar "N" se o comprimento da sequência não for múltiplo de 3
        if len(dna_seq) % 3 != 0:
            padding_length = 3 - (len(dna_seq) % 3)
            dna_seq += "N" * padding_length
            st.warning(f"A sequência foi preenchida com {padding_length} 'N' para garantir um comprimento múltiplo de 3.")

        # Traduzir a sequência de DNA para proteína
        protein = str(Seq(dna_seq).translate())
        return protein
    except Exception as e:
        return f"Erro: {str(e)}"
