import ViennaRNA


def calculate_mfe(sequence):
    try:
        rna_structure = ViennaRNA.fold(sequence)
    except:
        return "error"
    return rna_structure[1]
