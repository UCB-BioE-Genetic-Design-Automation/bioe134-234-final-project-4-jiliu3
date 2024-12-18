def reverse_complement(sequence):
    """
    Calculates the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The reverse complement of the DNA sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters.
    """
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def translate(sequence):
    """
    Translates a DNA sequence into a protein sequence based on the standard genetic code.

    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The corresponding protein sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters or is not a multiple of three.
    """
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")
    if len(sequence) % 3 != 0:
        raise ValueError("Length of DNA sequence is not a multiple of three, which is required for translation.")

    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        protein += codon_table.get(codon, '_')  # Using '_' for unknown or stop codons
    return protein


if __name__ == "__main__":
    # Example DNA sequence for demonstration
    dna_example = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    try:
        # Perform reverse complement
        rc_result = reverse_complement(dna_example)
        print(f"Reverse Complement of '{dna_example}': {rc_result}")

        # Perform translation
        translation_result = translate(dna_example)
        print(f"Translation of '{dna_example}': {translation_result}")
    except Exception as e:
        print(f"Error: {str(e)}")
