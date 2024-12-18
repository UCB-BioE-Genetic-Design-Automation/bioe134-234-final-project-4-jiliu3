from .polynucleotide import Polynucleotide, dsDNA
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio.Restriction import RestrictionBatch, AllEnzymes
from pydna.utils import rc
from Bio.Restriction import Restriction
import math
from construction_file import Step, ConstructionFile

def get_restriction_enzyme(enzyme_name):
    """
    Retrieve a restriction enzyme by name.
    """
    if enzyme_name in AllEnzymes:
        return Restriction.__dict__.get(enzyme_name, None)
    else:
        raise ValueError(f"Enzyme '{enzyme_name}' is not recognized.")

symbol_to_gene = {
    "lacz": "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC",
    "gfp": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCAGATGAACTTCAGGGTC",
    # Add other gene symbols and their sequences as needed
}

reagent_to_price = {
    # Enzyme prices (per unit)
    **{enzyme: 10 for enzyme in ['AarI', 'ApaI', 'BamHI', 'BbsI', 'BglII', 'BsaI', 'BseRI', 'BsmBI', 'ClaI', 'EcoRI', 'EcoRV', 'HindIII', 'I-CreI', 'I-SceI', 'KpnI', 'NcoI', 'NdeI', 'NotI', 'PstI', 'SacI', 'SalI', 'SapI', 'SmaI', 'SpeI', 'XbaI', 'XhoI']},
    
    # Antibiotic prices (per unit)
    'G418': 5, 'Hygro': 4, 'Nat': 6, 'Zeo': 8,

    # Other reagents
    'Taq Polymerase': 25,
    'dNTPs': 15,
    'DNA Ligase': 20,
    'Gibson Assembly Mix': 40,
    'Competent Cells': 50
}

def round_to_nearest_15(minutes):
    return math.ceil(minutes / 15) * 15

def simulate(constructionFile: ConstructionFile): 
    time = 0
    cost = 0
    for step in constructionFile:
        if step.operation == 'Digest':
            result = perform_digest(step)
        elif step.operation == 'PCR':
            result = perform_pcr(step)
        elif step.operation == 'Ligate':
            result = perform_ligate(step)
        elif step.operation == 'GoldenGate':
            result = perform_goldengate(step)
        elif step.operation == 'Gibson':
            result = perform_gibson(step)
        elif step.operation == 'Transform':
            result = perform_transform(step)
        time += result[2]
        cost +=result[1]
        symbol_to_gene[step.output] = result  
    return symbol_to_gene, time, cost

def resolve_gene_input(gene_input, dictionary):
    """
    Resolves a gene input into its full sequence.
    If the input is a gene symbol (lowercase), it looks up the sequence in the dictionary.
    If the input is a sequence (uppercase), it returns it as-is.
    """
    if gene_input.islower():
        if gene_input in dictionary:
            return dictionary[gene_input]
        else:
            raise ValueError(f"Gene symbol '{gene_input}' not found in symbol_to_gene dictionary.")
    elif gene_input.isupper():
        return gene_input
    else:
        raise ValueError(f"Invalid gene input '{gene_input}'. It must be either uppercase (sequence) or lowercase (symbol).")


def perform_pcr(pcr_step):
    """
    Perform a PCR operation with cost and time predictions.
    """
    forward_primer = resolve_gene_input(pcr_step.forward_oligo, symbol_to_gene)
    reverse_primer = resolve_gene_input(pcr_step.reverse_oligo, symbol_to_gene)
    template = resolve_gene_input(pcr_step.template, symbol_to_gene)
    output = pcr_step.output

    # Simulate PCR product creation
    product_sequence = forward_primer + template + reverse_primer
    product = Polynucleotide(
        product_sequence,
        ext5="",
        ext3="",
        is_double_stranded=True,
        is_circular=False,
        mod_ext5="phosphate",
        mod_ext3="phosphate"
    )

    # Predict cost
    cost = reagent_to_price['Taq Polymerase'] + reagent_to_price['dNTPs']

    # Predict time based on template length
    template_length = len(template)
    time = round_to_nearest_15(30 + template_length / 500 * 15)  # 30 minutes base, +15 per 500 bp

    return product, cost, time


from Bio.Restriction import RestrictionBatch, AllEnzymes

def perform_digest(digest_step):
    """
    Perform a Digest operation with cost and time predictions.
    """
    dna_sequence = resolve_gene_input(digest_step.dna, symbol_to_gene)
    enzymes = digest_step.enzymes
    frag_select = digest_step.fragSelect
    output = digest_step.output

    # Convert DNA to Dseqrecord
    dna_dseq = Dseqrecord(Dseq(dna_sequence))

    # Create enzyme batch and perform digestion
    enzyme_batch = RestrictionBatch([enz for enz in enzymes if enz in AllEnzymes])
    fragments = dna_dseq.cut(enzyme_batch)

    # Select fragment
    if 0 <= frag_select < len(fragments):
        selected_fragment = fragments[frag_select]
        product = Polynucleotide(
            str(selected_fragment.seq.watson),
            ext5="",
            ext3="",
            is_double_stranded=True,
            is_circular=False,
            mod_ext5="hydroxyl",
            mod_ext3="hydroxyl"
        )
    else:
        raise ValueError("Fragment selection index out of range.")

    # Predict cost
    cost = sum(reagent_to_price[enzyme] for enzyme in enzymes)

    # Predict time based on DNA length and number of enzymes
    dna_length = len(dna_sequence)
    time = round_to_nearest_15(15 + len(enzymes) * 10 + dna_length / 1000 * 10)  # Base time + enzyme + length

    return product, cost, time


def perform_gibson(gibson_step):
    """
    Perform a Gibson assembly operation with cost and time predictions.
    """
    fragments = [resolve_gene_input(fragment, symbol_to_gene) for fragment in gibson_step.dnas]
    output = gibson_step.output

    # Simulate Gibson Assembly
    gibson_product = ''.join(fragments)
    product = Polynucleotide(
        gibson_product,
        ext5="",
        ext3="",
        is_double_stranded=True,
        is_circular=False,
        mod_ext5="hydroxyl",
        mod_ext3="hydroxyl"
    )

    # Predict cost
    cost = reagent_to_price['Gibson Assembly Mix']

    # Predict time based on total length of fragments
    total_length = sum(len(frag) for frag in fragments)
    time = round_to_nearest_15(30 + total_length / 1000 * 15)  # Base time + 15 per 1 kb

    return product, cost, time


def perform_goldengate(gg_step):
    """
    Perform a Golden Gate operation with cost and time predictions.
    """
    # Resolve the fragments
    fragments = [resolve_gene_input(fragment, symbol_to_gene) for fragment in gg_step.dnas]
    enzyme_name = gg_step.enzyme
    output = gg_step.output

    # Convert DNA sequences to Dseqrecords
    dna_dseqs = [Dseqrecord(Dseq(dsDNA(seq).sequence)) for seq in fragments]

    # Restriction enzyme
    enzyme = get_restriction_enzyme(enzyme_name)
    if enzyme is None:
        raise ValueError(f"Enzyme '{enzyme_name}' not found.")

    print(f"Using enzyme: {enzyme_name}")

    # Simulate digestion
    digested_fragments = []
    for dna in dna_dseqs:
        # Check the fragment before digestion
        print(f"Original fragment: {str(dna.seq.watson)}")

        # Perform the cut (digestion)
        digested = dna.cut(enzyme)

        # If no digestion occurs, just add the fragment as is
        if not digested:
            print(f"No cut found for enzyme {enzyme_name} on fragment: {str(dna.seq.watson)}")
            digested_fragments.append(dna)
        else:
            # Check the digested fragments
            print(f"Digested fragments: {[str(fragment.seq.watson) for fragment in digested]}")
            digested_fragments.extend(digested)

    if not digested_fragments:
        raise ValueError("Digested fragments are empty, digestion failed.")

    # Ligate fragments: Combine the fragments to form the final ligated product
    ligated_sequence = ''.join(str(fragment.seq.watson) for fragment in digested_fragments)

    # Ensure the overhang (restriction site) is removed
    site_sequence = enzyme.site
    ligated_sequence = ligated_sequence.replace(site_sequence, "")

    # Create the final product sequence
    product = Polynucleotide(
        ligated_sequence,
        ext5="",
        ext3="",
        is_double_stranded=True,
        is_circular=False,
        mod_ext5="phosphate",
        mod_ext3="phosphate"
    )

    # Predict cost
    cost = reagent_to_price[enzyme_name] + reagent_to_price['DNA Ligase']

    # Predict time (Golden Gate typically takes around 30 minutes)
    time = round_to_nearest_15(30)  # Golden Gate time is fairly consistent, typically 30 minutes

    return product, cost, time





def perform_transform(transform_step):
    """
    Perform a Transformation operation with cost and time predictions.
    """
    plasmid = resolve_gene_input(transform_step.dna, symbol_to_gene)
    host = transform_step.strain
    antibiotics = transform_step.antibiotics
    temperature = transform_step.temperature
    output = transform_step.output

    # Simulate transformation
    transformation_result = f"Transformation successful: {output} now carries {plasmid} in host {host} at {temperature}°C."

    # Predict cost
    cost = reagent_to_price['Competent Cells'] + sum(reagent_to_price[antibiotic] for antibiotic in antibiotics)

    # Predict time
    time = round_to_nearest_15(45)  # Assume 45 minutes for transformation

    return transformation_result, cost, time

def perform_ligate(ligate_step):
    """
    Perform a Ligate operation with cost and time predictions.
    """
    fragments = [resolve_gene_input(fragment, symbol_to_gene) for fragment in ligate_step.dnas]
    output = ligate_step.output

    # Simulate ligation - concatenate sequences
    ligated_product = ''.join(fragments)
    product = Polynucleotide(
        ligated_product,
        ext5="",
        ext3="",
        is_double_stranded=True,
        is_circular=False,
        mod_ext5="phosphate",
        mod_ext3="phosphate"
    )

    # Predict cost
    cost = reagent_to_price['DNA Ligase']

    # Predict time based on total fragment length
    total_length = sum(len(fragment) for fragment in fragments)
    time = round_to_nearest_15(30 + total_length / 1000 * 10)  # Base time + 10 min per 1 kb

    return product, cost, time
