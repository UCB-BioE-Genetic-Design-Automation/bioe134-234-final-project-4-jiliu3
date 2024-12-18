import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pytest
from cf_simulator.construction_file import PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from cf_simulator.construction_file_simulator import perform_pcr, perform_digest, perform_ligate, perform_gibson, perform_goldengate

# Mock the symbol_to_gene dictionary
symbol_to_gene = {
    "lacz": "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC",
    "gfp": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCAGATGAACTTCAGGGTC",
}

### PCR Tests
def test_pcr_with_short_template():
    # Template with short sequence
    forward_primer = "ATGC"
    reverse_primer = "CGTA"
    template = "GATC"
    pcr_step = PCR(forward_primer, reverse_primer, template, "PCR_Product")
    product, _, _ = perform_pcr(pcr_step)
    assert product.sequence == "ATGCGATCCGTA", "PCR sequence mismatch for short template"

def test_pcr_with_symbol_template():
    # Using gene symbol as template
    forward_primer = "ATGC"
    reverse_primer = "CGTA"
    template = "lacz"
    pcr_step = PCR(forward_primer, reverse_primer, template, "PCR_Product")
    product, _, _ = perform_pcr(pcr_step)
    expected_sequence = "ATGC" + symbol_to_gene["lacz"] + "CGTA"
    assert product.sequence == expected_sequence, "PCR sequence mismatch for symbol template"


### Digest Tests
def test_digest_single_cut():
    # Simple digest with one EcoRI site
    dna_sequence = "ATGAATTCTAG"
    enzymes = ["EcoRI"]
    digest_step = Digest(dna_sequence, enzymes, 0, "Digest_Product")
    product, _, _ = perform_digest(digest_step)
    expected_sequence = "ATG"
    assert product.sequence == expected_sequence, "Digest sequence mismatch for single EcoRI cut"

### Ligate Tests
def test_ligate_two_fragments():
    # Ligating two fragments
    fragments = ["ATGC", "CGTA"]
    ligate_step = Ligate(fragments, "Ligate_Product")
    product, _, _ = perform_ligate(ligate_step)
    expected_sequence = "ATGCCGTA"
    assert product.sequence == expected_sequence, "Ligate sequence mismatch for two fragments"

def test_ligate_three_fragments():
    # Ligating three fragments
    fragments = ["ATGC", "GATC", "TACG"]
    ligate_step = Ligate(fragments, "Ligate_Product")
    product, _, _ = perform_ligate(ligate_step)
    expected_sequence = "ATGCGATCTACG"
    assert product.sequence == expected_sequence, "Ligate sequence mismatch for three fragments"


### Gibson Tests
def test_gibson_two_fragments():
    # Overlap join two fragments
    fragments = ["ATGC", "CGTA"]
    gibson_step = Gibson(fragments, "Gibson_Product")
    product, _, _ = perform_gibson(gibson_step)
    expected_sequence = "ATGCCGTA"
    assert product.sequence == expected_sequence, "Gibson sequence mismatch for two fragments"

def test_gibson_three_fragments():
    # Overlap join three fragments
    fragments = ["ATGC", "GCTA", "TAGC"]
    gibson_step = Gibson(fragments, "Gibson_Product")
    product, _, _ = perform_gibson(gibson_step)
    expected_sequence = "ATGCGCTATAGC"
    assert product.sequence == expected_sequence, "Gibson sequence mismatch for three fragments"


### Golden Gate Tests
def test_goldengate_with_type_iis():
    # Using a Type IIS enzyme (BsaI) with one recognition site
    fragments = ["ATGCGTACGTAGCT", "TACGCTAGGCTAGC"]
    enzyme = "BsaI"
    gg_step = GoldenGate(fragments, enzyme, "GoldenGate_Product")
    product, _, _ = perform_goldengate(gg_step)
    expected_sequence = "ATGCGTACGTAGCTTACGCTAGGCTAGC"  # Simplified, assumes ligation without sticky ends
    assert product.sequence == expected_sequence, "Golden Gate sequence mismatch for Type IIS enzyme"

def test_goldengate_multiple_sites():
    # Using a Type IIS enzyme (BsaI) with multiple recognition sites
    fragments = ["ATGCGAAGTGGTAG", "TACGTTTGTGGTAG"]
    enzyme = "BsaI"
    gg_step = GoldenGate(fragments, enzyme, "GoldenGate_Product")
    product, _, _ = perform_goldengate(gg_step)
    expected_sequence = "ATGCGAAGTGGTAGTACGTTTGTGGTAG"
    assert product.sequence == expected_sequence, "Golden Gate sequence mismatch for multiple sites"
