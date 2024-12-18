import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pytest
from cf_simulator.polynucleotide import Polynucleotide, dsDNA
from cf_simulator.construction_file import PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from cf_simulator.construction_file_simulator import perform_pcr, perform_digest, perform_ligate, perform_gibson, perform_goldengate, perform_transform

# Mock the symbol_to_gene dictionary
symbol_to_gene = {
    "lacZ": "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC",
    "gfp": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCAGATGAACTTCAGGGTC",
}

# Mock reagent_to_price dictionary
reagent_to_price = {
    'Taq Polymerase': 25,
    'dNTPs': 15,
    'DNA Ligase': 20,
    'Gibson Assembly Mix': 40,
    'Competent Cells': 50,
    'EcoRI': 10,
    'BamHI': 10,
}

### PCR Tests
def test_pcr_basic():
    pcr_step = PCR("ATGC", "CGTA", "GATC", "PCR_Product")
    product, cost, time = perform_pcr(pcr_step)
    assert product.sequence == "ATGCGATCCGTA", "PCR sequence mismatch"
    assert cost == 40, "PCR cost incorrect"
    assert time == 45, "PCR time incorrect"

def test_pcr_long_template():
    pcr_step = PCR("ATGC", "CGTA", "A" * 1000, "PCR_Product")
    product, cost, time = perform_pcr(pcr_step)
    assert len(product.sequence) == 1008, "PCR sequence length mismatch"

### Digest Tests
def test_digest_basic():
    digest_step = Digest("ATGCGAATTCGATCG", ["EcoRI"], 0, "Digest_Product")
    product, cost, time = perform_digest(digest_step)
    assert "GAATTC" not in product.sequence, "EcoRI site not removed in digest"

def test_digest_multiple_enzymes():
    digest_step = Digest("ATGCGAATTCGATCG", ["EcoRI", "BamHI"], 0, "Digest_Product")
    product, cost, time = perform_digest(digest_step)
    assert "GAATTC" not in product.sequence, "EcoRI site not removed in digest"

### Ligate Tests
def test_ligate_basic():
    ligate_step = Ligate(["ATGC", "CGTA"], "Ligate_Product")
    product, cost, time = perform_ligate(ligate_step)
    assert product.sequence == "ATGCCGTA", "Ligate sequence mismatch"

def test_ligate_long_fragments():
    ligate_step = Ligate(["A" * 1000, "T" * 1000], "Ligate_Product")
    product, cost, time = perform_ligate(ligate_step)
    assert len(product.sequence) == 2000, "Ligate sequence length mismatch"


### Gibson Tests
def test_gibson_basic():
    gibson_step = Gibson(["ATGC", "CGTA"], "Gibson_Product")
    product, cost, time = perform_gibson(gibson_step)
    assert product.sequence == "ATGCCGTA", "Gibson sequence mismatch"

def test_gibson_long_fragments():
    gibson_step = Gibson(["A" * 1000, "T" * 1000], "Gibson_Product")
    product, cost, time = perform_gibson(gibson_step)
    assert len(product.sequence) == 2000, "Gibson sequence length mismatch"

### Golden Gate Tests
def test_goldengate_basic():
    gg_step = GoldenGate(["ATGC", "CGTA"], "BsaI", "GoldenGate_Product")
    product, cost, time = perform_goldengate(gg_step)
    assert product.sequence == "ATGCCGTA", "Golden Gate sequence mismatch"


def test_goldengate_with_site():
    gg_step = GoldenGate(["GAATTCATGC", "CGTAGAATTC"], "EcoRI", "GoldenGate_Product")
    product, cost, time = perform_goldengate(gg_step)
    assert "GAATTC" not in product.sequence, "Golden Gate site not removed"
