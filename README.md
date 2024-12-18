
# BioE 134 Final Project Submission

## Project Overview

This project provides two core bioinformatics utilities: 

1. **Reverse Complement (revcomp)**: Calculates the reverse complement of a DNA sequence.
2. **Translate**: Translates a DNA sequence into a protein sequence according to the standard genetic code.

These functions are implemented in **Python** and are part of the broader bioinformatics toolset aimed at automating genetic sequence analysis tasks.

---

## Scope of Work

As part of the final project for BioE 134, I developed functions to validate construction files, as well as give an estimated cost of the reagents used as well as the time it will take to run the proccesses.

1. **perform_digest**: This function takes in a Step and performs the digest operation on it. It takes DNA and enzyme from the step and it uses the Bio library in python to perform the cut. It returns a Polynucleotide object as well as an estimated cost and time. 

2. **perform_gibson**: This function takes in a Step and performs gibson assembly by joining all the fragments given in the step. It returns a Polynucleotide object as well as an estimated cost and time of operation.

3. **perform_goldengate**: This function takes in a Step and performs golden gate on it. This also uses the Bio library to preform the cuts with enzymes. It returns a Polynucleotide object as well as an estimated cost and time of operation.
   
4. **perform_transform**: This function takes in a Step and simulates a transform operation. It outputs a string detailing the conditions, such as temperature, of the transformation as well as an estimated price and time.

5. **perform_pcr**: This function takes in a Step and simulates a PCR operation. It appends the forward, template, and reverse sequences and outputs it as a Polynucleotide. The time and cost are also outputed.

6. **perform_ligation**: This function takes in a Step and simulates a ligation operation. It appends the fragements together and outputs it as a Polynucleotide. The time and cost are also outputed.

7. **simulate**: This function runs all the operations in the input Construction File and outputs all of the products with their genetic sequence in a dictionary where the name is the key and the sequence is the value. The total estimated time and cost for all the operations is also included.


---

## Function Descriptions

### 1. Perform Gibson (`perform_gibson`)

- **Description**: This function simulates a Gibson assembly operation, combining multiple DNA fragments into a single sequence. The function calculates the cost and time of the assembly process and generates the final DNA product.  
- **Input**: A list of DNA fragment sequences or gene symbols to combine, and the output product name.  
- **Output**: A `Polynucleotide` object representing the assembled DNA sequence, along with cost and time predictions.  

**Example**:  
```python
perform_gibson(gibson_step)
# Returns: (Polynucleotide(sequence="ATGCGAATTCGCG"), 40, 45)
```

### 2. Perform PCR (`perform_pcr`)

- **Description**: This function simulates a PCR operation, amplifying a DNA template using specified forward and reverse primers. The function calculates the cost and time required to perform the PCR.
- **Input**: A Step object with forward primer, reverse primer, template sequence, and the output product name.
- **Output**: A `Polynucleotide` object representing the assembled DNA sequence, along with cost and time predictions.  

**Example**:  
```python
perform_pcr(pcr_step)
# Returns: (Polynucleotide(sequence="ATGCGAATTCGCG"), 40, 60)
```

### 3. Perform Golden Gate (`perform_goldengate`)

- **Description**: This function performs a Golden Gate assembly operation by digesting DNA fragments using a specific restriction enzyme and then ligating them into a single product. The function calculates cost and time predictions for the process.
- **Input**: A step with a list of DNA fragments or gene symbols, the restriction enzyme to use, and the output product name.
- **Output**: A `Polynucleotide` object representing the assembled DNA sequence, along with cost and time predictions.  

**Example**:  
```python
perform_goldengate(gg_step)
# Returns: (Polynucleotide(sequence="ATGCGAATTCGCG"), 30, 45)
```

### 4. Perform Ligate (`perform_ligate`)

- **Description**: This function simulates a ligation operation, combining multiple DNA fragments into a single linear sequence using DNA ligase. The function predicts the cost and time required for the ligation process.
- **Input**: A Step containing a list of DNA fragments or gene symbols to ligate and the output product name.
- **Output**: A Polynucleotide object representing the ligated DNA sequence, along with cost and time predictions.

**Example**:  
```python
perform_ligate(ligate_step)
# Returns: (Polynucleotide(sequence="ATGCGAATTCGCG"), 20, 45)
```

### 5. Perform Digest (`perform_digest`)

- **Description**: This function simulates a DNA digestion operation using one or more restriction enzymes. It calculates the cost and time required for the digestion and selects a specific fragment based on the given index.
- **Input**: A Step containing DNA sequence or gene symbol, a list of restriction enzymes, the fragment selection index, and the output product name.
- **Output**: A Polynucleotide object representing the digested DNA fragment, along with cost and time predictions.

**Example**:  
```python
perform_digest(digest_step)
# Returns: (Polynucleotide(sequence="GAATTC"), 20, 30)
```
### 6. Perform Transform (`perform_transform`)

- **Description**: This function simulates a transformation operation, introducing a plasmid into a host organism. It calculates the cost and time required for the transformation.
- **Input**: A Step containing a DNA sequence or gene symbol, host strain, antibiotics, temperature, and the output product name.
- **Output**: A success message describing the transformation, along with cost and time predictions.

**Example**:  
```python
perform_transform(transform_step)
# Returns: ("Transformation successful: Host carries plasmid.", 55, 45)
```
### 7. Simulate_CF (`simulate`)

- **Description**: This function simulates a series of molecular biology operations, such as PCR, Gibson assembly, Golden Gate assembly, digestion, ligation, and transformation. It processes each step sequentially, calculating the cumulative cost and time for the entire workflow while keeping track of the intermediate and final products.  
- **Input**: A `ConstructionFile` object containing a list of sequential steps with their respective operation types, inputs, and parameters.  
- **Output**: A dictionary mapping gene symbols to sequences, the total time required for the workflow, and the total cost of the workflow.  

**Example**:  
```python
simulate(constructionFile)
# Returns: ({"gene1": "ATGCGAATTCGCG", "gene2": "GAATTCGCGTAC"}, 300, 200)
```

## Error Handling

### 1. PCR (`perform_pcr`)
- Raises `ValueError` if any of the input sequences (forward primer, reverse primer, or template) are invalid or not found in the `symbol_to_gene` dictionary.  

### 2. Digest (`perform_digest`)
- Raises `ValueError` if any enzyme in the input list is not recognized or not found in `AllEnzymes`.  
- Raises `ValueError` if the fragment selection index is out of range.  

### 3. Gibson Assembly (`perform_gibson`)
- Raises `ValueError` if any of the DNA fragments are invalid or not found in the `symbol_to_gene` dictionary.  

### 4. Golden Gate Assembly (`perform_goldengate`)
- Raises `ValueError` if the specified restriction enzyme is not recognized.  
- Raises `ValueError` if the digestion step fails or results in no fragments.  

### 5. Ligation (`perform_ligate`)
- Raises `ValueError` if any of the DNA fragments are invalid or not found in the `symbol_to_gene` dictionary.  

### 6. Transformation (`perform_transform`)
- Raises `ValueError` if the input plasmid or antibiotics are invalid or not found in the `symbol_to_gene` or `reagent_to_price` dictionaries.  

### 7. Simulate_CF (`simulate-cf`)
- Raises `ValueError` if any step in the `ConstructionFile` contains invalid parameters or an unrecognized operation type.  
- Propagates errors raised by the individual functions (`perform_pcr`, `perform_digest`, etc.) when those operations fail.  


---

## Testing

Both functions have been tested with standard, edge, and invalid input cases. A comprehensive suite of tests has been implemented using **pytest**.

- **Test File**: `tests/test_operations.py` and `tests/test_basic_operations.py`

The tests include:
- example sequences both short and long
- length validations
- tests for functionality of each operation

## Conclusion

These functions provide operations for working construction files. 

