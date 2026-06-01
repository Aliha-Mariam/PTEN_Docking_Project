# PTEN Gene Analysis, Structure Validation and Molecular Docking of SF1670

## Introduction

PTEN (Phosphatase and Tensin Homolog) is one of the most important tumor suppressor genes in humans. Located on chromosome 10q23.31, PTEN encodes a dual-specificity phosphatase protein that negatively regulates the PI3K/AKT signaling pathway, a key pathway involved in cell growth, proliferation, survival, and metabolism. Loss of PTEN function through mutations, deletions, or reduced expression has been associated with numerous cancers, including breast cancer, prostate cancer, glioblastoma, and endometrial cancer.

Understanding the sequence, structure, and molecular interactions of PTEN is essential for studying its biological role and its involvement in disease progression. Advances in bioinformatics provide powerful computational methods to investigate genes and proteins without the need for extensive laboratory experiments. These approaches enable researchers to retrieve biological sequences, analyze evolutionary conservation, validate protein structures, and predict molecular interactions with potential therapeutic compounds.

In this project, a comprehensive bioinformatics pipeline was developed to study the PTEN gene from sequence to structure and molecular docking analysis. This project demonstrates the integration of sequence analysis, structural bioinformatics, protein validation, and molecular docking into a single computational workflow, providing insights into PTEN biology and its interaction with the inhibitor SF1670.

---

## Project Objectives

- Retrieve PTEN DNA sequence from NCBI
- Extract the Coding DNA Sequence (CDS)
- Transcribe DNA into mRNA
- Translate mRNA into protein
- Validate nucleotide and protein sequences using BLAST
- Retrieve experimental PTEN structure from PDB
- Retrieve AlphaFold predicted PTEN structure
- Compare experimental and AlphaFold structures
- Validate AlphaFold prediction using pLDDT scores
- Retrieve PTEN inhibitor SF1670 from PubChem
- Prepare ligand for docking
- Perform molecular docking using GNINA
- Analyze docking scores
- Visualize the best PTEN–SF1670 docking pose

---

## Software and Tools Used

| Tool | Description | Purpose in This Project |
|--------|-------------|-------------------------|
| Biopython | An open-source Python library for computational molecular biology and bioinformatics. It provides modules for sequence analysis, database access, BLAST searches, and structural biology. | Used for retrieving PTEN sequences from NCBI, reading FASTA and GenBank files, transcription, translation, BLAST analysis, and protein structure manipulation. |
| NCBI Entrez | A programmatic interface provided by the National Center for Biotechnology Information (NCBI) for accessing biological databases. | Used to retrieve PTEN nucleotide and protein sequences directly from NCBI databases. |
| BLAST (Basic Local Alignment Search Tool) | A sequence similarity search algorithm used to compare nucleotide or protein sequences against large biological databases. | Used to validate PTEN mRNA and protein sequences by identifying similar sequences and confirming sequence identity. |
| RCSB Protein Data Bank (PDB) | The world's largest repository of experimentally determined three-dimensional protein structures obtained through X-ray crystallography, NMR spectroscopy, and cryo-EM. | Used to download the experimentally determined PTEN structure (PDB ID: 1D5R) for structural analysis and molecular docking. |
| AlphaFold Protein Structure Database | A database of highly accurate protein structure predictions generated using DeepMind's AlphaFold artificial intelligence system. | Used to obtain the predicted PTEN structure for comparison with the experimental structure and structural validation. |
| Open Babel | A chemical informatics toolkit that converts molecular structures between different file formats and generates 3D molecular coordinates. | Used to convert the SF1670 ligand from SDF format into PDB and PDBQT formats required for docking simulations. |
| GNINA | A deep learning-enhanced molecular docking software built upon AutoDock Vina. It combines traditional docking algorithms with convolutional neural network (CNN) scoring functions. | Used to predict the binding orientation, binding affinity, and interaction of SF1670 with the PTEN protein active site. |
| Py3Dmol | A Python wrapper for 3Dmol.js that enables interactive visualization of molecular structures directly within Jupyter notebooks and Google Colab. | Used to visualize PTEN structures, AlphaFold models, structural superimposition results, and PTEN–SF1670 docking complexes. |
| NumPy | A scientific computing library for Python that provides efficient numerical operations and multidimensional array handling. | Used for coordinate calculations, active-site center determination, and structural data processing. |
| Pandas | A data analysis and manipulation library for Python that provides DataFrame structures for organizing tabular data. | Used to create BLAST result tables, docking score tables, and organize computational results. |
| Matplotlib | A data visualization library for Python used to generate scientific plots and graphs. | Used to generate pLDDT confidence plots and docking affinity score visualizations. |
| PubChem | A public chemical database maintained by the National Institutes of Health (NIH) that provides information about chemical compounds, structures, and biological activities. | Used to retrieve the PTEN inhibitor SF1670 (PubChem CID: 6442177) for ligand preparation and docking studies. |
---

## Dataset Information

### Gene

**PTEN (Phosphatase and Tensin Homolog)**

### NCBI RefSeq mRNA
NM_000314.8

### Protein Accession
NP_000305.3

### Experimental Structure
PDB ID: 1D5R

### AlphaFold Structure
AF-P60484-F1-model_v6

### Ligand
SF1670

### PubChem CID
6442177

---

# Pipeline Description

## 1. PTEN DNA Retrieval

```python
# NCBI sequence retrieval
handle = Entrez.efetch(
    db="nucleotide",
    id=ncbi_mrna_id,
    rettype="gb",
    retmode="text"
)

record = SeqIO.read(handle, "genbank")

handle.close()

print("Downloaded:", record.id)
print("Length:", len(record.seq))
```

This code retrieves the PTEN nucleotide sequence from the NCBI database using Biopython's Entrez module. The `Entrez.efetch()` function connects to the NCBI nucleotide database and downloads the record corresponding to the accession number stored in `ncbi_mrna_id` (NM_000314.8). The sequence is requested in GenBank format (`rettype="gb"`), which contains both the nucleotide sequence and detailed biological annotations such as coding regions, gene information, and transcript features. The downloaded data are then parsed using `SeqIO.read()` and stored in the variable `record` as a GenBank object. After retrieval, the connection to the NCBI server is closed using `handle.close()`. Finally, the accession number and total sequence length are printed to verify that the PTEN sequence was successfully downloaded.

---

## 2. Coding Sequence Extraction

```python
dna = None

for feature in record.features:

    if feature.type == "CDS":

        dna = feature.extract(record.seq)

        break

if dna is None:
    raise ValueError("CDS not found")

print("DNA Length:", len(dna))
```


The GenBank record downloaded from NCBI contains several annotated features, including genes, untranslated regions (UTRs), and coding sequences (CDS). Since only the coding sequence is translated into a protein, the CDS must be extracted before proceeding to transcription and translation.

This code iterates through all annotated features in the PTEN GenBank record using a `for` loop. For each feature, it checks whether the feature type is `"CDS"`. Once the coding sequence is found, the `extract()` function retrieves the exact nucleotide sequence corresponding to the PTEN coding region and stores it in the variable `dna`. The loop then stops using the `break` statement because only one CDS is required.

An error-checking step is included using:

```python
if dna is None:
    raise ValueError("CDS not found")
```

This ensures that the program stops if no coding sequence is detected in the record.

---

## 3. DNA Transcription

```python
# Read DNA from FASTA file
record = SeqIO.read("PTEN_DNA.fasta", "fasta")
dna = record.seq

# Transcribe DNA → mRNA
mrna = dna.transcribe()

print("mRNA generated successfully")
print("Length of mRNA:", len(mrna))

# Preview
print("\nFirst 100 bases of mRNA:\n")
print(mrna[:100])
```

This code first reads the previously saved PTEN DNA sequence from the FASTA file using the `SeqIO.read()` function and stores it in the variable `dna`. The `transcribe()` method is then applied to convert the DNA sequence into its corresponding mRNA sequence. During transcription, thymine (T) nucleotides present in DNA are replaced by uracil (U) nucleotides in RNA, while the remaining bases remain unchanged according to standard transcription rules.

Once transcription is completed, the program prints a confirmation message, the total length of the generated mRNA sequence, and the first 100 nucleotides as a preview. This allows verification that the transcription process was performed successfully.

---

## 4. mRNA Sequence Validation Using BLASTN


```python
# Read mRNA from FASTA file
record = SeqIO.read("PTEN_mRNA.fasta", "fasta")
mrna_seq = record.seq

print("Running BLASTN on PTEN mRNA...")

# Run BLASTN
blast_result = NCBIWWW.qblast(
    program="blastn",
    database="nt",
    sequence=str(mrna_seq),
    hitlist_size=5
)

# Save BLAST output
blast_xml = blast_result.read()

with open("PTEN_mRNA_BLAST.xml", "w") as f:
    f.write(blast_xml)

print("BLASTN completed and saved.")
```

The code begins by reading the PTEN mRNA sequence from the FASTA file and storing it in the variable `mrna_seq`. The `NCBIWWW.qblast()` function is then used to submit the sequence to the online BLAST service. The parameter `program="blastn"` specifies nucleotide-to-nucleotide alignment, while `database="nt"` instructs BLAST to search against the comprehensive NCBI nucleotide database. The `hitlist_size=5` parameter limits the output to the top five matching sequences.

The BLAST results are returned in XML format, which contains detailed alignment information including sequence identities, alignment lengths, scores, and E-values. The XML output is saved as `PTEN_mRNA_BLAST.xml` for further analysis and interpretation.

---

## 5. Analysis of BLASTN Results

```python
# Open saved BLAST file
with open("PTEN_mRNA_BLAST.xml") as f:
    blast_record = NCBIXML.read(f)

rows = []

# Take top 5 hits
for i, hit in enumerate(blast_record.alignments[:5]):

    hsp = hit.hsps[0]

    identity = (hsp.identities / hsp.align_length) * 100

    rows.append([
        i + 1,
        hit.title,
        hsp.align_length,
        round(identity, 2),
        hsp.expect
    ])

blast_table = pd.DataFrame(
    rows,
    columns=["Rank", "Hit", "Alignment Length", "Identity %", "E-value"]
)

blast_table
```

After performing the BLASTN search, the resulting XML file contains detailed information about all sequence alignments returned by NCBI. This code parses the XML file and extracts the most relevant information from the top five BLAST hits.

The XML file is first loaded using the `NCBIXML.read()` function, which converts the BLAST output into a structured Biopython object. An empty list called `rows` is then created to store information about each alignment.

The code iterates through the first five alignments returned by BLAST. For each alignment, the first High-Scoring Pair (HSP) is selected because it represents the best alignment segment between the query sequence and the database sequence. The percentage sequence identity is calculated using:

```python
identity = (hsp.identities / hsp.align_length) * 100
```

Several important alignment statistics are extracted:

- Rank of the hit
- Description of the matching sequence
- Alignment length
- Percentage identity
- E-value

These values are stored in a list and later converted into a Pandas DataFrame, creating a clean and readable results table.


---

## 6. Translation of mRNA into Protein Sequence


```python
# Read mRNA sequence
record = SeqIO.read("PTEN_mRNA.fasta", "fasta")

mrna_seq = record.seq

# Translate mRNA into protein
protein = mrna_seq.translate(to_stop=True)

print("Protein Length:", len(protein))
print("\nFirst 100 amino acids:\n")
print(protein[:100])
```

Following transcription, the generated PTEN mRNA sequence was translated into its corresponding protein sequence. 
The code begins by reading the PTEN mRNA sequence from the FASTA file using the `SeqIO.read()` function. The sequence is stored in the variable `mrna_seq`.The `translate()` method is then applied to convert the mRNA sequence into a protein sequence. During this process, every three nucleotides (codon) are interpreted according to the standard genetic code and converted into the corresponding amino acid.
The parameter:

```python
to_stop=True
```

instructs Biopython to stop translation when the first stop codon is encountered. This ensures that only the biologically relevant protein sequence is generated and excludes any downstream untranslated regions.
After translation, the code displays:
- Total protein length
- First 100 amino acids of the translated protein sequence

---

## 7. Protein Sequence Validation Using BLASTP

After translating the PTEN mRNA sequence into a protein sequence, the resulting protein was validated in two ways.

First, the official PTEN protein sequence (Accession: **NP_000305.3**) was downloaded from the NCBI Protein database using Entrez. The translated protein generated from the mRNA sequence was then compared directly with the official PTEN protein sequence. Matching sequence lengths and identical amino acid sequences confirmed that the translation process was accurate and produced the correct PTEN protein.

To further verify the biological identity of the protein, a **BLASTP** search was performed against the NCBI non-redundant protein database (**nr**). BLASTP compares protein sequences and identifies similar proteins across different organisms. The top five hits were extracted and summarized in a table containing alignment length, percentage identity, and E-values.


---

## 8. Retrieval and Characterization of the Experimental PTEN Structure

```python
# Download PDB structure
pdb_id = "1D5R"

url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

response = requests.get(url)

if response.status_code == 200:
    with open("PTEN_1D5R.pdb", "w") as f:
        f.write(response.text)
    print("PDB downloaded successfully:", pdb_id)
else:
    print("Download failed")


file_path = "PTEN_1D5R.pdb"

with open(file_path) as f:
    lines = f.readlines()

info = {
    "PDB ID": None,
    "Title": None,
    "Method": None,
    "Resolution": None,
}

for line in lines:

    if line.startswith("HEADER"):
        info["PDB ID"] = line.split()[-1]

    if line.startswith("TITLE"):
        info["Title"] = line[10:].strip()

    if line.startswith("EXPDTA"):
        info["Method"] = line[10:].strip()

    if "RESOLUTION." in line:
        parts = line.split()
        info["Resolution"] = parts[-2] + " Å"

print(info)
```

To perform structural analysis and molecular docking, an experimentally determined three-dimensional structure of PTEN was required. The PTEN crystal structure with PDB ID **1D5R** was downloaded directly from the Protein Data Bank (PDB) using the Requests library.

The downloaded structure was saved as `PTEN_1D5R.pdb`, which contains atomic coordinates and structural annotations for the protein. After downloading, the PDB file was parsed line by line to extract important metadata describing the structure.

The code retrieves:

- **PDB ID** – Unique identifier of the structure.
- **Title** – Description of the deposited structure.
- **Experimental Method** – Technique used to determine the structure.
- **Resolution** – Quality measure of the experimentally determined structure.

These details were stored in a dictionary and displayed for inspection.


### Output

```python
{
 'PDB ID': '1D5R',
 'Title': 'CRYSTAL STRUCTURE OF HUMAN PTEN',
 'Method': 'X-RAY DIFFRACTION',
 'Resolution': '2.1 Å'
}
```

---

## 9. Identification of Heteroatoms and Receptor Preparation for Docking

```python
# Extract Ligands / HETATM (PTEN structure contents)

het_residues = []

with open("PTEN_1D5R.pdb", "r") as f:
    for line in f:
        if line.startswith("HETATM"):
            het_residues.append(line[17:20].strip())

import pandas as pd

het_table = pd.Series(het_residues).value_counts().reset_index()
het_table.columns = ["Residue", "Atom_count"]

het_table
```

### Output

| Residue | Atom_count |
|----------|------------|
| HOH | 382 |
| TLA | 10 |

---

```python
# Clean receptor for docking

pdb_file = "PTEN_1D5R.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("PTEN", pdb_file)

class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "

class WaterSelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() == "HOH"

class LigandSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] != " " and residue.get_resname() != "HOH"

io = PDBIO()
io.set_structure(structure)

io.save("PTEN_cleanprotein.pdb", ProteinSelect())
io.save("PTEN_waters.pdb", WaterSelect())
io.save("PTEN_ligands.pdb", LigandSelect())

print("Files saved successfully:")
print("PTEN_cleanprotein.pdb")
print("PTEN_waters.pdb")
print("PTEN_ligands.pdb")
```

Before molecular docking can be performed, the experimental protein structure must be carefully prepared. Protein Data Bank files often contain additional molecules besides the protein itself, including crystallographic water molecules, ions, cofactors, and bound ligands. These components are represented as **HETATM** records within the PDB file.

The first part of the analysis scanned the PTEN structure and extracted all HETATM residues. Two types of non-protein molecules were identified:
- **HOH** – Water molecules
- **TLA** – Tartaric acid molecules present in the crystal structure
The frequency of each residue was calculated and displayed as a summary table.
After identifying the non-protein components, the structure was separated into three independent files using Biopython's `PDBIO` and `Select` classes.
### Generated Files
| File | Contents |
|--------|----------|
| PTEN_cleanprotein.pdb | Protein atoms only |
| PTEN_waters.pdb | Water molecules (HOH) |
| PTEN_ligands.pdb | Non-water heteroatoms and ligands (TLA) |

By separating the structure into protein, water, and ligand components, a clean receptor model is generated for docking studies. The resulting file:

```text
PTEN_cleanprotein.pdb
```
contains only the PTEN protein and serves as the receptor structure for all subsequent docking analyses with the inhibitor SF1670.


---

## 10. Retrieval of the AlphaFold Predicted PTEN Structure

```python
# Download AlphaFold structure

url = "https://alphafold.ebi.ac.uk/files/AF-P60484-F1-model_v6.pdb"

r = requests.get(url)

if r.status_code == 200 and "ATOM" in r.text:
    with open("PTEN_alphafold.pdb", "w") as f:
        f.write(r.text)
    print("AlphaFold PTEN structure downloaded successfully")
else:
    print("Download failed — use manual download from AlphaFold entry page")
    print("https://alphafold.ebi.ac.uk/entry/P60484")
```

In addition to the experimentally determined PTEN structure obtained from the Protein Data Bank, a computationally predicted PTEN structure was retrieved from the AlphaFold Protein Structure Database.

AlphaFold is an artificial intelligence system  that predicts highly accurate three-dimensional protein structures directly from amino acid sequences. The AlphaFold Protein Structure Database provides predicted structures for millions of proteins across different organisms.
The code downloads the AlphaFold model corresponding to the human PTEN protein (**UniProt ID: P60484**) and saves it locally as:

```text
PTEN_alphafold.pdb
```

Before saving, the code verifies that:
- The download request was successful (`status_code == 200`)
- The file contains atomic coordinate records (`ATOM`)
This ensures that a valid protein structure file has been retrieved.

---

## 11.Visualization of the Experimental PTEN Structure

```python
display(HTML("<h3>Experimental PTEN Structure (Clean Protein)</h3>"))

with open("PTEN_cleanprotein.pdb", "r") as f:
    exp_pdb = f.read()

view = py3Dmol.view(width=850, height=450)
view.addModel(exp_pdb, "pdb")
view.setStyle({"cartoon": {"color": "spectrum"}})
view.zoomTo()
view.show()

display(HTML("<hr><br>"))
```

In this step, the cleaned experimental PTEN protein structure is visualized using the Py3Dmol library. The file `PTEN_cleanprotein.pdb`, which was generated after removing water molecules and non-protein components, is loaded and rendered in an interactive 3D format. The structure is displayed in cartoon representation, which is commonly used in structural biology to clearly show the overall folding pattern of proteins, including helices, sheets, and loop regions.The color scheme “spectrum” is applied so that different regions of the protein are visually distinguished, making it easier to interpret structural organization. The function `zoomTo()` automatically adjusts the view so that the entire protein fits within the display window. This step helps in visually confirming that the receptor protein is correctly prepared and structurally intact before proceeding to docking analysis. 

---

## 12. Visualization of AlphaFold Predicted PTEN Structure


```python
# Visualize AlphaFold PTEN

display(HTML("<h3>AlphaFold Predicted PTEN Structure</h3>"))

with open("PTEN_alphafold.pdb", "r") as f:
    af_pdb = f.read()

view = py3Dmol.view(width=850, height=450)

view.addModel(af_pdb, "pdb")
view.setStyle({"cartoon": {"color": "spectrum"}})

view.zoomTo()
view.show()

display(HTML("<hr><br>"))
```

In this step, the computationally predicted PTEN structure obtained from the AlphaFold database is visualized using Py3Dmol. The file `PTEN_alphafold.pdb` contains the three-dimensional coordinates of the PTEN protein predicted by AlphaFold based on its amino acid sequence.The structure is displayed in a cartoon representation, which helps in clearly visualizing the secondary structural elements such as alpha helices, beta sheets, and loop regions. The spectrum color scheme is applied to differentiate different regions of the protein and make the structure more interpretable.This visualization allows comparison of the predicted PTEN structure with the experimentally determined structure, helping to assess structural similarity and confidence in the prediction. The view is automatically centered and scaled using `zoomTo()` so that the entire protein is visible in a single frame.

---

## 13. Structural Alignment of Experimental and AlphaFold PTEN

```python
# LOAD STRUCTURES
exp_file = "PTEN_cleanprotein.pdb"
af_file = "PTEN_alphafold.pdb"

parser = PDBParser(QUIET=True)
exp_structure = parser.get_structure("EXP_PTEN", exp_file)
af_structure = parser.get_structure("AF_PTEN", af_file)

print("Structures loaded successfully")


# EXTRACT C-ALPHA ATOMS
exp_atoms = []
af_atoms = []

for exp_res in exp_structure.get_residues():
    exp_id = exp_res.get_id()

    for af_res in af_structure.get_residues():
        af_id = af_res.get_id()

        # match based on residue number + chain
        if exp_id[1] == af_id[1]:
            if "CA" in exp_res and "CA" in af_res:
                exp_atoms.append(exp_res["CA"])
                af_atoms.append(af_res["CA"])
            break

print("Matched CA atoms:", len(exp_atoms))


# ALIGN STRUCTURES
sup = Superimposer()

n = len(exp_atoms)

sup.set_atoms(exp_atoms, af_atoms)
sup.apply(af_structure.get_atoms())

print("Aligned residues:", n)
print("RMSD (Å):", round(sup.rms, 3))


# VISUAL OVERLAY
display(HTML("<h3>PTEN Structural Superimposition (Experimental vs AlphaFold)</h3>"))

with open("PTEN_cleanprotein.pdb") as f:
    exp_pdb = f.read()

with open("PTEN_alphafold.pdb") as f:
    af_pdb = f.read()

view = py3Dmol.view(width=900, height=600)

# Experimental (red)
view.addModel(exp_pdb, "pdb")
view.setStyle({"model": 0}, {"cartoon": {"color": "red"}})

# AlphaFold (blue)
view.addModel(af_pdb, "pdb")
view.setStyle({"model": 1}, {"cartoon": {"color": "blue"}})

view.zoomTo()
view.show()
```

### Output

- Matched C-alpha atoms: **307**
- Aligned residues: **307**
- RMSD (Å): **1.555**

In this step, both the experimentally determined PTEN structure and the AlphaFold predicted structure were aligned using C-alpha (Cα) atoms, which represent the backbone of the protein. A total of 307 residues were successfully matched and used for structural alignment, indicating strong overlap between the two structures.

The Superimposer algorithm calculated an RMSD (Root Mean Square Deviation) value of **1.555 Å**, which represents the average distance between corresponding atoms after optimal superposition. Since RMSD values below ~2 Å generally indicate good structural similarity, this result suggests that the AlphaFold predicted PTEN structure closely resembles the experimentally determined structure.

Finally, both structures were visually superimposed, with the experimental structure shown in red and the AlphaFold structure shown in blue. This allows clear visual confirmation of structural similarity and minor deviations between the two models.


---

## 14. AlphaFold pLDDT Confidence Analysis
```python
# pLDDT CONFIDENCE PLOT
parser = PDBParser(QUIET=True)
structure = parser.get_structure("AF_PTEN", "PTEN_alphafold.pdb")

plddt = []

for model in structure:
    for chain in model:
        for residue in chain:
            if "CA" in residue:
                plddt.append(residue["CA"].get_bfactor())

plt.figure(figsize=(12,4))
plt.plot(plddt, linewidth=1)

plt.axhline(90, linestyle="--", label="Very high confidence")
plt.axhline(70, linestyle="--", label="Confident")
plt.axhline(50, linestyle="--", label="Low confidence")

plt.title("PTEN AlphaFold pLDDT Confidence Score")
plt.xlabel("Residue Index")
plt.ylabel("pLDDT Score")
plt.legend()
plt.tight_layout()
plt.show()


# CALCULATE pLDDT STATISTICS
average_plddt = sum(plddt) / len(plddt)

high_confidence = sum(1 for x in plddt if x > 90)
good_confidence = sum(1 for x in plddt if x > 70)
low_confidence = sum(1 for x in plddt if x < 50)

print("Structure Validation Summary")
print("----------------------------")
print("Number of residues:", len(plddt))
print("Average pLDDT:", round(average_plddt, 2))
print("Residues > 90:", high_confidence)
print("Residues > 70:", good_confidence)
print("Residues < 50:", low_confidence)
```

### Results

- Number of residues: **403**
- Average pLDDT: **83.03**
- Residues > 90 (very high confidence): **277**
- Residues > 70 (confident): **313**
- Residues < 50 (low confidence): **69**
In this step, the confidence of the AlphaFold predicted PTEN structure is evaluated using the pLDDT (predicted Local Distance Difference Test) score. This score is stored in the B-factor field of the AlphaFold PDB file and provides a residue-wise measure of prediction reliability.
The code extracts pLDDT values for each amino acid by reading the C-alpha atoms from the structure. These values are then plotted to visualize how confidence varies across the protein sequence. Threshold lines are added at 90, 70, and 50 to classify regions into very high confidence, confident, and low-confidence regions respectively.
Statistical analysis shows that the protein has an average pLDDT score of **83.03**, which indicates generally good model reliability. A large number of residues (277) have very high confidence (>90), while only a smaller portion (69 residues) fall into the low-confidence region (<50), typically corresponding to flexible or disordered regions.
Overall, this analysis confirms that the AlphaFold PTEN structure is reliable for downstream applications such as structural comparison and molecular docking.

---

## 15. Ligand Selection and Preparation for Molecular Docking

```python
# Download ligand
cid = "6442177"  # SF1670
sdf_file = "SF1670.sdf"

urls = [
    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d",
    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
]

downloaded = False

for url in urls:
    print("Trying:", url)
    r = requests.get(url, timeout=30)

    if r.status_code == 200 and "M  END" in r.text:
        with open(sdf_file, "w") as f:
            f.write(r.text)

        print("Ligand downloaded successfully from:", url)
        downloaded = True
        break

if not downloaded:
    raise Exception("Ligand download failed")


# prepare ligand for docking
!obabel SF1670.sdf -O SF1670.pdb
!obabel SF1670.sdf -O SF1670.pdbqt
```

In this step, the ligand molecule SF1670 is downloaded from the PubChem database and prepared for molecular docking. The ligand is retrieved using its PubChem Compound ID (CID: 6442177), which ensures that the correct chemical structure is obtained from a reliable public chemical repository.
The code first tries two different PubChem API links: one requesting the 3D structure and another fallback option. The response is checked to ensure it contains valid molecular structure data (indicated by the presence of `"M  END"` in the SDF file format). Once confirmed, the ligand structure is saved locally as `SF1670.sdf`.

After downloading, Open Babel is used to convert the ligand into two different formats:
- **PDB format (`SF1670.pdb`)**: used for visualization and structural inspection  
- **PDBQT format (`SF1670.pdbqt`)**: required for docking simulations in GNINA
The PDBQT format includes additional chemical information such as charges and rotatable bonds, which are essential for flexible docking.

## Why SF1670 was selected as ligand?
SF1670 was chosen because it is a biologically relevant small-molecule inhibitor of PTEN. It is widely used in research studies related to PTEN regulation and cancer signaling pathways. Since PTEN is a tumor suppressor protein, SF1670 is useful for studying how inhibition affects PTEN activity.
This makes SF1670 an ideal ligand for this project because:
- It has known biological interaction with PTEN  
- It is available in public chemical databases (PubChem)  
- It provides a realistic case study for molecular docking  
- It helps validate the docking pipeline using a meaningful inhibitor
  
---

## 16. Identification of Docking (Active Site) Center

```python
# FIND DOCKING CENTER

pocket_residues = ["HIS", "CYS", "ARG", "LYS", "GLY"]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("PTEN", "PTEN_cleanprotein.pdb")

coords = []

for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() in pocket_residues:
                for atom in residue:
                    coords.append(atom.coord)

coords = np.array(coords)

center_x, center_y, center_z = coords.mean(axis=0)

print("NEW ACTIVE SITE CENTER:")
print(center_x, center_y, center_z)
```

### Output

```
NEW ACTIVE SITE CENTER:
35.06179 83.71493 30.46342
```
In this step, the active site (docking center) of the PTEN protein is identified. Molecular docking requires a defined 3D region (grid box center) where the ligand will be searched and positioned during simulation. Therefore, it is necessary to estimate the coordinates of the binding pocket.
The code selects amino acid residues that are commonly involved in protein–ligand interactions, including:
- Histidine (HIS)  
- Cysteine (CYS)  
- Arginine (ARG)  
- Lysine (LYS)  
- Glycine (GLY)
These residues are often found in or near functional binding sites due to their chemical properties such as charge, flexibility, and ability to form hydrogen bonds. After selecting these residues from the cleaned PTEN protein structure, the coordinates of all their atoms are collected. The mean (average) position of these coordinates is then calculated using NumPy. This average point represents the estimated center of the binding pocket.

## Final Docking Center

The computed active site center for docking is:
- **X = 35.06179**
- **Y = 83.71493**
- **Z = 30.46342**

This step is important because docking accuracy depends heavily on correctly defining the binding region. By selecting functionally important residues and calculating their spatial center, the ligand is guided toward the most biologically relevant region of PTEN.
This ensures that the docking simulation focuses on a realistic active site rather than random regions of the protein surface.

---

## 17. Molecular Docking using GNINA


```python
# RUN GNINA DOCKING
!./gnina \
-r PTEN_cleanprotein.pdb \
-l SF1670.pdbqt \
--center_x 35.06179 \
--center_y 83.71493 \
--center_z 30.46342 \
--size_x 28 \
--size_y 28 \
--size_z 30 \
--exhaustiveness 10 \
--num_modes 10 \
--seed 42 \
-o PTEN_SF1670_docked.sdf


# EXTRACT BEST POSE
!obabel PTEN_SF1670_docked.sdf -O PTEN_best_pose.pdb -f 1 -l 1


# SHOW RESULT SUMMARY
print("Docking completed successfully")
print("Best pose file created: PTEN_best_pose.pdb")
```

In this step, molecular docking was performed using GNINA, a deep learning–enhanced docking tool based on AutoDock Vina. The cleaned PTEN protein structure is used as the receptor, and SF1670 (prepared in PDBQT format) is used as the ligand.

The docking box is defined using the previously calculated active site center:
- X = 35.06179  
- Y = 83.71493  
- Z = 30.46342  

This ensures that the ligand searches for binding poses only within the biologically relevant region of the protein.
The parameters used in GNINA include:

- **Exhaustiveness = 10** → controls search intensity (higher = more accurate but slower)  
- **Num_modes = 10** → returns top 10 binding poses  
- **Seed = 42** → ensures reproducibility of results  
- **Box size (28 Å)** → defines search space around active site  

After docking, GNINA generates multiple ligand poses with predicted binding affinities. The pose with the lowest (most negative) binding energy is considered the best binding conformation.
Finally, Open Babel is used to extract the best docking pose and convert it into a PDB file for visualization.

### Results

| Mode | Affinity (kcal/mol) | Intramol (kcal/mol) | CNN Pose Score | CNN Affinity |
|------|----------------------|----------------------|----------------|---------------|
| 1    | -10.47               | 265.96               | 0.3087         | 4.364         |
| 2    | -8.12                | 70.40                | 0.3050         | 4.485         |
| 3    | -9.18                | 72.33                | 0.2862         | 4.434         |
| 4    | -9.94                | 73.55                | 0.2844         | 4.374         |
| 5    | -8.83                | 71.84                | 0.2735         | 4.256         |
| 6    | -8.61                | 72.97                | 0.2607         | 4.471         |
| 7    | -7.56                | 69.20                | 0.2549         | 4.594         |
| 8    | -8.36                | 73.60                | 0.2522         | 4.417         |
| 9    | -9.04                | 71.08                | 0.2511         | 4.394         |
| 10   | -7.98                | 72.48                | 0.2499         | 4.475         |

- Best binding affinity: **-10.47 kcal/mol**
- Best pose successfully extracted for structural analysis and visualization

---
## 18. Docking Score Table and Visualization

```python
# SCORE TABLE
scores = [-10.47, -8.12, -9.18, -9.94, -8.83,
          -8.61, -7.56, -8.36, -9.04, -7.98]

df = pd.DataFrame({
    "Pose": list(range(1, 11)),
    "Affinity (kcal/mol)": scores
})

df


# PLOT DOCKING RESULTS
plt.figure(figsize=(6,4))
plt.bar(df["Pose"], df["Affinity (kcal/mol)"])
plt.gca().invert_yaxis()
plt.xlabel("Pose")
plt.ylabel("Binding Affinity (kcal/mol)")
plt.title("PTEN–SF1670 Docking Results")
plt.show()
```

In this step, the docking results obtained from GNINA are organized into a structured table using the Pandas library. Each docking pose (from 1 to 10) is assigned its corresponding binding affinity value in kcal/mol.

The table helps in clearly comparing how different ligand conformations interact with the PTEN protein. Since binding affinity is a key indicator of interaction strength, lower (more negative) values represent stronger binding.

| Pose | Affinity (kcal/mol) |
|------|----------------------|
| 1    | -10.47               |
| 2    | -8.12                |
| 3    | -9.18                |
| 4    | -9.94                |
| 5    | -8.83                |
| 6    | -8.61                |
| 7    | -7.56                |
| 8    | -8.36                |
| 9    | -9.04                |
| 10   | -7.98                |

After creating the table, a bar plot is generated using Matplotlib to visually represent the docking results. The x-axis represents different docking poses, while the y-axis represents binding affinity values. The y-axis is inverted so that stronger (more negative) binding affinities appear higher in the graph, making interpretation more intuitive.


## Key Insight

The graph clearly shows that **Pose 1 has the strongest binding affinity (-10.47 kcal/mol)** among all generated poses, indicating it is the most favorable docking conformation of SF1670 with PTEN.

---

## 19: Visualization of Best Docked Complex

```python
# visualization of the best docking result
display(HTML("<h3>PTEN–SF1670 Docked Complex</h3>"))

with open("PTEN_cleanprotein.pdb") as f:
    protein = f.read()

with open("PTEN_best_pose.pdb") as f:
    ligand = f.read()

view = py3Dmol.view(width=850, height=500)

view.addModel(protein, "pdb")
view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})

view.addModel(ligand, "pdb")
view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})

view.zoomTo()
view.show()
```

In this step, the final docked complex of PTEN protein with the SF1670 ligand is visualized using Py3Dmol, a molecular visualization tool.
First, the cleaned PTEN protein structure is loaded from the PDB file, and the best docking pose of the ligand (SF1670) is also loaded from the GNINA output. These two structures represent the receptor (protein) and ligand (drug molecule) after docking.
The protein is displayed in a **cartoon representation** with a spectrum color scheme, which helps in clearly visualizing the overall protein structure and secondary structural elements such as helices and sheets.The ligand is displayed in a **stick representation** with green coloring, which highlights its position and orientation inside the protein binding pocket.Finally, both structures are shown together in a single interactive 3D view, allowing clear visualization of how SF1670 binds within the PTEN active site.



This visualization confirms the docking result by showing the physical binding location of SF1670 within the PTEN protein, supporting the predicted binding affinity and docking analysis.


---

## Conclusion

This project successfully implemented a complete bioinformatics and molecular docking pipeline for the PTEN protein and SF1670 ligand.
The PTEN gene was retrieved from NCBI, converted into mRNA and protein, and validated using BLAST, confirming sequence accuracy. Structural analysis using experimental (PDB: 1D5R) and AlphaFold models showed strong similarity with an RMSD of **1.555 Å** and good AlphaFold confidence (average pLDDT = **83.03**).Molecular docking using GNINA identified SF1670 as a strong binder to PTEN, with a best binding affinity of **-10.47 kcal/mol**. Visualization confirmed stable interaction within the active site.The results demonstrate a successful end-to-end pipeline from gene retrieval to protein–ligand docking, highlighting SF1670 as a potential PTEN-interacting compound and validating the use of computational tools in drug discovery.

---

## Workflow

```text
PTEN DNA Sequence Retrieval
            ↓
Coding Sequence Extraction
            ↓
mRNA Transcription
            ↓
Protein Translation
            ↓
BLAST Validation
            ↓
Experimental Structure Retrieval
            ↓
AlphaFold Structure Retrieval
            ↓
Structure Comparison & Validation
            ↓
Ligand Retrieval (SF1670)
            ↓
Ligand Preparation
            ↓
Molecular Docking (GNINA)
            ↓
Docking Analysis
            ↓
Protein-Ligand Visualization
```
---

## References

1. Biopython Documentation — https://biopython.org
2. NCBI Database — https://www.ncbi.nlm.nih.gov
3. Protein Data Bank — https://www.rcsb.org
4. AlphaFold Database — https://alphafold.ebi.ac.uk
5. PubChem Database — https://pubchem.ncbi.nlm.nih.gov
6. Open Babel — https://openbabel.org
7. GNINA — https://github.com/gnina/gnina
8. PTEN Structure (PDB ID: 1D5R)
9. SF1670 (PubChem CID: 6442177)

---
