# Bioinformatics Master 1 project - UE HAU805I Stage M1 - Align program

## Introduction

Project done for the 1st year of Master degree in Bioinformatics' intership - Program in C++ able to calculate evolutionary distances between amino acid sequences from an aligned FASTA file and create a distance matrice using 5 methods:
- Distance estimation [Default].
- Jukes-Cantor model for amino acid [1].
- Poisson model for amino acid [2].
- Estimation of Kimura for PAM model [3].
- Estimation models for evolutionary distances between amino acid sequences: Poisson Correction and Equal-Input from Thomas Bigot and al., article [4].

## Installation

Clone the repository then compile the Align program with g++ or another C++ compiler.

```
git clone https://github.com/noeliepalermo/Align
g++ main.cpp -o align
```

## Usage

```
./align [evolutionary distances method option] aligned FASTA file [output file option]
```
The FASTA file must contain aligned proteins sequences.

You can use only one option for the evolutionary distances method and the output file.

### Evolutionary distances methods options:

```-d```,```--divergence```: Distance estimation [Default].

```-h```,```--help```: Display informations about the program.

```-jc```,```--jukescantor```: Jukes-Cantor model for amino acid.

```-k```, ```--kimura```: Estimation of Kimura for PAM model.

```-p```,```--poisson```: Poisson model for amino acid.

```-pc```,```--poissoncorrection```: Poisson Correction method from Thomas Bigot and al., article.

```-ei```,```--equalinput```: Equal-Input method from Thomas Bigot and al., article.

The two methods (Poisson Correction and Equal-Input) from Thomas Bigot and al., article, estimate evolutionary distances for 27 amino acid substitution models:
```AB```, ```BLOSUM62```, ```cpREV64```, ```cpREV```, ```Dayhoff``` [Default], ```DCMut-Dayhoff```, ```DCMut-JTT```, ```DEN```, ```FLU```, ```gcpREV```, ```HIVb```, ```HIVw```, ```JTT```, ```LG```, ```mtART```, ```mtInv```, ```mtMAM```, ```mtMet```, ```mtREV```, ```mtVer```, ```mtZOA```, ```PMB```, ```rtREV```, ```stmtREV```, ```VT```, ```WAG``` and ```WAG*```.

### Output file options:

```-m, --matrice```: Output file ```mat.dist```, with a triangular distance matrice in PHYLIP format.

```-o, --output```: Output file ```seqs.dist```, with a distance matrice and other informations: number of sequences, evolutionary distances and header sequences.

## Quick Demo

For testing the program Align, you can use the ```test_align.fasta``` file, which contains 26 proteins sequences from the PhylomeDB. Ignore gaps between all columns of the alignment for generate the expected results. Command to execute the test file:
```
./align -d test_align.fasta -m
```
The expected results are in the ```test_result_align.dist```.

## References
[1] Swofford D.L., Olsen G.J., W. P. and D.M., H. (1996). Phylogenetic inference. in : Hillis, d.m., moritz, c. and mable, b.k., eds. Molecular Systematics, 2nd Edition, Sinauer Associates, Sunderland (MA), (2) :407–514.

[2] Zuckerkandl, E., & Pauling, L. (1965). Evolutionary divergence and convergence in proteins. In Evolving genes and proteins (pp. 97-166). Academic Press.

[3] Kimura, M. (1991). The neutral theory of molecular evolution: a review of recent evidence. The Japanese Journal of Genetics, 66(4), 367-386.

[4] Bigot, T., Guglielmini, J., & Criscuolo, A. (2019). Simulation data for the estimation of numerical constants for approximating pairwise evolutionary distances between amino acid sequences. Data in brief, 25, 104212.
