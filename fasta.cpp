/*
    Project in C++ done for the 1st year of Master degree in Bioinformatics' intership - University of Montpellier, France (2022-2024)
    Program Align in C++ able to calculate evolutionary distances between amino acids sequences from an aligned FASTA file and create a distance matrice using 5 methods.

    Class Fasta: Read FASTA file, checking amino acids sequences alignement and get informations (headers and sequences).

    Author: Noëlie PALERMO

    Contact: palermo.n@live.fr

    Version: "1.0"

    Date: 09/06/2023

    Licence: "This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <algorithm>

#include "fasta.hpp"

using namespace std;

// Function to print the Help manual
void Fasta::usage(int argc, char **argv)
{
    cout << "Usage: "<< argv[0] << " [evolutionary distances method option] aligned FASTA file [output file option] \n"
        << "Program in C++ able to calculate evolutionary distances between amino acids sequences from an aligned FASTA file and create a distance matrice.\n"
        << "Evolutionary distances methods options:\n"
        << "-d, --divergence         Distance estimation [Default].\n"
        << "-h, --help               Display informations about the program.\n"
        << "-jc, --jukescantor       Jukes-Cantor model for amino acids.\n"
        << "-k, --kimura             Kimura estimation for PAM model.\n"
        << "-p, --poisson            Poisson model for amino acids.\n"
        << "-pc, --poissoncorrection Poisson Correction method from Thomas Bigot and al., (2019) article.\n"
        << "-ei, --equalinput        Equal-Input method from Thomas Bigot and al., (2019) article.\n"
        << "For the Poisson Correction and Equal-Input methods, you can choose between one of the following 27 amino acids substitution models: \n"
        << "AB\n"
        << "BLOSUM62\n"
        << "cpREV64\n"
        << "cpREV\n"
        << "Dayhoff [Default]\n"
        << "DCMut-Dayhoff\n"
        << "DCMut-JTT\n"
        << "DEN\n"
        << "FLU\n"
        << "gcpREV\n"
        << "HIVb\n"
        << "HIVw\n"
        << "JTT\n"
        << "LG\n"
        << "mtART\n"
        << "mtInv\n"
        << "mtMAM\n"
        << "mtMet\n"
        << "mtREV\n"
        << "mtVer\n"
        << "mtZOA\n"
        << "PMB\n"
        << "rtREV\n"
        << "stmtREV\n"
        << "VT\n"
        << "WAG\n"
        << "WAG*\n"
        << "\n"
        << "Output file options:\n"
        << "-m, --matrice            Output file mat.dist, with a triangular distance matrice in PHYLIP format.\n"
        << "-o, --output             Output file seqs.dist, with a distance matrice and other informations: number of sequences, evolutionary distances and header sequences.\n"
        << endl;
}

//  Function to check if the FASTA file exist
bool Fasta::existe(int argc, char **argv)
{
    ifstream fichier(argv[2]);
    // If FASTA file doesn't exist: exit program
    if(! fichier.is_open()){
        cerr << "Error: the FASTA file can't be open.\n";
        return false;
        exit(-1);
    }
    fichier.close();
    return true; //  If sequences are aligned return true
}

// Function to stock FASTA file informations in struct fasta vector
vector<fasta> Fasta::vecteurFasta(int argc, char **argv)
{
    // FASTA file variable
    ifstream fichierFasta(argv[2]);

    // Variables for file's lines, headers and sequences
    string ligne, entete, sequence;

    // Variable to number headers
    int numerotation = 0;
    
    // Variable fasta to stock each informations in the struct fasta vector
    fasta stockage;

    // Variable to stock FASTA file's informations in a vector
    vector<fasta> vecFasta;

    while (getline(fichierFasta, ligne)) 
    {
        // If a line is empty, begin with a ">" and the header is not empty, the sequence is stocked.
        if( ligne.empty() || ligne[0] == '>' ){
            if( !entete.empty() ){
                stockage.s = sequence;
                vecFasta.push_back(stockage);
                entete.clear();
            }
            // If the line is not empty the header is stock
            if( !ligne.empty() ){
                entete = ligne.substr(1);
                stockage.e = entete;
                // Sequence number increase by 1
                stockage.n = numerotation++;
            }
            // Clear sequence variable
            sequence.clear();
        // Else if header is not empty and sequences continye on multiple lines, they are concatenated
       } else if( !entete.empty() ){
            if( ligne.find(' ') != std::string::npos ){ 
                entete.clear();
                sequence.clear();
            } else {
                sequence += ligne;
            }
        }
    }
    // The last sequence of the file is stock
    if( !entete.empty() ){
        stockage.s = sequence;
        vecFasta.push_back(stockage);
    }

    return vecFasta; // Return struct fasta vector
}

// Function to stock struct fasta length vector
int Fasta::tailleFasta(vector<fasta> vecFasta)
{
    int tailleVecteur = vecFasta.size();
    return tailleVecteur; // Return struct fasta length vector
}

// Function to check if the number of sequences in the FASTA file is strictly superior to 3
bool Fasta::superieurAtrois(vector<fasta> vecFasta)
{
    // Variable to count the number of sequences in struct fasta
    int compteur;

    for(int i = 0; i < vecFasta.size(); i++)
    {
        compteur++;
    }
    if (compteur < 3){
        cerr << "Error: The number of sequences must be equal or superior to 3 to create an evolutinary distance matrice.\n";
        return false;
        exit(-1);
    }
    
    return true;
}

// Function to check if sequences are aligned
bool Fasta::tailleSequence(vector<fasta> vecFasta, int tailleVecteur)
{
    for (int i = 0; i < tailleVecteur; i++)
    {
        for (int j = 1; j < tailleVecteur; j++)
        {
            if (vecFasta[i].s.size() != vecFasta[j].s.size())
            {
                // If sequences haven't the same length, they are not aligned: exit program
                cerr << "Error: Amino acids sequences are not aligned.\n";
                return false;
                exit(-1);
            }
            else
            {
                continue;
            }
        }
    }
    return true; // If sequences are aligned return true
}

// Function to create a new alignment with no gaps in columns
vector<fasta> Fasta::ignoreAllGaps(vector<fasta> vecFasta, int tailleVecteur)
{
    // Iterate gaps loop
    int g = 0;

    // Gaps positions vector
    vector<int> gaps;

    // New vector and variable for sequences alignement
    string newSequence;
    vector<string> newAlignement;

    // Get gaps positions into alignment columns
    while (g <= tailleVecteur-1){
        string sequenceComparee = vecFasta[g].s; // Iteration through compared sequence
        int tailleSequence = sequenceComparee.length(); // Variable length of comapared sequence
        for (int k = 0; k <= tailleSequence; k++) 
        {
            // If it's a gap site, it's alignement position is stocked in gaps' vector
            if (sequenceComparee[k] == '-') 
            {
                // Check if gap position is not already in the vector
                if (find(gaps.begin(), gaps.end(), k) == gaps.end()){
                    gaps.push_back(k);
                }
            }
        }
        g++;
    }

    // New alignment's loop iterator 
    int a = 0;

    // Creation of the new alignment without gaps
    while (a <= tailleVecteur-1){
        string sequenceComparee = vecFasta[a].s;  // Iteration through compared sequence
        int tailleSequence = sequenceComparee.length(); // Variable length of comapared sequence
        for (int k = 0; k < tailleSequence; k++) 
        {
            // If site position is already in gaps vector, it's ignored
            if (find(gaps.begin(), gaps.end(), k) != gaps.end()){
                continue;
            }else{
                // Else add the site to the new sequence
                newSequence += sequenceComparee[k];
            } 
        }
        // Add new sequence in new alignment vector
        newAlignement.push_back(newSequence);
        newSequence.clear();
        a++;
    }


    // Remplace sequences without gaps in struct fasta
    for (int i = 0; i < newAlignement.size(); i++)
    {
        vecFasta[i].s = newAlignement[i];
    }

    return vecFasta; // Return new alignment with no gaps in columns
}

// Function to create 1D dynamic which contains elements of struct fasta
fasta* Fasta::tableauFasta(vector<fasta> vecFasta, int tailleVecteur)
{
    fasta *tableauStruct = NULL;
    tableauStruct = new fasta[tailleVecteur];

    // If no memory can be allocated, program exit
    if (tableauStruct == NULL)
    {
        cerr << "Error: impossible to allocate memory.\n";
        exit(0); 
    }

    // Creation of the 1D dynamic with number, header and sequences
    for (int i = 0; i < tailleVecteur; i++)
    {
        tableauStruct[i].n = vecFasta[i].n;
        tableauStruct[i].e = vecFasta[i].e;
        tableauStruct[i].s = vecFasta[i].s;  
    }
    // Return 1D dynamic with elements of struct fasta
    return tableauStruct;
}