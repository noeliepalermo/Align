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
#include <vector>
#include <string.h>

#ifndef FASTA_HPP
#define FASTA_HPP

/*
 Structure pour contenir les informations d'un fichier FASTA (entête et séquences protéiques alignées)
 Une numérotation des entêtes a été ajoutée pour faciliter l'identification des séquences si les entêtes ne sont pas au format du NCBI
*/
struct fasta
{
  std::string e; // Header
  int n; // Header number
  std::string s; // Sequence
};


class Fasta
{
  public:
    // Fasta Class constructor
    Fasta()
    {
      std::cout << "Fasta Class constructor.\n";
    };

    // Fasta Class destructor
    ~Fasta()
    {
      std::cout << "Fasta Class destructor.\n";
    };

    // Function to print the Help manual.
    void usage(int argc, char **argv);

    // Function to check if the FASTA file exist.
    bool existe(int argc, char **argv);
    
    // Function to stock FASTA file informations into struct fasta vector
    std::vector<fasta> vecteurFasta(int argc, char **argv);

    // Function to stock struct fasta length vector
    int tailleFasta(std::vector<fasta> vecFasta);

    // Function to check if the number of sequences in the FASTA file is strictly superior to 3
    bool superieurAtrois(std::vector<fasta> vecFasta);

    // Function to check if sequences are aligned
    bool tailleSequence(std::vector<fasta> vecFasta, int tailleVecteur);

    // Function to create a new alignment with no gaps in columns
    std::vector<fasta> ignoreAllGaps(std::vector<fasta> vecFasta, int tailleVecteur);

    // Function to create 1D dynamic which contains elements of struct fasta
    fasta* tableauFasta(std::vector<fasta> vecFasta, int tailleVecteur);
};
#endif