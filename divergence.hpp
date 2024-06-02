/*
    Project in C++ done for the 1st year of Master degree in Bioinformatics' intership - University of Montpellier, France (2022-2024)
    Program Align in C++ able to calculate evolutionary distances between amino acids sequences from an aligned FASTA file and create a distance matrice using 5 methods.

    Class Divergence: Calculcation of distances estimation, creation of the evolutinary distances matrice and output file (seqs.dist or mat.dist).

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

#include "fasta.hpp" // fasta.hpp inclusion to use it's functions (inheritance)

#ifndef DIVERGENCE_HPP
#define DIVERGENCE_HPP

class Divergence
{
public:

  // Divergence Class constructor
  Divergence()
  {
    std::cout << "Divergence Class constructor.\n";
  };

  // Divergence Class destructor
  ~Divergence()
  {
    std::cout << "Divergence Class destructor.\n";
  };

  // Function to calculate distances estimation between two sequences and stock them into a vector
  std::vector<double> vecteurDivergences(std::vector<fasta> vecFasta, int tailleVecteur);

  // Function to stock the lentgh of distances estimation vector
  int tailleDivergenceObservee(std::vector<double> vecteurDivergenceObservee);

  // Function to stock evolutinary distances into a matrice (2D array)
  double **distance(std::vector<double> vecteurDistances, int tailleVecteur);

  /*
    Function to create seqs.dist output file (3rd argument option "-o" or "--output"):
      Number of amino acids sequences
      Triangular matrice of evolutinary distances
      Header sequences
      Compared sequences and evolutinary distances
  */
  std::ofstream fichierDist(double **distanceEvolutive, std::vector<double> vecteurDistances, std::vector<fasta> vecFasta, int tailleVecteur);

  /*
    Function to create mat.dist (3rd argument option "-m" or "--matrice"): 
      Triangular matrice in PHYLIP format
  */
  std::ofstream fichierMat(double **distanceEvolutive, std::vector<double> vecteurDistances, std::vector<fasta> vecFasta, int tailleVecteur);
};

#endif