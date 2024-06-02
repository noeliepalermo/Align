/*
    Project in C++ done for the 1st year of Master degree in Bioinformatics' intership - University of Montpellier, France (2022-2024)
    Program Align in C++ able to calculate evolutionary distances between amino acids sequences from an aligned FASTA file and create a distance matrice using 5 methods.

    Class Divergence: Calculcate distances estimation, creation of the evolutinary distances matrice and output file (seqs.dist or mat.dist).

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
#include <iomanip>
#include <math.h>
#include <iomanip>

#include "divergence.hpp"

using namespace std;

// Function to calculate distances estimation between two amino acids sequences and stock them in a vector
vector<double> Divergence::vecteurDivergences(vector<fasta> vecFasta, int tailleVecteur)
{
    // Variable to count the number of substitutions between two sequences
    double substitution = 0;

    double match = 0;
    double total = 0;
    // Vector to stock temporarily distances estimation
    vector<double> vecteurDivergenceObservee;

    // Variable for distances estimation: p = n/l, n number of substitution and l number of compared homologous sites
    double divergence = 0;

    // Variable for indels, gaps in sequences
    double indel = 0;

    // Variable to compare homologous sites (length of comapred sequence - indels)
    double site = 0;

    // Iterator to calculate distances estimation
    int j = 0;
    int i = 1;

    // Sequences comparison by "matrix" order (ex. between A/B A/C A/D A/E then B/C, B/D, B/E and D/E)
    while (j <= tailleVecteur-1){
        for (i = 1; i < tailleVecteur-j; i++){
            string seq1 = vecFasta[i+j].s; // Iteration through sequences 1
            string seq2 = vecFasta[j].s; // Iteration through sequences 2
            int tailleSeq1 = seq1.length(); // Variable of sequence 2 length
            for (int k = 0; k <= tailleSeq1; k++) 
            {
                // If compared sites are different not gaps or unknwon amino acids "X": they are substitution
                if ((seq1[k] != seq2[k]) && (seq2[k] != seq1[k]) && ((seq1[k] != '-') && (seq2[k] != '-'))  && ((seq1[k] != 'X') && (seq2[k] != 'X')) && ((seq1[k] != 'x') && (seq2[k] != 'x'))) 
                {
                    substitution++;
                }
                // If compared sites are different gaps or unknwon amino acids "X": they are count in the indels variable
                else if ((seq1[k] == '-') || (seq2[k] == '-') || (seq1[k] == 'X') || (seq2[k] == 'X') || (seq1[k] == 'x') || (seq2[k] == 'x'))  
                {
                    indel++;
                }else{
                // Else sites are the same and ignored
                    continue;
                }   
            }
            // Number of compared homologous sites
            site = tailleSeq1 - indel;
            // Calculate distances estimation
            divergence = substitution/site;
            // If distances estimation are not between the same sequence, they are not stocked
            vecteurDivergenceObservee.push_back(divergence);

            // Reset variables
            substitution = 0;
            indel = 0;
            divergence =0;
            site =0;
            match = 0;
        }
        j++;
    }
   return vecteurDivergenceObservee; // Return distances estimation vector
}

// Function to stock the lentgh of distances estimation vector
int Divergence::tailleDivergenceObservee(vector<double> vecteurDivergenceObservee)
{
    int tailleDivergenceObservee = vecteurDivergenceObservee.size();
    return tailleDivergenceObservee; // Return the lentgh of distances estimation vector
}

// Function to stock evolutinary distances in a matrice (2D array)
double** Divergence::distance(vector<double> distances, int tailleVecteur)
{
    // Matrice size: lentgh of vector struct fasta - 1
    int tailleMatrice = tailleVecteur - 1;

    // Initialisation of the 2D array
    double **distanceEvolutive = new double *[tailleMatrice];

    for(int i = 0; i < tailleMatrice; i++)
    {
        distanceEvolutive[i] = new double[tailleMatrice]{}; // Initialisate matrice at 0
    }
    // If no memory can be allocated, program exit
    if (distanceEvolutive == NULL)
    {
        cerr << "Error: impossible to allocate memory.\n";
        exit(0);
    }
    // Variable to stop the iteration and get a triangular matrice
    int passage = 0;
    // Variable to iterate into evolutinary distances matrice
    int vecteur = 0;

    // Function to create (upper) evolutinary distances triangular matrice
    for(int i = 0; i < tailleMatrice; i++)
    {
        for (int j = 0+passage; j < tailleMatrice; j++)
        {
            // Copy evolutinary distances vector values in upper matrice
            distanceEvolutive[i][j] = distances[vecteur++];
        }
        passage++;
    }  
   return distanceEvolutive; // Return evolutinary distances matrice
}

/*
    Function to create seqs.dist output file (3rd argument option "-o" or "--output"):
        Number of amino acids sequences
        Triangular matrice of evolutinary distances
        Header sequences
        Compared sequences and evolutinary distances
*/
ofstream Divergence::fichierDist(double** distanceEvolutive, vector<double> vecteurDistances, vector<fasta> vecFasta, int tailleVecteur)
{
    ofstream fichier("seqs.dist");
    if (fichier.is_open())
    {
        fichier << "#distances order: d(1,2),...,d(1,n) <new line> d(2,3),...,d(2,n) <new line>...\n"
        // Number of amino acids sequences
        << tailleVecteur;
        fichier << "\n";

        // Triangular matrice of evolutinary distances
        for(int i = 0; i < tailleVecteur - 1; i++)
        {
            for (int j = 0; j < tailleVecteur - 1; j++)
            {
                if (i <= j){
                    cout << fixed;
                    // Decimal value to 6
                    cout << setprecision(6) << showpoint;
                    fichier << distanceEvolutive[i][j] << " ";
                }
            }
            fichier << "\n";
        }

       // Vector to stock pairwise headers
       vector<string> paires;

        // Vector to stock only headers
        vector<string> enteteSeule;

        // Vector to stock pairwise headers (calculate evolutinary distances between two sequences)
        int m = 0;
        int n = 1;
        while (m <= tailleVecteur-1){
            // Stock the first word of the header separate by a blank space (" ")
            int partie2 = vecFasta[m].e.find(" ");
            string entete2 = vecFasta[m].e.substr(0, partie2);
            enteteSeule.push_back(entete2+ "");
            // Stock only header
            fichier << enteteSeule[m] << " "; 

            // Creation of pairwise headers in function of the order of compared sequences
            for (n = 1; n < tailleVecteur-m; n++){
                int partie1 = vecFasta[n+m].e.find(" ");
                string entete1 = vecFasta[n+m].e.substr(0, partie1);
                paires.push_back(entete2+","+entete1+": ");
            }
            m++;
        }

        fichier << "\n";
        fichier << "\n";
        fichier << "#pairwise distances\n";

        // Match between comapared sequences and evolutinary distances
        for (int i = 0; i < paires.size(); ++i) 
        {
            fichier << paires[i] << vecteurDistances[i] << "\n";    
        }
        fichier.close();
    }else{
        cout << "The file can't be write\n";
    }
    // Close file
    fichier.close();
    return fichier; // Return seqs.dist file
}


/*
    Function to create mat.dist (3rd argument option "-m" or "--matrice"): 
        Triangular matrice in PHYLIP format
*/
ofstream Divergence::fichierMat (double** distanceEvolutive, vector<double> vecteurDistances, vector<fasta> vecFasta, int tailleVecteur)
{
    ofstream fichier("mat.dist");
    if (fichier.is_open())
    {
        // Number of amino acids sequences
        fichier << tailleVecteur << "\n";

        // New triangular matrice (upper)
        double temporaire[tailleVecteur][tailleVecteur]{};

        // Variable to stop the iteration and get new triangular matrice
        int vecteur = 0;

        // Variable to iterate into evolutinary distances vector
        int passage = 0;

        for(int i = 0; i < tailleVecteur; i++)
        {
            for (int k = 0+passage; k < tailleVecteur; k++)
            {
                if (i < k)
                {
                    // Copy values of the evolutinary distances vector into upper matrice
                    temporaire[i][k] = vecteurDistances[vecteur++];
                }
            }
            passage++;
        }

        // Creation of the matrice in PHYLIP format
        for(int i = 0; i < tailleVecteur; i++)
        {
            // Stock headers
            if ( i <= tailleVecteur -1){
                // Stock the first word of the header separate by a blank space (" ")
                int entete = vecFasta[i].e.find(" ");
                string subEntete = vecFasta[i].e.substr(0, entete );
                // Stock only header
                fichier << subEntete << " ";
            }

            for (int k = 0; k < tailleVecteur; k++)
            {
                if (k >= i)
                {
                    // Copy values of the upper triangular matrice to the bottom
                    temporaire[k][i] =  temporaire[i][k];
                }
                fichier << fixed;
                fichier << setprecision(6) << showpoint;
                // Matice in PHYLIP format with diagonal equal to 0
               fichier << temporaire[i][k] << "\t";
            }
            fichier << "\n";
        }
        fichier << "\n";
        fichier.close();
    }else{
        cout << "The file can't be write\n";
    }
    // Close file
    fichier.close();
    return fichier; // Return mat.dist file
}
