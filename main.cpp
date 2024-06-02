/*
    Project in C++ done for the 1st year of Master degree in Bioinformatics' intership - University of Montpellier, France (2022-2024)
    Program Align in C++ able to calculate evolutionary distances between amino acids sequences from an aligned FASTA file and create a distance matrice using 5 methods:
        * Distance estimation.
        * Jukes-Cantor model for amino acids.
        * Poisson model for amino acids.
        * Kimura estimation for PAM model.
        * Estimation models for evolutionary distances between amino acids sequences: Poisson Correction and Equal-Input (27 amino acids substitution models).
    
    ./align [evolutionary distances method option] aligned FASTA file [output file option]

    Two options are available for the output file:
        * Output file mat.dist, with a triangular distance matrice in PHYLIP format.
        * Output file seqs.dist, with a distance matrice and other informations: number of sequences, evolutionary distances and header sequences.

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
#include <iomanip>
#include <vector>
#include <string.h>
#include <bits/stdc++.h>

#include "fasta.cpp"
#include "divergence.cpp"
#include "methode.cpp"

using namespace std;

int main(int argc, char** argv){
    // Calculate execution time of the program
    clock_t start, end;
    start = clock();

    Fasta fichier; // Object class Fasta 

    Divergence divergence; // Object class Divergence 

    Methode methode; // Object class Methode 

    struct alphaPC aPC; // Class Methode: Struct alphaPC

    struct alphaEI aEI; // Class Methode: Struct alphaEI

    struct betaEI bEI; // Class Methode: Struct betaEI

    string modele; // Variable for substitution models

    string gaps; // Variable for ignoring gaps in all alignment columns

    bool verifier, superieurTrois, taille; // Variable for boolean in functions: existe, superieurAtrois et tailleSequence 

    vector<fasta> vecFasta; // Variable to stock FASTA file informations in struct fasta vector

    int tailleVecteur, tailleDivergenceObservee; // Variable to stock the number of fasta vector elements and estimate distances vector length 

    fasta* tableauStruct; // Variable to stock 1D array informations of struct fasta 

    vector<double> vecteurDivergenceObservee, vecteurDistancesEvolutives; // Variable to stock distance estimation and evolutionary distances 

    double** matriceDistancesEvolutives; // Variable for distance matrice
    
    double alpha, beta; // Distance estimations parameters alpha and beta
    cout << endl;

    // Print Help manual if program arguments are inferior or equel to 3

    if ((argc <= 3) || (strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)) {
        fichier.usage(argc, argv);
        exit(0);
    }else{
        verifier = fichier.existe(2, argv); // Checking existence of the FASTA file
        cout << "Checking existence of the FASTA file...\n";
        // If FASTA file exists, then stock sequences and headers
        if (verifier != 0){
            cout << "The FASTA file exists.\n";

            vecFasta = fichier.vecteurFasta(2, argv); // Struct fasta vector of FASTA file informations
            cout << endl;
            /*
            Checking if the number of sequences in the FASTA file is strictly superior to 3
            */
            superieurTrois = fichier.superieurAtrois(vecFasta); 

            cout << "The number of amino acids sequences is sufficient to construct distance matrice.\n";
            tailleVecteur = fichier.tailleFasta(vecFasta); // Struct fasta vector length

            cout << endl;
    
            taille = fichier.tailleSequence(vecFasta, tailleVecteur); // Checking if sequences are aligned 

            cout << "Checking sequences alignement...\n";

            if (taille != 0){
                cout << "Amino acids sequences are aligned.\n";
                cout << endl;

                /*
                User must choose to keep or not gaps in sequences alignement.
                If the answer is yes, then alignment is recreated without gaps
                */
                cout << "Ignore gaps in all alignment columns?\n"
                << "Y/n.\n";
                cin >> gaps;
                if ((gaps == "Y") || (gaps == "Yes") || (gaps == "y") || (gaps == "yes") || (gaps == "YES"))
                {
                    cout << "Remove gaps in the alignment.\n";
                    vecFasta = fichier.ignoreAllGaps(vecFasta, tailleVecteur); // Remove gaps in the alignment and creation of the new one
                    taille = fichier.tailleSequence(vecFasta, tailleVecteur); // Stock new length of sequences 
                }else{
                    cout << "Default: keeping gaps.\n"; // By default gaps are keep
                }
                cout << endl;

                tableauStruct = fichier.tableauFasta(vecFasta, tailleVecteur); // Creation of an 1D array with struct fasta elements

                cout << "Calculate distances estimation between sequences...\n";
                vecteurDivergenceObservee = divergence.vecteurDivergences(vecFasta, tailleVecteur); // Creation of distances estimation vector
                tailleDivergenceObservee = divergence.tailleDivergenceObservee(vecteurDivergenceObservee); // Stock distances estimation vector length
                cout << "Distances estimation are calculate.\n";

                cout << endl;
                cout << "Checking evolutionary distances method...\n"; // Creation of evolutionary distances vector in function of the method
                if ((strcmp(argv[1], "-p") == 0) || (strcmp(argv[1], "--poisson") == 0))
                {
                    vecteurDistancesEvolutives = methode.poisson(vecteurDivergenceObservee, tailleDivergenceObservee); // Poisson model for amino acids 
                    cout << "Method: Poisson model for amino acids.\n";
                }else if ((strcmp(argv[1], "-k") == 0) || (strcmp(argv[1], "--kimura") == 0))
                {
                   vecteurDistancesEvolutives = methode.kimura(vecteurDivergenceObservee, tailleDivergenceObservee); // Kimura estimation for PAM model
                   cout << "Method: Kimura estimation for PAM model..\n";
                }else if ((strcmp(argv[1], "-jc") == 0) || (strcmp(argv[1], "--jukescantor") == 0))
                {
                   vecteurDistancesEvolutives = methode.jukesCantor(vecteurDivergenceObservee, tailleDivergenceObservee); // Jukes-Cantor model for amino acids
                   cout << "Method: Jukes-Cantor model for amino acids..\n";
                }else if ((strcmp(argv[1], "-pc") == 0) || (strcmp(argv[1], "--poissoncorrection") == 0)) // Estimation model: Poisson-Correction 
                {
                    cout << "Estimation model: Poisson-Correction.\n";
                    cout << "Please enter the amino acids substitution model: \n";
                    cin >> modele; // User must enter manually the amino acids substitution model between 27 options available
                    // Alpha variable for Poisson-Correction 
                    alpha = aPC.setaPC(modele);
                    beta = 1.00000; // Fixed Beta variable for Poisson-Correction
                    vecteurDistancesEvolutives = methode.estimationGu(vecteurDivergenceObservee, tailleDivergenceObservee, alpha, beta); 
                }else if ((strcmp(argv[1], "-ei") == 0) || (strcmp(argv[1], "--equalinput") == 0)) // Estimation model: Equal-Input
                {
                    cout << "Estimation model: Equal-Input.\n";
                    cout << "Please enter the amino acids substitution model: \n";
                    cin >> modele; // User must enter manually the amino acids substitution model between 27 options available
                    alpha = aEI.setaEI(modele); // Alpha variable for Equal-Input 
                    beta = bEI.setbEI(modele); // Beta variable for Equal-Input 
                    vecteurDistancesEvolutives = methode.estimationGu(vecteurDivergenceObservee, tailleDivergenceObservee, alpha, beta);
                }else{

                    cout << "Default method: Distance estimation.\n";
                    vecteurDistancesEvolutives = vecteurDivergenceObservee; // Default method: Distance estimation
                }
                
                cout << endl;

                matriceDistancesEvolutives = methode.distance(vecteurDistancesEvolutives, tailleVecteur); // Creation of evolutionary distances matrice
                cout << "Creation of evolutionary distances matrice.\n";
                
                /*
                Creation of the output file
                */

                if ((strcmp(argv[3], "-o") == 0) || (strcmp(argv[3], "--output") == 0))
                {
                    divergence.fichierDist(matriceDistancesEvolutives, vecteurDistancesEvolutives, vecFasta, tailleVecteur);
                    cout << "Creation of seqs.dist file (evolutionary distances matrice informations).\n";
                }else if ((strcmp(argv[3], "-m") == 0) || (strcmp(argv[3], "--matrice") == 0))
                {
                    divergence.fichierMat(matriceDistancesEvolutives, vecteurDistancesEvolutives, vecFasta, tailleVecteur);
                    cout << "Creation of mat.dist file (evolutionary distances matrice, PHYLIP format).\n";
                }

                // Destruction of evolutionary distances matrice
                
                for (int i = 0; i < tailleVecteur-1; i++)
                {
                    delete[] matriceDistancesEvolutives[i];
                }
                delete[] matriceDistancesEvolutives;
                
                // Destruction stock 1D array struct fasta 
                
                delete [] tableauStruct;
                }else{
                    cerr << "Error: Amino acids sequences are not aligned.\n";
                    exit(-1);
                }
        }else{
            cerr << "Error: the FASTA file can't be open.\n";
            exit(-1); 
        }
    }
    end = clock(); // End of program
    cout << endl;
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); // Calculate execution time 
    cout << "End of program.\n";
    cout << endl;
    cout << "Execution time: " << fixed
         << time_taken << setprecision(5); // Print execution time
    cout << " secondes " << endl;
    cout << endl;
  return 0;
}