/*
    Project in C++ done for the 1st year of Master degree in Bioinformatics' intership - University of Montpellier, France (2022-2024)
    Program Align in C++ able to calculate evolutionary distances between amino acids sequences from an aligned FASTA file and create a distance matrice using 5 methods.

    Class Methode: Calculcation evolutinary distances using one of the 5 methods: Distance estimation, Jukes-Cantor model for amino acids, Poisson model for amino acids,
    Kimura estimation for PAM model, estimation models (Poisson Correction and Equal-Input).

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
#include <iomanip>
#include <math.h> 

#include "methode.hpp"

using namespace std;

// Function to calculate evolutinary distances with Poisson model for amino acids (options "-p" or "--poisson")
vector<double> Methode::poisson(vector<double> divergenceObservee, int tailleDivergence)
{
    // Vector to stock evolutinary distances
    vector<double> p;

    // Evolutinary distances variable
    double distance = 0;

    for (int i = 0; i < tailleDivergence; i++) 
    {
        // t = -ln(1-p)
        distance = -log(1.0-divergenceObservee[i]);
        p.push_back(distance);
        distance = 0;
    }
    return p;
}

// Function to calculate evolutinary distances with Kimura estimation for PAM model (options "-k" or "--kimura")
vector<double> Methode::kimura(vector<double> divergenceObservee, int tailleDivergence)
{
    // Vector to stock evolutinary distances
    vector<double> k;

    // Evolutinary distances variable
    double distance = 0;

    for (int i = 0; i < tailleDivergence; i++) 
    {
        // t = -ln(1-p-0.2*p²)
        distance = -log(1.0-divergenceObservee[i]-0.2*(pow(divergenceObservee[i],2.0)));
        k.push_back(distance);
        distance = 0;
    }
    return k;
}

// Function to calculate evolutinary distances with Jukes-Cantor model for amino acids (options "-jc" or "--jukescantor")
vector<double> Methode::jukesCantor(vector<double> divergenceObservee, int tailleDivergence)
{
     // Vector to stock evolutinary distances
    vector<double> jc;

    // Evolutinary distances variable
    double distance = 0;

    for (int i = 0; i < tailleDivergence; i++) 
    {
        // t = -19/20*n(1-20/19*p)
        distance = -19.0/20.0*log(1.0-20.0/19.0*divergenceObservee[i]);
        jc.push_back(distance);
        distance = 0;
    }
    return jc;
}

// Function to calculate evolutinary distances with estimation models (Poisson Correction or Equal-Input) (options "-pc" or "--poissoncorection" / "-ei" or "--equal-input")
vector<double> Methode::estimationGu(vector<double> divergenceObservee, int tailleDivergence, double alpha, double beta)
{
    // Vector to stock evolutinary distances
    vector<double> gu;

    // Evolutinary distances variable
    double distance = 0;

    // Alpha and Beta variables
    double a, b;
    for (int i = 0; i < tailleDivergence; i++) 
    {
        // t = a*b*((1-p/b)^-1a - 1)
        distance = alpha*beta*(pow(1-divergenceObservee[i]/beta, -1/alpha)-1);
        gu.push_back(distance);
        distance = 0;
    }
    return gu;
}