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
#include <math.h> 
#include <iomanip>


#include "divergence.hpp" // divergence.hpp inclusion to use it's functions (inheritance)

#ifndef METHODE_HPP
#define METHODE_HPP

// Struct of Poisson-Correction method "a" parameter for 27 amino acids substitution models (beta fixed to 1)
struct alphaPC
{
  private:
    const double AB = 1.71521;
    const double BLOSUM62 = 3.24334;
    const double cpREV64 = 2.63503;
    const double cpREV = 1.98628;
    const double Dayhoff = 1.99924;
    const double DCMutDayhoff	= 2.01070;
    const double DCMutJTT = 2.55191;
    const double DEN = 2.12834;
    const double FLU = 1.52820;
    const double gcpREV = 1.76147;
    const double HIVb = 1.83588;
    const double HIVw = 1.62839;
    const double JTT = 2.57163;
    const double LG = 2.21046;
    const double mtArt = 0.93628;
    const double mtInv = 1.57997;
    const double mtMam = 0.90348;
    const double mtMet = 1.40469;
    const double mtREV = 1.23867;
    const double mtVer = 1.15596;
    const double mtZoa = 1.05466;
    const double PMB = 3.45924;
    const double rtREV = 2.08011;
    const double stmtREV = 2.03813;
    const double VT = 3.41801;
    const double WAG = 2.69788;
    double WAGstar = 2.80430;
  public :
    std::string modeleSubstitution;

    // Function to give the value for the "a" parameter of Poisson-Correction method.
    double setaPC(std::string modeleSubstitution)
    {
      double aPC; // "a" parameter variable (Poisson-Correction model)
      if ((modeleSubstitution == "AB") || (modeleSubstitution == "ab") || (modeleSubstitution == "Ab"))
      { 
        aPC = AB;
      }else if ((modeleSubstitution == "BLOSUM62") || (modeleSubstitution == "blosum62") || (modeleSubstitution == "Blosum62")){
        aPC = BLOSUM62;
      }else if ((modeleSubstitution == "CPREV64") || (modeleSubstitution == "cprev64") || (modeleSubstitution == "cpREV64")){
        aPC = cpREV64;
      }else if ((modeleSubstitution == "CPREV") || (modeleSubstitution == "cprev") || (modeleSubstitution == "cpREV")){
        aPC = cpREV;
      }else if ((modeleSubstitution == "DCMutDayhoff") || (modeleSubstitution == "Dcmutdayhoff") || (modeleSubstitution == "DcmutDayhoff") || (modeleSubstitution == "dcmutDayhoff") || (modeleSubstitution == "dcmutdayhoff") || (modeleSubstitution == "DCMUTDAYHOFF")){
        aPC = DCMutDayhoff;
      }else if ((modeleSubstitution == "DCMutJTT") || (modeleSubstitution == "Dcmutjtt") || (modeleSubstitution == "DcmutJtt") || (modeleSubstitution == "dcmutJtt") || (modeleSubstitution == "dcmutjtt") || (modeleSubstitution == "DCMUTJTT")){
        aPC = DCMutJTT;
      }else if ((modeleSubstitution == "DEN") || (modeleSubstitution == "Den") || (modeleSubstitution == "den")){
        aPC = DEN;
      }else if ((modeleSubstitution == "FLU") || (modeleSubstitution == "Flu") || (modeleSubstitution == "flu")){
        aPC = FLU;
      }else if ((modeleSubstitution == "GCPREV") || (modeleSubstitution == "gcpREV") || (modeleSubstitution == "GcpRev") || (modeleSubstitution == "gcprev")){
        aPC = gcpREV;
      }else if ((modeleSubstitution == "HIVB") || (modeleSubstitution == "Hivb") || (modeleSubstitution == "hivb") || (modeleSubstitution == "HIVb")){
        aPC = HIVb;
      }else if ((modeleSubstitution == "HIVW") || (modeleSubstitution == "Hivw") || (modeleSubstitution == "hivw")|| (modeleSubstitution == "HIVw")){
        aPC = HIVw;
      }else if ((modeleSubstitution == "JTT") || (modeleSubstitution == "Jtt") || (modeleSubstitution == "jtt")){
        aPC = JTT;
      }else if ((modeleSubstitution == "LG") || (modeleSubstitution == "Lg") || (modeleSubstitution == "lg")){
        aPC = LG;
      }else if ((modeleSubstitution == "MTART") || (modeleSubstitution == "mtArt") || (modeleSubstitution == "MTart") || (modeleSubstitution == "mtART") || (modeleSubstitution == "mtart")){
        aPC = mtArt;
      }else if ((modeleSubstitution == "MTINV") || (modeleSubstitution == "mtInv") || (modeleSubstitution == "MTinv") || (modeleSubstitution == "mtINV") || (modeleSubstitution == "mtinv")){
        aPC = mtInv;
      }else if ((modeleSubstitution == "MTMAM") || (modeleSubstitution == "mtMam") || (modeleSubstitution == "MTmam") || (modeleSubstitution == "mtMAM") || (modeleSubstitution == "mtmam")){
        aPC = mtMam;
      }else if ((modeleSubstitution == "MTMET") || (modeleSubstitution == "mtMet") || (modeleSubstitution == "MTmet") || (modeleSubstitution == "mtMET") || (modeleSubstitution == "mtmet")){
        aPC = mtMet;
      }else if ((modeleSubstitution == "MTREV") || (modeleSubstitution == "mtREV") || (modeleSubstitution == "MTrev") || (modeleSubstitution == "Mtrev") || (modeleSubstitution == "mtrev")){
        aPC = mtREV;
      }else if ((modeleSubstitution == "MTVER") || (modeleSubstitution == "mtVer") || (modeleSubstitution == "MTver") || (modeleSubstitution == "mtVER") || (modeleSubstitution == "mtver")){
        aPC = mtVer;
      }else if ((modeleSubstitution == "MTZOA") || (modeleSubstitution == "mtZOA") || (modeleSubstitution == "MTzoa") || (modeleSubstitution == "Mtzoa") || (modeleSubstitution == "mtzoa")){
        aPC = mtZoa;
      }else if ((modeleSubstitution == "PMB") || (modeleSubstitution == "Pmb") || (modeleSubstitution == "pmb")){
        aPC = PMB;
      }else if ((modeleSubstitution == "RTREV") || (modeleSubstitution == "rtREV") || (modeleSubstitution == "RTrev") || (modeleSubstitution == "Rtrev") || (modeleSubstitution == "rtrev")){
        aPC = rtREV;
      }else if ((modeleSubstitution == "STMTREV") || (modeleSubstitution == "stmtREV") || (modeleSubstitution == "STMTrev") || (modeleSubstitution == "Stmtrev") || (modeleSubstitution == "stmtrev")){
        aPC = stmtREV;
      }else if ((modeleSubstitution == "VT") || (modeleSubstitution == "Vt") || (modeleSubstitution == "vt")){
        aPC = VT;
      }else if ((modeleSubstitution == "WAG") || (modeleSubstitution == "Wag") || (modeleSubstitution == "wag")){
        aPC = WAG;
      }else if ((modeleSubstitution == "WAG*") || (modeleSubstitution == "Wag*") || (modeleSubstitution == "WAGSTAR") || (modeleSubstitution == "wag*")  || (modeleSubstitution == "Wagstar") || (modeleSubstitution == "wagstar")){
        aPC = WAGstar;        
      }else{
          std::cout << "Default model for alpha parameter: Dayhoff.\n";
          // By default "a" parameter get value of the Dayhoff model
          aPC = 1.99924;
      }
      std::cout << "Alpha for Poisson Correction method: " << aPC << std::endl;
      std::cout << "Beta for Poisson Correction method: " << std::fixed << std::setprecision(5) << 1.00000 << std::endl;
      return aPC;
    }
};

// Struct of Equal-Input methode "a" parameter for 27 amino acids substitution models
struct alphaEI
{
  private:
    const double ABalpha = 2.78549;
    const double BLOSUM62alpha = 6.32690;
    const double cpREV64alpha = 4.64357;
    const double cpREValpha = 3.14971;
    const double DayhoffAlpha = 3.14582;
    const double DCMutDayhoffAlpha	= 3.16983;
    const double DCMutJTTalpha = 4.36663;
    const double DENalpha = 3.34672;
    const double FLUalpha = 2.22717;
    const double gcpREValpha = 2.72778;
    const double HIVbAlpha = 2.77572;
    const double HIVwAlpha = 2.45611;
    const double JTTalpha = 4.39688;
    const double LGalpha = 3.56820;
    const double mtArtAlpha = 1.35206;
    const double mtInvAlpha = 2.85866;
    const double mtMamAlpha = 1.30527;
    const double mtMetAlpha = 2.34419;
    const double mtREValpha = 1.95601;
    const double mtVerAlpha = 1.91274;
    const double mtZoaAlpha = 1.57251;
    const double PMBalpha = 7.10575;
    const double rtREValpha = 3.30578;
    const double stmtREValpha = 3.77358;
    const double VTalpha = 6.96847;
    const double WAGalpha = 4.81653;
    const double WAGstarAlpha = 5.01598;
    public :
    std::string modeleSubstitution;

    // Function to give the value for the "a" parameter of Equal-Input method.
    double setaEI(std::string modeleSubstitution)
    {
      double aEI; // "a" parameter variable (Equal-Input)
      if ((modeleSubstitution == "AB") || (modeleSubstitution == "ab") || (modeleSubstitution == "Ab"))
      { 
        aEI = ABalpha;
      }else if ((modeleSubstitution == "BLOSUM62") || (modeleSubstitution == "blosum62") || (modeleSubstitution == "Blosum62")){
        aEI = BLOSUM62alpha;
      }else if ((modeleSubstitution == "CPREV64") || (modeleSubstitution == "cprev64") || (modeleSubstitution == "cpREV64")){
        aEI = cpREV64alpha;
      }else if ((modeleSubstitution == "CPREV") || (modeleSubstitution == "cprev") || (modeleSubstitution == "cpREV")){
        aEI = cpREValpha;
      }else if ((modeleSubstitution == "DCMutDayhoff") || (modeleSubstitution == "Dcmutdayhoff") || (modeleSubstitution == "DcmutDayhoff") || (modeleSubstitution == "dcmutDayhoff") || (modeleSubstitution == "dcmutdayhoff") || (modeleSubstitution == "DCMUTDAYHOFF")){
        aEI = DCMutDayhoffAlpha;
      }else if ((modeleSubstitution == "DCMutJTT") || (modeleSubstitution == "Dcmutjtt") || (modeleSubstitution == "DcmutJtt") || (modeleSubstitution == "dcmutJtt") || (modeleSubstitution == "dcmutjtt") || (modeleSubstitution == "DCMUTJTT")){
        aEI = DCMutJTTalpha;
      }else if ((modeleSubstitution == "DEN") || (modeleSubstitution == "Den") || (modeleSubstitution == "den")){
        aEI = DENalpha;
      }else if ((modeleSubstitution == "FLU") || (modeleSubstitution == "Flu") || (modeleSubstitution == "flu")){
        aEI = FLUalpha;
      }else if ((modeleSubstitution == "GCPREV") || (modeleSubstitution == "gcpREV") || (modeleSubstitution == "GcpRev") || (modeleSubstitution == "gcprev")){
        aEI = gcpREValpha;
      }else if ((modeleSubstitution == "HIVB") || (modeleSubstitution == "Hivb") || (modeleSubstitution == "hivb") || (modeleSubstitution == "HIVb")){
        aEI = HIVbAlpha;
      }else if ((modeleSubstitution == "HIVW") || (modeleSubstitution == "Hivw") || (modeleSubstitution == "hivw")|| (modeleSubstitution == "HIVw")){
        aEI = HIVwAlpha;
      }else if ((modeleSubstitution == "JTT") || (modeleSubstitution == "Jtt") || (modeleSubstitution == "jtt")){
        aEI = JTTalpha;
      }else if ((modeleSubstitution == "LG") || (modeleSubstitution == "Lg") || (modeleSubstitution == "lg")){
        aEI = LGalpha;
      }else if ((modeleSubstitution == "MTART") || (modeleSubstitution == "mtArt") || (modeleSubstitution == "MTart") || (modeleSubstitution == "mtART") || (modeleSubstitution == "mtart")){
        aEI = mtArtAlpha;
      }else if ((modeleSubstitution == "MTINV") || (modeleSubstitution == "mtInv") || (modeleSubstitution == "MTinv") || (modeleSubstitution == "mtINV") || (modeleSubstitution == "mtinv")){
        aEI = mtInvAlpha;
      }else if ((modeleSubstitution == "MTMAM") || (modeleSubstitution == "mtMam") || (modeleSubstitution == "MTmam") || (modeleSubstitution == "mtMAM") || (modeleSubstitution == "mtmam")){
        aEI = mtMamAlpha;
      }else if ((modeleSubstitution == "MTMET") || (modeleSubstitution == "mtMet") || (modeleSubstitution == "MTmet") || (modeleSubstitution == "mtMET") || (modeleSubstitution == "mtmet")){
        aEI = mtMetAlpha;
      }else if ((modeleSubstitution == "MTREV") || (modeleSubstitution == "mtREV") || (modeleSubstitution == "MTrev") || (modeleSubstitution == "Mtrev") || (modeleSubstitution == "mtrev")){
        aEI = mtREValpha;
      }else if ((modeleSubstitution == "MTVER") || (modeleSubstitution == "mtVer") || (modeleSubstitution == "MTver") || (modeleSubstitution == "mtVER") || (modeleSubstitution == "mtver")){
        aEI = mtVerAlpha;
      }else if ((modeleSubstitution == "MTZOA") || (modeleSubstitution == "mtZOA") || (modeleSubstitution == "MTzoa") || (modeleSubstitution == "Mtzoa") || (modeleSubstitution == "mtzoa")){
        aEI = mtZoaAlpha;
      }else if ((modeleSubstitution == "PMB") || (modeleSubstitution == "Pmb") || (modeleSubstitution == "pmb")){
        aEI = PMBalpha;
      }else if ((modeleSubstitution == "RTREV") || (modeleSubstitution == "rtREV") || (modeleSubstitution == "RTrev") || (modeleSubstitution == "Rtrev") || (modeleSubstitution == "rtrev")){
        aEI = rtREValpha;
      }else if ((modeleSubstitution == "STMTREV") || (modeleSubstitution == "stmtREV") || (modeleSubstitution == "STMTrev") || (modeleSubstitution == "Stmtrev") || (modeleSubstitution == "stmtrev")){
        aEI = stmtREValpha;
      }else if ((modeleSubstitution == "VT") || (modeleSubstitution == "Vt") || (modeleSubstitution == "vt")){
        aEI = VTalpha;
      }else if ((modeleSubstitution == "WAG") || (modeleSubstitution == "Wag") || (modeleSubstitution == "wag")){
        aEI = WAGalpha;
      }else if ((modeleSubstitution == "WAG*") || (modeleSubstitution == "Wag*") || (modeleSubstitution == "WAGSTAR") || (modeleSubstitution == "wag*")  || (modeleSubstitution == "Wagstar") || (modeleSubstitution == "wagstar")){
        aEI = WAGstarAlpha;        
      }else{
        std::cout << "Default model for alpha parameter: Dayhoff.\n";
        // By default "a" parameter get value of the Dayhoff model
        aEI = 3.14582;
      }
    std::cout << "Alpha for Equal-input method: " << aEI << std::endl;
    return aEI;
  }
};

// Struct of Equal-Input "b" parameter for 27 amino acids substitution models
struct betaEI
{
  private:
    const double ABbeta = 0.93407;
    const double BLOSUM62beta = 0.94151;
    const double cpREV64beta = 0.93948;
    const double cpREVbeta = 0.93916;
    const double DayhoffBeta = 0.93993;
    const double DCMutDayhoffBeta	= 0.93993;
    const double DCMutJTTbeta = 0.94193;
    const double DENbeta = 0.94143;
    const double FLUbeta = 0.94110;
    const double gcpREVbeta = 0.93745;
    const double HIVbBeta = 0.94179;
    const double HIVwBeta = 0.93819;
    const double JTTbeta = 0.94191;
    const double LGbeta = 0.94051;
    const double mtArtBeta = 0.92743;
    const double mtInvBeta = 0.92211;
    const double mtMamBeta = 0.92473;
    const double mtMetBeta = 0.92546;
    const double mtREVbeta = 0.92467;
    const double mtVerBeta = 0.92052;
    const double mtZoaBeta = 0.92686;
    const double PMBbeta = 0.94195;
    const double rtREVbeta = 0.94024;
    const double stmtREVbeta = 0.92778;
    const double VTbeta = 0.94092;
    const double WAGbeta = 0.94055;
    const double WAGstarBeta = 0.94055;
  public:
  // Function to give the value for the "b" parameter of Equal-Input model.
   double setbEI(std::string modeleSubstitution)
    {
      double bEI; // "b" parameter variable (Equal-Input)
      if ((modeleSubstitution == "AB") || (modeleSubstitution == "ab") || (modeleSubstitution == "Ab"))
      { 
        bEI = ABbeta;
      }else if ((modeleSubstitution == "BLOSUM62") || (modeleSubstitution == "blosum62") || (modeleSubstitution == "Blosum62")){
        bEI = BLOSUM62beta;
      }else if ((modeleSubstitution == "CPREV64") || (modeleSubstitution == "cprev64") || (modeleSubstitution == "cpREV64")){
        bEI = cpREV64beta;
      }else if ((modeleSubstitution == "CPREV") || (modeleSubstitution == "cprev") || (modeleSubstitution == "cpREV")){
        bEI = cpREVbeta;
      }else if ((modeleSubstitution == "DCMutDayhoff") || (modeleSubstitution == "Dcmutdayhoff") || (modeleSubstitution == "DcmutDayhoff") || (modeleSubstitution == "dcmutDayhoff") || (modeleSubstitution == "dcmutdayhoff") || (modeleSubstitution == "DCMUTDAYHOFF")){
        bEI = DCMutDayhoffBeta;
      }else if ((modeleSubstitution == "DCMutJTT") || (modeleSubstitution == "Dcmutjtt") || (modeleSubstitution == "DcmutJtt") || (modeleSubstitution == "dcmutJtt") || (modeleSubstitution == "dcmutjtt") || (modeleSubstitution == "DCMUTJTT")){
        bEI = DCMutJTTbeta;
      }else if ((modeleSubstitution == "DEN") || (modeleSubstitution == "Den") || (modeleSubstitution == "den")){
        bEI = DENbeta;
      }else if ((modeleSubstitution == "FLU") || (modeleSubstitution == "Flu") || (modeleSubstitution == "flu")){
        bEI = FLUbeta;
      }else if ((modeleSubstitution == "GCPREV") || (modeleSubstitution == "gcpREV") || (modeleSubstitution == "GcpRev") || (modeleSubstitution == "gcprev")){
        bEI = gcpREVbeta;
      }else if ((modeleSubstitution == "HIVB") || (modeleSubstitution == "Hivb") || (modeleSubstitution == "hivb") || (modeleSubstitution == "HIVb")){
        bEI = HIVbBeta;
      }else if ((modeleSubstitution == "HIVW") || (modeleSubstitution == "Hivw") || (modeleSubstitution == "hivw")|| (modeleSubstitution == "HIVw")){
        bEI = HIVwBeta;
      }else if ((modeleSubstitution == "JTT") || (modeleSubstitution == "Jtt") || (modeleSubstitution == "jtt")){
        bEI = JTTbeta;
      }else if ((modeleSubstitution == "LG") || (modeleSubstitution == "Lg") || (modeleSubstitution == "lg")){
        bEI = LGbeta;
      }else if ((modeleSubstitution == "MTART") || (modeleSubstitution == "mtArt") || (modeleSubstitution == "MTart") || (modeleSubstitution == "mtART") || (modeleSubstitution == "mtart")){
        bEI = mtArtBeta;
      }else if ((modeleSubstitution == "MTINV") || (modeleSubstitution == "mtInv") || (modeleSubstitution == "MTinv") || (modeleSubstitution == "mtINV") || (modeleSubstitution == "mtinv")){
        bEI = mtInvBeta;
      }else if ((modeleSubstitution == "MTMAM") || (modeleSubstitution == "mtMam") || (modeleSubstitution == "MTmam") || (modeleSubstitution == "mtMAM") || (modeleSubstitution == "mtmam")){
        bEI = mtMamBeta;
      }else if ((modeleSubstitution == "MTMET") || (modeleSubstitution == "mtMet") || (modeleSubstitution == "MTmet") || (modeleSubstitution == "mtMET") || (modeleSubstitution == "mtmet")){
        bEI = mtMetBeta;
      }else if ((modeleSubstitution == "MTREV") || (modeleSubstitution == "mtREV") || (modeleSubstitution == "MTrev") || (modeleSubstitution == "Mtrev") || (modeleSubstitution == "mtrev")){
        bEI = mtREVbeta;
      }else if ((modeleSubstitution == "MTVER") || (modeleSubstitution == "mtVer") || (modeleSubstitution == "MTver") || (modeleSubstitution == "mtVER") || (modeleSubstitution == "mtver")){
        bEI = mtVerBeta;
      }else if ((modeleSubstitution == "MTZOA") || (modeleSubstitution == "mtZOA") || (modeleSubstitution == "MTzoa") || (modeleSubstitution == "Mtzoa") || (modeleSubstitution == "mtzoa")){
        bEI = mtZoaBeta;
      }else if ((modeleSubstitution == "PMB") || (modeleSubstitution == "Pmb") || (modeleSubstitution == "pmb")){
        bEI = PMBbeta;
      }else if ((modeleSubstitution == "RTREV") || (modeleSubstitution == "rtREV") || (modeleSubstitution == "RTrev") || (modeleSubstitution == "Rtrev") || (modeleSubstitution == "rtrev")){
        bEI = rtREVbeta;
      }else if ((modeleSubstitution == "STMTREV") || (modeleSubstitution == "stmtREV") || (modeleSubstitution == "STMTrev") || (modeleSubstitution == "Stmtrev") || (modeleSubstitution == "stmtrev")){
        bEI = stmtREVbeta;
      }else if ((modeleSubstitution == "VT") || (modeleSubstitution == "Vt") || (modeleSubstitution == "vt")){
        bEI = VTbeta;
      }else if ((modeleSubstitution == "WAG") || (modeleSubstitution == "Wag") || (modeleSubstitution == "wag")){
        bEI = WAGbeta;
      }else if ((modeleSubstitution == "WAG*") || (modeleSubstitution == "Wag*") || (modeleSubstitution == "WAGSTAR") || (modeleSubstitution == "wag*")  || (modeleSubstitution == "Wagstar") || (modeleSubstitution == "wagstar")){
        bEI = WAGstarBeta;        
      }else{
          std::cout << "Default model for beta parameter: Dayhoff.\n";
          // Par défaut le paramètre b du modèle de Dayhoff
          bEI = 0.93993;
      }
      std::cout << "Beta for Equal-input method: " << bEI << std::endl;
      return bEI;
    }
};

class Methode : public Divergence
{
  private:
    std::vector<double> vecteurDivergenceObservee; // Estimate distances vector

    int tailleDivergence; // Estimate distances length vector
  
  public:
    // Methode Class constructor
    Methode()
    {
      std::cout << "Methode Class constructor.\n";
    };

    // Methode Class destructor
    ~Methode()
    {
      std::cout << "Methode Class destructor.\n";
    };

    // Function to calculate evolutinary distances with Poisson model for amino acids (options "-p" or "--poisson")
    std::vector<double> poisson(std::vector<double> divergenceObservee, int tailleDivergence);

    // Function to calculate evolutinary distances with Kimura estimation for PAM model (options "-k" or "--kimura")
    std::vector<double> kimura(std::vector<double> divergenceObservee, int tailleDivergence);

    // Function to calculate evolutinary distances with Jukes-Cantor model for amino acids (options "-jc" or "--jukescantor")
    std::vector<double> jukesCantor(std::vector<double> divergenceObservee, int tailleDivergence);

    // Function to calculate evolutinary distances with estimation models (Poisson Correction or Equal-Input) (options "-pc" or "--poissoncorection" / "-ei" or "--equal-input")
    std::vector<double> estimationGu(std::vector<double> divergenceObservee, int tailleDivergence, double alpha, double beta);
};
#endif