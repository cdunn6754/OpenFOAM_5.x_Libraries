/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EulerImplicitSystem.H"


// Constructor
Foam::EulerImplicitSystem::EulerImplicitSystem
(
    const psiReactionThermo& thermo, 
    const basicSpecieMixture& composition,
    const volScalarField& Ns,
    const fvMesh& mesh
)
    :
    thermo_(thermo),
    composition_(composition),
    cellState_(thermo, composition, Ns, mesh),
    Y_source_dims(1,0,-1,0,0),// Source dimensions (already multipied by volume)
    N_source_dims(0,0,-1,0,0),// Source dimensions (already multipied by volume)
    MW_(9),
    N_source(Ns.size(),0.0),
    speciesSources(9),
    newton_method_failures(0),
    C_(7)
{
    Info<< "Creating Soot Model Solver \n\n" << endl;
    
    // Initialize all specie sources to zero
    speciesSources.insert("SOOT", scalarField(Ns.size(),0.0));
    speciesSources.insert("C2H2", scalarField(Ns.size(),0.0));
    speciesSources.insert("O2", scalarField(Ns.size(),0.0));
    speciesSources.insert("OH", scalarField(Ns.size(),0.0));
    speciesSources.insert("H2", scalarField(Ns.size(),0.0));
    speciesSources.insert("CO", scalarField(Ns.size(),0.0));
    speciesSources.insert("H", scalarField(Ns.size(),0.0));
    speciesSources.insert("CO2", scalarField(Ns.size(),0.0));
    speciesSources.insert("H2O", scalarField(Ns.size(),0.0));

    // set up MW table
    MW_.insert("SOOT", composition.W(composition.species()["SOOT"]));
    MW_.insert("C2H2", composition.W(composition.species()["C2H2"]));
    MW_.insert("O2", composition.W(composition.species()["O2"]));
    MW_.insert("OH", composition.W(composition.species()["OH"]));
    MW_.insert("H2", composition.W(composition.species()["H2"]));
    MW_.insert("CO", composition.W(composition.species()["CO"]));
    MW_.insert("H", composition.W(composition.species()["H"]));
    MW_.insert("CO2", composition.W(composition.species()["CO2"]));
    MW_.insert("H2O", composition.W(composition.species()["H2O"]));

    // set up constants
    C_.insert("A", 0.0);
    C_.insert("1_1", 0.0);
    C_.insert("1_2", 0.0);
    C_.insert("1_3", 0.0);
    C_.insert("1_4", 0.0);
    C_.insert("2_1", 0.0);
    C_.insert("2_2", 0.0);
}



/// Member functions

// // Set the constants based on new time step species/rho/T
// void Foam::EulerImplicitSystem::setConstants(const label& cellNumber) 
// {
//     // Set these as the current field values.
//     const scalar& Y_c2h2_ = this->cellState_.frozenSpecieMassFractions()["C2H2"];
//     const scalar& Y_o2_ = this->cellState_.frozenSpecieMassFractions()["O2"];
//     const scalar& Y_oh_ = this->cellState_.frozenSpecieMassFractions()["OH"];
            
//     setConstants(cellNumber, Y_c2h2_, Y_o2_, Y_oh_);
// }

// void Foam::EulerImplicitSystem::setConstants
// (
//     const label& cellNumber,
//     const scalar Y_c2h2_,
//     const scalar Y_o2_,
//     const scalar Y_oh_
// ) 
// {
//     const scalar& rho_ = this->cellState_.thermoProperties()["rho"];
//     const scalar& T_ = this->cellState_.thermoProperties()["T"];
            
//     // Particle surface area constant, used in both equations
//     this->C_["A"] = rho_ * pi * Foam::pow(6/(pi * rho_s), (2./3.));

//     // Soot mass fraction equation (Y*rho_/MW) = concentrations [kmol/m^3]
//     this->C_["1_1"] = 2. * (6.3e3* Foam::exp(-21000./T_)) * 
//         (Y_c2h2_*rho_/MW_["C2H2"]) * MW_["SOOT"];
//     this->C_["1_2"] = 2. * (7.5e2 * Foam::exp(-12100./T_)) * C_["A"] * 
//         (Y_c2h2_*rho_/MW_["C2H2"]) * MW_["SOOT"];
//     this->C_["1_3"] = (7.15e2 * Foam::sqrt(T_) * Foam::exp(-19800./T_)) * 
//         C_["A"] * (Y_o2_*rho_/MW_["O2"]) * MW_["SOOT"];
//     this->C_["1_4"] = 0.36 * Foam::sqrt(T_) * C_["A"] * 
//         (Y_oh_*rho_/MW_["OH"]) * MW_["SOOT"];

//     // Soot particle number density equation
//     this->C_["2_1"] = 2. * (6.3e3 * Foam::exp(-21000./T_)) *
//         (Y_c2h2_*rho_/MW_["C2H2"]) * (Na/n_c);
//     this->C_["2_2"] = 2. * Ca * Foam::pow(6./(pi * rho_s), (1./6.))
//         * Foam::sqrt((6. * sigma * T_)/rho_s) 
//         * Foam::sqr(rho_);
// }


// Foam::label Foam::EulerImplicitSystem::nEqns() const
// {
//     return 2;
// }


// // The value of the derivatives given Y and N
// // I.E. dY/dt = f(Y,N) with f taken from Kronenburg 
// // CMC paper
// void Foam::EulerImplicitSystem::derivatives
// (
//     const scalar& Y, 
//     const scalar& N,
//     scalarField& derivative
// )
// {
            
//     // Evolution of soot mass fraction
//     derivative[0] = (C_["1_1"] + 
//     (C_["1_2"] - C_["1_3"] - C_["1_4"]) * (Foam::pow(Y,(2./3.)) * Foam::pow(N,(1./3.))));

//     // Evolution of soot number density (particles/kg)
//     derivative[1] = (C_["2_1"] - 
//     C_["2_2"] * (Foam::pow(Y,(1./6.)) * Foam::pow(N,(11./6.))));
// }

// Evaluate the newton equations for the euler implicit form of the 
// two source equations i.e. 'Y^i - Y^(i-1) - dt(dY/dt^i) = 0
// for a single cell.
// Assumes
// Y_0 and N_0 hold the values at the integration start time.
// 
// Y,N - values for Y_soot and N_soot (guessed values most often)
// dt - Euler implicit integration time step
// fY, fN - values of the newton method equations evaluated at Y and N
void Foam::EulerImplicitSystem::newtonEquations
(
    const scalar& Y,
    const scalar& N,
    const scalar& Y_0,
    const scalar& N_0,
    const scalar dt,
    scalar& fY,
    scalar& fN
) const
{
    // Evolution of soot mass fraction y[0] equation
    fY = Y - Y_0 - dt * (C_["1_1"] + 
    (C_["1_2"] - C_["1_3"] - C_["1_4"]) * (Foam::pow(Y,(2./3.)) * 
    Foam::pow(N,(1./3.))));

    // Evolution of soot number density (particles/m^3) equation
    fN = N - N_0 - dt * (C_["2_1"] - 
    C_["2_2"] * (Foam::pow(Y,(1./6.)) * Foam::pow(N,(11./6.))));

}

// // Find the Jacobian of the newton method system
// // i.e. the derivatives of the functions from
// // newtonEquations
// // again this peforms the calculation for a single cell.
// void Foam::EulerImplicitSystem::jacobian
// (
//     const scalar& Y,
//     const scalar& N,
//     const scalar dt,
//     scalarSquareMatrix& J,
//     scalarSquareMatrix& invJ
// ) const
// {   

//     // Versions to be used in denominators
//     const scalar Ysafe = Foam::max(Y + Foam::ROOTVSMALL, Foam::ROOTVSMALL);        
//     const scalar Nsafe = Foam::max(N + Foam::ROOTVSMALL, Foam::ROOTVSMALL);

//     // set the jacobian entries
//     J(0, 0) = 1 - dt*(C_["1_2"] - C_["1_3"] - C_["1_4"]) *
//         (2./3.) * Foam::pow(Ysafe,(-1./3.)) * Foam::pow(N,(1./3.));

//     J(0, 1) = - dt*(C_["1_2"] - C_["1_3"] - C_["1_4"]) *
//         Foam::pow(Y,(2./3.)) * (1./3.) * Foam::pow(Nsafe,(-2./3.));

//     J(1, 0) = - dt * C_["2_2"] * 
//         (1./6.) * Foam::pow(Ysafe,(-5./6.)) * Foam::pow(N,(11./6.));

//     J(1, 1) = 1 - dt * C_["2_2"] * 
//         Foam::pow(Y,(1./6.)) * (11./6.) * Foam::pow(Nsafe,(5./6.));
         
//     // now invert the jacobian
//     this->invertJacobian(J,invJ);

// }

// void Foam::EulerImplicitSystem::invertJacobian
// (
//     const scalarSquareMatrix& J,
//     scalarSquareMatrix& invJ
// ) const
// {
//     // get inverse of 2x2 jacobian with adjugate/determinant method
//     scalar detJ = J(0,0)*J(1,1) - J(0,1)*J(1,0);

//     if (fabs(detJ) < Foam::SMALL)
//     {
//         Info << "Small determinant: " << detJ << endl;
//         detJ = Foam::SMALL;
//     }

//     invJ(0,0) = J(1,1)/detJ;
//     invJ(0,1) = -J(0,1)/detJ;
//     invJ(1,0) = -J(1,0)/detJ;
//     invJ(1,1) = J(0,0)/detJ;
            
// }

// // Perform a newton method integration on a single cell
// // returns the value of the Y_s and N_s at time
// // + dt
// // in Y_1 and N_1
// bool Foam::EulerImplicitSystem::newtonMethodStep
// (
//     const scalar dt,
//     const scalar Y_0,
//     const scalar N_0,
//     scalar& Y_1,
//     scalar& N_1,
//     const scalar relTol
// )
// {
             
//     // set intial values for this cell
//     scalar Y_i = Y_0;
//     scalar N_i = Y_0;
                
//     // N.M. equations value storage
//     scalar fY(0.0);
//     scalar fN(0.0);

//     // Jacobian storage
//     scalarSquareMatrix J(this->nEqns());
//     scalarSquareMatrix invJ(this->nEqns());

//     // temporary variables for Newton's method iterations
//     scalar Y_next(0.0);
//     scalar N_next(0.0);
//     label converged(0);

//     label newton_iteration(0);

//     // Newton method iterations
//     while( !converged && (newton_iteration < 10)) 
//     {
//         // get jacobian and inverse
//         this->jacobian(Y_i, N_i, dt, J, invJ);
                  
//         // newton method function values
//         this->newtonEquations(Y_i, N_i, Y_0, 
//         N_0, dt, fY, fN);

//         // Newton method step
//         Y_next = Y_i - (invJ(0,0) * fY + invJ(0,1) * fN);
//         N_next = N_i - (invJ(1,0) * fY + invJ(1,1) * fN);

//         // Force them to be zero if one is negative
//         if (Y_next <= 0. || N_next <= 0.)
//         {
//             Y_next = 0.0;
//             N_next = 0.0;
//         }
//         // else check for convergence
//         else if 
//             (fabs(Y_next - Y_i)/fabs(Y_next + Foam::SMALL) < relTol 
//             && 
//             fabs(N_next - N_i)/fabs(N_next + Foam::SMALL) < relTol)
//         {
//             converged = 1;
//         }

//         Y_i = Y_next;
//         N_i = N_next;

//         // increment the iteration counter
//         newton_iteration++;
//     }// end Newton method loop

//     // If we did not converge in 10 iterations
//     if (newton_iteration >= 10)
//     {
//         // increment failure counter
//         this->newton_method_failures ++;

//         // if not converged set it as it was previously
//         // bad idea? 
//         Y_i = Y_0;
//         N_i = N_0;

//         return false;

//     }
    
//     // set the references to the 
//     // correct values for return
//     Y_1 = Y_i;
//     N_1 = N_i;
    
//     return true;

// } // end newtonMethodStep function

void Foam::EulerImplicitSystem::newRatesOfChange
(
    scalar& r_growth,
    scalar& r_oxidation_O2,
    scalar& r_oxidation_OH,
    scalar& r_gasification_H2O,
    scalar& r_gasification_CO2,
    scalar& r_agglomeration
)
{
    // get relevant state variables for the cell
    scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    scalar cell_T = this->cellState_.thermoProperties()["T"];
    scalar cell_c2h2 = this->cellState_.frozenSpecieMassFractions()["C2H2"];
    scalar cell_Ys = this->cellState_.Ysoot();
    scalar cell_Ns = this->cellState_.Nsoot();
    // Partial pressures
    scalar P_co2 = this->cellState_.frozenSpeciePPressures()["CO2"];
    scalar P_o2 = this->cellState_.frozenSpeciePPressures()["O2"];
    scalar P_h2o = this->cellState_.frozenSpeciePPressures()["H2O"];
    scalar P_oh = this->cellState_.frozenSpeciePPressures()["OH"];
    // Specie concentrations
    const scalar C_c2h2(cell_rho * cell_c2h2 * (1/this->MW_["C2H2"]));
    // Molar Mass of soot
    const scalar M_s = this->MW_["SOOT"];
    // current cell volume
    const scalar V_ = this->cellState_.cellVolume();


    // Soot particle diameter and surface area for cell
    // (these are considered uniform in the cell)
    scalar cell_dp(0.0);
    scalar cell_As(0.0);
    if (cell_Ns > Foam::SMALL)
    {
        cell_dp = 
            Foam::pow
            (
                (6.*cell_Ys)/(this->pi * this->rho_s * cell_Ns),
                (1./3.)
            );

        cell_As = this->pi * Foam::sqr(cell_dp) *
            cell_rho * cell_Ns;
    }

    // if (cellState_.cellNumber() == 72)
    // {
    //     Info << "In ROC, As, dp: "<< cell_As << ",  " << cell_dp << nl <<endl;

    //     Info << "Cell Volume: " << cellState_.cellVolume() << nl << endl;
    // }
      
    // Set constants from papers
    // Kronenburg(2000)
    scalar A_growth(7.5e2);
    scalar T_growth(12100.0);
    // Fletcher(2017)
    scalar A_o2(7.98e-1);
    scalar A_oh(1.89e-3);
    scalar A_co2(3.06e-17);
    scalar A_h2o(6.27e4);
    scalar E_o2(1.77e5);
    scalar E_co2(5.56e3);
    scalar E_h20(2.95e5);
    scalar n_(0.13);
    scalar R_(8.314); //[J/(mol*K)]

    // Now calculate the rates
    // Soot Mass fraction [kmol/(m^3*s)]
    if( cellState_.cellNumber() == 48)
    {
        // Info << "0made it here" << endl;
        // Info <<"PP_H2O: " << P_h2o << endl;
        // Info <<"Pressure: " << cellState_.thermoProperties()["p"] << endl;
        // Info <<"Y_h2o: " << cellState_.frozenSpecieMassFractions()["H2O"] << endl;
        
    }

    r_growth = A_growth * Foam::exp(-T_growth/cell_T) * cell_As * C_c2h2;
    r_oxidation_O2 = Foam::pow(cell_T, -0.5) * 
        A_o2 * P_o2 * Foam::exp(-E_o2/(R_*cell_T)) *
        (cell_As/(V_ * M_s));
    r_oxidation_OH = Foam::pow(cell_T, -0.5) * A_oh * P_oh * (cell_As/(V_ * M_s));
        
    r_gasification_CO2 = A_co2 * Foam::pow(P_co2,0.5) * Foam::sqr(cell_T) *
        Foam::exp(-E_co2/(R_*cell_T)) * (cell_As/(V_ * M_s));
       
    r_gasification_H2O = A_h2o * Foam::pow(P_h2o,n_) * Foam::pow(cell_T, -0.5) *
        Foam::exp(-E_h20/(R_ * cell_T)) * (cell_As/(V_ * M_s));
    r_agglomeration = 2 * this->Ca * Foam::sqrt(cell_dp) * 
        Foam::sqrt((6.0 * this->sigma * cell_T)/(this->rho_s)) * 
        Foam::sqr(cell_rho * cell_Ns);


    
    
} //end newRatesOfChange


// void Foam::EulerImplicitSystem::RatesOfChange
// (
//     scalar& rate_3n,
//     scalar& rate_3g,
//     scalar& rate_4,
//     scalar& rate_5,
//     const label cellNumber
// )
// {
//     // get relevant state variables for the cell
//     const scalar& cell_rho = this->cellState_.thermoProperties()["rho"];
//     const scalar& cell_T = this->cellState_.thermoProperties()["T"];
//     const scalar& cell_c2h2 = this->cellState_.frozenSpecieMassFractions()["C2H2"];
//     const scalar& cell_o2 = this->cellState_.frozenSpecieMassFractions()["O2"];
//     const scalar& cell_oh = this->cellState_.frozenSpecieMassFractions()["OH"];
//     const scalar& cell_Ys = this->cellState_.Ysoot();
//     const scalar& cell_Ns = this->cellState_.Nsoot();

//     // Specie concentrations
//     const scalar C_c2h2(cell_rho * cell_c2h2 * (1/this->MW_["C2H2"]));
//     const scalar C_oh(cell_rho * cell_oh * (1/this->MW_["OH"]));
//     const scalar C_o2(cell_rho * cell_o2 * (1/this->MW_["O2"]));
//     // Soot particle diameter and surface area for cell
//     // (these are considered uniform in the cell)
//     scalar cell_dp(0.0);
//     scalar cell_As(0.0);
//     if (cell_Ns > Foam::SMALL)
//     {
//         cell_dp = 
//             Foam::pow
//             (
//                 (6.*cell_Ys)/(this->pi * this->rho_s * cell_Ns),
//                 (1./3.)
//             );

//         cell_As = this->pi * Foam::sqr(cell_dp) *
//             cell_rho * cell_Ns;
//     }
      
//     // Now actually find the rate constants
//     // First reaction rate constants, names taken from Kronenburg
//     const scalar k_3n = 0.63e4 * Foam::exp(-21000./cell_T);
//     const scalar k_3g = 0.75e3 * Foam::exp(-12100./cell_T);
//     const scalar k_4  = 7.15e2 * Foam::pow(cell_T,0.5) * 
//         Foam::exp(-19800./cell_T);
//     const scalar k_5 = 0.36 * Foam::pow(cell_T,0.5);

//     // Form the reaction rates [kmol/(s*m^3)] (i think)
//     rate_3n = k_3n * C_c2h2;
//     rate_3g = k_3g * C_c2h2 * cell_As;
//     rate_4  = k_4  * C_o2 * cell_As;
//     rate_5 = k_5 * C_oh * cell_As;
// } //end RatesOfChange

void Foam::EulerImplicitSystem::explicitEuler
(
    const scalar subdt,
    const scalar Ysoot_current,
    const scalar N_current,
    scalar& Ysoot_1,
    scalar& N_1
)
{
    // TODO: this is a crappy way to do it, we have to call RatesOfChange
    // for both this and the frozen species, the two functions should be 
    // combined in the future probably. Or we could make the rates a
    // member variable maybe.

    scalar cell_rho = this->cellState_.thermoProperties()["rho"];

    // Soot per volume [kg_soot/m^3]
    const scalar massSoot = cell_rho * Ysoot_current;

    // rates for soot mass fraction reactions [kmol/(m^3 * s)]
    scalar r_growth(0.0); // growth of soot particle from c2h2 not nucleation
    scalar r_oxidation_O2(0.0);
    scalar r_oxidation_OH(0.0);
    scalar r_gasification_H2O(0.0);    
    scalar r_gasification_CO2(0.0);
    // rate for soot particle number density reaction[1/(m^3 * s)]
    scalar r_agglomeration(0.0);
    // calculate the rates
    this->newRatesOfChange(r_growth, r_oxidation_O2,r_oxidation_OH, 
    r_gasification_H2O, r_gasification_CO2, r_agglomeration);

    // total rate of soot consumption for conveneince
    scalar r_soot_consumption = r_oxidation_O2 + r_oxidation_OH + r_gasification_CO2
        + r_gasification_H2O;

    // Take explicit Euler step
    scalar massSootFinal = massSoot + this->MW_["SOOT"] *
        (2*r_growth  - r_soot_consumption) * subdt;

    // convert back to mass fraction
    Ysoot_1 = massSootFinal/cell_rho;
    
    // Need to divide by rho because this rate is for d(N*rho)/dt
    // and we only want to evolve N
    N_1 = N_current - (1.0/cell_rho) * r_agglomeration * subdt;


    // Force them to be positive
    if (Ysoot_1 <= 0. || N_1 <= 0.)
    {
        Ysoot_1 = 0.0;
        N_1 = 0.0;
    }
                
}

void Foam::EulerImplicitSystem::AdvanceFrozenSpecies
(
    const scalarField& Y_initial,
    scalarField& Y_final,
    const label cellNumber,
    const scalar subdt,
    const bool advanceProductSpecies
)
{
    // density in this cell
    const scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    // get relevant state variables for the cell
    // start with per volume mass in cell [kg/m^3]
    const scalar massC2H2(Y_initial[0]*cell_rho);
    const scalar massO2(Y_initial[1]*cell_rho);
    const scalar massOH(Y_initial[2]*cell_rho);
    const scalar massCO2(Y_initial[6]*cell_rho);
    const scalar massH2O(Y_initial[7]*cell_rho);
    
    // rates for soot mass fraction reactions [kmol/(m^3 * s)]
    scalar r_growth(0.0); // growth of soot particle from c2h2 not nucleation
    scalar r_oxidation_O2(0.0);
    scalar r_oxidation_OH(0.0);
    scalar r_gasification_H2O(0.0);    
    scalar r_gasification_CO2(0.0);
    // rate for soot particle number density reaction[1/(m^3 * s)] (not used here)
    scalar r_agglomeration(0.0);
    // calculate the rates
    this->newRatesOfChange(r_growth, r_oxidation_O2,r_oxidation_OH, 
    r_gasification_H2O, r_gasification_CO2, r_agglomeration);
    
    // Explicit time step
    // Reactant species
    scalar massC2H2Final = massC2H2 - this->MW_["C2H2"]*r_growth* subdt;
    scalar massO2Final = massO2 - this->MW_["O2"]*r_oxidation_O2*subdt;
    scalar massOHFinal = massOH - this->MW_["OH"]*r_oxidation_OH*subdt;
    scalar massCO2Final = massCO2 - this->MW_["CO2"]*r_gasification_CO2*subdt;
    scalar massH2OFinal = massH2O - this->MW_["H2O"]*r_gasification_H2O*subdt;


    // Assume the overall cell density has not changed 
    // convert back to mass fractions
    Y_final[0]= massC2H2Final / cell_rho;
    Y_final[1]= massO2Final / cell_rho;
    Y_final[2]= massOHFinal / cell_rho;

    // Product species
    if (advanceProductSpecies)
    {
        // get the mass /volume of the product species
        const scalar massH2(Y_initial[3]*cell_rho);
        const scalar massCO(Y_initial[4]*cell_rho);
        const scalar massH(Y_initial[5]*cell_rho);
        
        
        // Add the mass/volume from the current explict step
        // (it's always a source for theses species)
        scalar massH2Final = massH2 + this->MW_["H2"]*(r_growth + r_gasification_H2O)
            *subdt;
        scalar massCOFinal = massCO + this->MW_["CO"]*(r_oxidation_OH +
            r_gasification_H2O + 2*r_gasification_CO2) *subdt;
        scalar massHFinal = massH + this->MW_["H"]*(r_oxidation_OH)*subdt;
        // CO2 is also a reactant in gasification
        massCO2Final += this->MW_["CO2"]*(r_oxidation_O2)*subdt;
        

        // Convert back to mass fraction
        Y_final[3] = massH2Final/cell_rho;
        Y_final[4] = massCOFinal/cell_rho;
        Y_final[5] = massHFinal/cell_rho;
        Y_final[6] = massCO2Final/cell_rho;
        Y_final[7] = massH2OFinal/cell_rho;
    }

    // DEBUG
    // if (cellNumber == 72)
    // {
    //     Info << "Initial Mass Fractions: \n" << 
    //         cellState_.frozenSpecieMassFractions() << 
    //         "\nRates: " << r_growth << ", " <<
    //         r_oxidation_O2 << ", " <<
    //         r_oxidation_OH << ", " <<
    //         r_gasification_H2O << ", " <<
    //         r_gasification_CO2 << ", " << endl;

    //     Info << "Pressure, temp, density and PP: " << cellState_.thermoProperties() 
    //         << "\n"
    //         << cellState_.frozenSpeciePPressures() << endl;

    //     Info << "Time step: " << subdt << nl << endl;

    //     Info << "Final H2O mass fraction: " << Y_final[7] << endl;

    //     FatalErrorInFunction
    //         << "Inside reactions"
    //             << abort(FatalError);
    // }

}// end AdvanceFrozenspecies

// bool Foam::EulerImplicitSystem::CheckFrozenSpecies
// (
//         const scalarField& Y_initial,
//         scalarField& Y_final,
//         const label cellNumber,
//         const scalar subdt
// )
// {

//     // Flag to indicate species are over-consumed
//     bool overConsumed = false;

//     // Use an explicit step to determine the
//     // frozen specie mass fractions at the end of the time step
//     this->AdvanceFrozenSpecies(Y_initial, Y_final, cellNumber, subdt);

//     // reactant mass fractions should always decrease in this context
//     scalarField massFractionConsumed = Y_final - Y_initial;
    
//     // make sure that for each of the three species
//     // a significant mass fractions remains
//     forAll(Y_final, specie)
//     {
//         if (Y_final[specie] < 0.0)
//         {
//             overConsumed = true;
//         }
//     }
//     return overConsumed;
// }// end CheckFrozenSpecies

// Foam::scalar Foam::EulerImplicitSystem::GetSafeSubdt
// (
//     const scalarField& Y_initial,
//     scalarField& Y_final,
//     const label cellNumber,
//     const scalar subdt 
// )
// {
//     // safe sub dt to be determined and returned.
//     scalar safeSubdt(subdt);

//     // presume that the current subdt is not safe.
//     bool dtTooLarge = true;

//     while (dtTooLarge)
//     {
//         // See if the current time step is safe
//         dtTooLarge = this->CheckFrozenSpecies(Y_initial,Y_final,cellNumber,safeSubdt);
        
//         // if it doesnt work, then halve it
//         if (dtTooLarge)
//         {
//             safeSubdt *= 0.5;
//         }

//         // If the substep becomes ridiculously small 
//         // we will just return -1.
//         // The thought is that in this scenario there is so little
//         // mass fraction (for one of the specie at least)
//         // that we can't integrate without
//         // losing it all then just skip integrating at all this cell/timestep.
//         // NOTE: This needs to be improved but
//         // I think this happens very rarely if ever. We should do 
//         // speciated time steps and take the lowest probably.
//         if (safeSubdt <= subdt / 32.0 || safeSubdt < 1e-10)
//         {
//             safeSubdt = -1;
//             break;
//         }
//     }

//     return safeSubdt;
// }

//****************** Public member functions ********************//

Foam::tmp<Foam::fvScalarMatrix> Foam::EulerImplicitSystem::sourceN
(
    const volScalarField& N_soot
)
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(N_soot, N_source_dims));

    // get reference to the fvmatrix inside the tmp
    fvScalarMatrix& fvm = tfvm.ref();

    fvm.source() = -this->N_source;

    // reset field to zero in preparation for next time step.
    this->N_source = 0.0;

    return tfvm;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::EulerImplicitSystem::sourceY
(
    const volScalarField& Y_field
)
{
    word specieName = Y_field.name();

    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Y_field, Y_source_dims));

    // get the reference to fvMatrix that tmp is wrapped around
    fvScalarMatrix& fvm = tfvm.ref();

    if (this->speciesSources.found(specieName))
    {
        fvm.source() = -this->speciesSources[specieName];
        
        // reset field to zero in preparation for next time step.
        this->speciesSources[specieName] = 0.0;
        return tfvm;
    }
    else
    {
        // return empty source because this specie is unrelated to the soot model
        return tfvm;
    }
}

// Main function for calculating the sources.
// Manages substeps of the main time step dt
void Foam::EulerImplicitSystem::updateSources
(
    const scalar dt,
    const scalarField meshVolumes,
    const label nSubSteps,
    const scalar relTol
)
{   


    // Update the fields stored in cellState to new time
    this->cellState_.updateCellStateFields();

    // max (of the cells) number of times 
    // that the subdt needed to be decreased
    label maxNumberOfIters(0);
    // Loop through cells
    forAll(this->thermo_.rho().ref(), cellNumber)
    {
        maxNumberOfIters = 0;
       
        // set default sub time step for this cell
        scalar subdt = dt/nSubSteps;

        // Grab the current field data for this cell
        // this will update the cellState thermo data too
        this->cellState_.updateCellState(cellNumber);
        
        // Local storage of frozen specie mass fractions.
        // the first 3 are reactants in the soot equation and are
        // used when determining a safe sub dt.
        // The second 3 are all products of those reactions.
        // they will be not be used in determining a safe dt
        scalarField Y_current(8,0.0);
        Y_current[0] = this->cellState_.frozenSpecieMassFractions()["C2H2"];
        Y_current[1] = this->cellState_.frozenSpecieMassFractions()["O2"];
        Y_current[2] = this->cellState_.frozenSpecieMassFractions()["OH"];
        Y_current[3] = this->cellState_.frozenSpecieMassFractions()["H2"];
        Y_current[4] = this->cellState_.frozenSpecieMassFractions()["CO"];
        Y_current[5] = this->cellState_.frozenSpecieMassFractions()["H"];
        Y_current[6] = this->cellState_.frozenSpecieMassFractions()["CO2"];
        Y_current[7] = this->cellState_.frozenSpecieMassFractions()["H2O"];

        // Index this list with the names of the species
        List<word> frozenSpecieNames(8,"none");
        frozenSpecieNames[0] = word("C2H2");
        frozenSpecieNames[1] = word("O2");
        frozenSpecieNames[2] = word("OH");
        frozenSpecieNames[3] = word("H2");
        frozenSpecieNames[4] = word("CO");
        frozenSpecieNames[5] = word("H");
        frozenSpecieNames[6] = word("CO2");
        frozenSpecieNames[7] = word("H2O");

        // Store a constant version of these initial mass fractions
        const scalarField Y_initial(Y_current);

        // Create another field for storage of Mass Fractions
        // at the end of the sub time step
        scalarField Y_final(Y_current.size(), 0.0);

        // Reset the frozen specie mass fraction (in this cell)
        // now that the safe sub step size has been determined.
        // Below we will evolve these mass fractions in concert with 
        // changes in the soot mass fraction.
        Y_current = Y_initial;
        Y_final = Y_initial;
        
        // Get the soot variables in this cell at the beginning of the
        // time step and save as const
        const scalar Ysoot_0 = this->cellState_.Ysoot();
        const scalar N_0 = this->cellState_.Nsoot();

        // Intermediate temporary soot variables for storage during
        // sub step time iterations
        scalar Ysoot_current(Ysoot_0);
        scalar N_current(N_0);

        // variables at end of next overall (i.e. cfd) time step
        // initialized to the current values.
        scalar Ysoot_1(Ysoot_0);
        scalar N_1(N_0);

        // bool to specify if the explicit method over shot and we need a 
        // new smaller time step (i.e. negative mass fractions)
        bool explicitStepSuccess(false);
        while(! explicitStepSuccess)
        {
            // loop through advancing the soot variables and 
            // local copy of frozen species 
            // over the sub time step
            scalar totalTime(0.0);
            bool notNegative(true);
            while(totalTime < dt and notNegative)
            {
                // Advance soot variables with explicit step.
                this->explicitEuler(subdt, Ysoot_current, N_current,
                Ysoot_1, N_1);
                
                // Advance, explicitly, local cell copy of frozen species to be used
                // in the next iteration for calculation of rate constants.
                this->AdvanceFrozenSpecies(Y_current, Y_final, cellNumber, subdt, true);

                // Update both the soot variables and the frozen species for
                // the next iteration
                Y_current = Y_final;
                Ysoot_current = Ysoot_1;
                N_current = N_1;

                // before the rates are recalculated on the next iteration update the
                // cell state but only frozen species and soot mass fractions,
                // and soot particle number density. 
                //The thermo/cell volume shouldn't be changing anyway.
                cellState_.updateCellState(Y_current, frozenSpecieNames,
                Ysoot_current, N_current);

                if (min(Y_current) < 0.0)
                {
                    notNegative = false;
                }
            
                totalTime += subdt;
            } // end while totalTime < dt

            // Now check to see if the iterations with the previous time step
            // produced physical results (i.e. didnt over consume reactants)
            if (min(Y_current) < 0.0)
            {
                // decrease time step
                maxNumberOfIters++;
                subdt = subdt/2.0;

                // reset the frozen species and the soot varibles
                Y_current = Y_initial;
                Y_final = Y_initial;
                Ysoot_current = Ysoot_0; 
                Ysoot_1 = Ysoot_0;
                N_current = N_0;
                N_1 = N_0;

                // And re-read the cell state from the main fields
                // (that are not changed by EulerImplicitSystem).
                // Dont bother with the thermo, it hasn't changed anyway.
                this->cellState_.updateCellState(cellNumber, false);
                
                // If the timestep become ridiculously small 
                // exit loop and leave the variables unchanged.
                // Results in sources of zero.
                if (subdt < 1e-10)
                {
                    //Info << "Timestep too small" << endl;
                    explicitStepSuccess = true;
                }
            }
            else
            {
                // If the species mass fractions seem okay then signal success
                explicitStepSuccess = true;
            }
            
        } // end while(!explicitStepSuccess)

        // density in this cell (assumed constant over these steps)
        const scalar cellRho = this->cellState_.thermoProperties()["rho"];


        // Now set the source terms given the changes in soot mass fraction,
        // frozen specie mass fraction  and soot number density. These 
        // source terms will be used when the respective transport
        // equations are solved
        // NOTE: Multiply by density because the transport
        //  equations are d(rho*Y)/dt not just d(Y)/dt, 
        //  we take rho to be constant here so just pull it out.
        const scalar cellVolume(meshVolumes[cellNumber]);
        // Soot number density
        this->N_source[cellNumber] = (N_1 - N_0)*cellVolume*cellRho/dt;
        // Soot mass fraction
        this->speciesSources["SOOT"][cellNumber] = 
            (Ysoot_1 - Ysoot_0)*cellVolume*cellRho/dt;
        // Other species
        // Reactant species 
        // NOTE: Multiply by density because the transport
        //  equations are d(rho*Y)/dt not just d(Y)/dt
        // and rho is constant over this step
        this->speciesSources["C2H2"][cellNumber] = 
            (Y_final[0] - Y_initial[0])*cellVolume*cellRho/dt;
        this->speciesSources["O2"][cellNumber] = 
            (Y_final[1] - Y_initial[1])*cellVolume*cellRho/dt;
        this->speciesSources["OH"][cellNumber] = 
            (Y_final[2] - Y_initial[2])*cellVolume*cellRho/dt;
        this->speciesSources["H2O"][cellNumber] = 
            (Y_final[7] - Y_initial[7])*cellVolume*cellRho/dt;
        // Product species
        this->speciesSources["H2"][cellNumber] = 
            (Y_final[3] - Y_initial[3])*cellVolume*cellRho/dt;
        this->speciesSources["CO"][cellNumber] = 
            (Y_final[4] - Y_initial[4])*cellVolume*cellRho/dt;
        this->speciesSources["H"][cellNumber] = 
            (Y_final[5] - Y_initial[5])*cellVolume*cellRho/dt;
        this->speciesSources["CO2"][cellNumber] = 
            (Y_final[6] - Y_initial[6])*cellVolume*cellRho/dt;

    }// end loop through cells
    Info << "The maximum number of outer iterations is " << maxNumberOfIters << endl;
}

// **************** Access Functions ************************//
