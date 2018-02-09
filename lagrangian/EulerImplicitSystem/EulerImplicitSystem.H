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

// TODO//
// 1) make the Y_current, etc. scalarField into a hashtable. I think it will make it
//    much more readable and make some other things easier.
// 2) combine the explicit step function for soot and the frozen variables so
//    we only need to calculate the reaction rates once per time step.
//    maybe add them to the cellState.
// 3) related to 2), maybe add a reation class that can hold all of the reactions
//    and return the reation rates. That way it would be easier to change reactions
//    in the future maybe even at runtime.

#ifndef EulerImplicitSystem_H
#define EulerImplicitSystem_H

#include <math.h>
#include "scalar.H"
#include "scalarField.H"
#include "dimensionSet.H"
#include "psiReactionThermo.H"
#include "fvScalarMatrix.H"
#include "label.H"
#include "tmp.H"
#include "HashTable.H"
#include "List.H"
#include "CellState.H"
#include "SootReactions.H"
#include <time.h>
#include <stdio.h>

namespace Foam
{

class EulerImplicitSystem
{

    // true constants
    const scalar rho_s = 2000.; // soot density [kg/m^3]
    const scalar Na = 6.232e26; //Avagadros number [particles/kMol]
    const scalar sigma = 1.38e-23; // boltzman [J/K]
    const scalar Ca = 9.; //agglomeration constant [-]
    const scalar n_c = 60.; // minimum particle number for agglomeration [-]
    const scalar pi = 3.1415; // ratio of the area of a circle and the square of its
                              // radius

    // Const ref to the thermo
    const psiReactionThermo& thermo_;
    // Const ref to the composition
    const basicSpecieMixture& composition_;
    // CellState object 
    CellState cellState_;
    // make the proper dimension sets for the sources
    dimensionSet Y_source_dims; //[kg/(m^3 * s)]
    dimensionSet N_source_dims; // [1/(m^3 * s)]

    // table for relevant specie molecular weights [kg/kmol]
    HashTable<scalar> MW_;

    // Soot source fields for two soot transport equation number density
    scalarField N_source;

    // Sources for the other related species
    HashTable<scalarField> speciesSources;

    // tracks number of cells in which 
    // the newtons method iterations failed 
    // to converge
    scalar newton_method_failures;
    
    // Constants are determined from examination of 
    // soot equation source equations while
    // only considering Y_s and N_s to be varible.
    // Keys are of the form "equation#_term#"
    HashTable<scalar> C_;


    // Set the constants for cellNumber cell 
    // based on current species/rho/T from main fields
    void setConstants(const label& cellNumber);

    // Set cell constants based on the global thermo fields (rho and T)
    // but use provided cell mass fractions.
    // This is used in the inner iterations when the frozen species
    // are updated but the other thermo fields are not.
    void setConstants
    (
        const label& cellNumber,
        const scalar Y_c2h2_,
        const scalar Y_o2_,
        const scalar Y_oh_
    );


    label nEqns() const;

    // The value of the derivatives given Y and N
    // I.E. dY/dt = f(Y,N) with f taken from Kronenburg 
    // CMC paper
    void derivatives
    (
        const scalar& Y, 
        const scalar& N,
        scalarField& derivative
    );

    // Evaluate the newton equations for the euler implicit form of the 
    // two source equations i.e. 'Y^i - Y^(i-1) - dt(dY/dt^i) = 0
    // for a single cell.
    // Assumes
    // Y_i and N_i hold the values at the integration start time.
    // 
    // Y,N - values for Y_soot and N_soot (guessed values most often)
    // dt - Euler implicit integration time step
    // fY, fN - values of the newton method equations evaluated at Y and N
    void newtonEquations
    (
        const scalar& Y,
        const scalar& N,
        const scalar& Y_i,
        const scalar& N_i,
        const scalar dt,
        scalar& fY,
        scalar& fN
    )const ;

    // Find the Jacobian of the newton method system
    // i.e. the derivatives of the functions from
    // newtonEquations
    // again this peforms the calculation for a single cell.
    void jacobian
    (
        const scalar& Y,
        const scalar& N,
        const scalar dt,
        scalarSquareMatrix& J,
        scalarSquareMatrix& invJ
    ) const;

    void invertJacobian
    (
        const scalarSquareMatrix& J,
        scalarSquareMatrix& invJ
    ) const;

    void explicitEuler
    (
        const scalar subdt,
        const scalar Ysoot_current,
        const scalar N_current,
        scalar& Ysoot_1,
        scalar& N_1
    );


    // Integrate over timestep dt.
    // Once Y^(t+dt),N ^(t + dt) are found use them to find the 
    // source terms. Update sources.
    bool newtonMethodStep
    (
        const scalar dt,
        const scalar Y_0,
        const scalar N_0,
        scalar& Y_1,
        scalar& N_1,
        const scalar relTol = 1e-6
    );

    // Given the current thermophysical/chemical state
    // in a cell determine the reaction rates of the 
    // reactions in the soot mass fraction equation
    // from Kronenburg(2000) and Fletcher(2017)
    void newRatesOfChange
    (
        scalar& r_growth,
        scalar& r_oxidation_O2,
        scalar& r_oxidation_OH,
        scalar& r_gasification_H2O,
        scalar& r_gasification_CO2,
        scalar& r_agglomeration
    );

    // Given the current thermophysical/chemical state
    // in a cell determine the reaction rates of the 
    // reactions in the soot mass fraction equation
    // from Kronenburg CMC paper.
    // Note: weird names (e.g. rate_3g) taken from same paper
    void RatesOfChange
    (
    scalar& rate_3n,
    scalar& rate_3g,
    scalar& rate_4,
    scalar& rate_5,
    const label cellNumber
    );


    // For a cell perform a single exlicit step
    // with the provided sub time step to calculate the
    // mass fractions at the end of subdt. 
    void AdvanceFrozenSpecies
    (
        const scalarField& Y_initial,
        scalarField& Y_final,
        const label cellNumber,
        const scalar subdt,
        const bool advanceProductSpecies = false
    );

    // Determine, for this cell, if the any of 
    // frozen species will be entirely consumed.
    // If they are then true is returned. If not 
    // and the advancement is safe then false is returned.
    bool CheckFrozenSpecies
    (
        const scalarField& Y_initial,
        scalarField& Y_final,
        const label cellNumber,
        const scalar subdt
    );

    // Determine an appropriate/safe sub dt step for the 
    // source integration such that the 'frozen' species
    // are not overconsumed.
    // subdt is the initial default time step (1/2 CFD dt)
    // if the substep becomes too small it returns -1
    scalar GetSafeSubdt
    (
        const scalarField& Y_initial,
        scalarField& Y_final,
        const label cellNumber,
        const scalar subdt 
    );

public: 

    // Constructor
    EulerImplicitSystem
    (
        const psiReactionThermo& thermo, 
        const basicSpecieMixture& composition,
        const volScalarField& Ns,
        const fvMesh& mesh
    );


    // Main function for calculating the sources.
    // Manages substeps of the main time step dt
    void updateSources
    (
        const scalar dt,
        const scalarField meshVolumes,
        const label nSubSteps = 2, 
        const scalar relTol = 1e-6
    );

    // **************** Access Functions ************************//
    // Return the number of failures for the most recent run
    scalar get_newton_method_failures();
    
    // Return the source for the N_s equation
    tmp<fvScalarMatrix> sourceN
    (
        const volScalarField& N_soot
    );
    
    // This will be called in all species equations
    // it returns the appropriate sources depending
    // on the specie mass fraction argument.
    tmp<fvScalarMatrix> sourceY
    (
        const volScalarField& Y_field
    );


}; // end of class EulerImlicitSystem

} // End namepace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
