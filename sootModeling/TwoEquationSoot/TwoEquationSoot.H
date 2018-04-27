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
// 3) related to 2), maybe add a reation class that can hold all of the reactions
//    and return the reaction rates. That way it would be easier to change reactions
//    in the future maybe even at runtime.

#ifndef TwoEquationSoot_H
#define TwoEquationSoot_H

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
#include "ODESystem.H"
#include "localEulerDdtScheme.H"


namespace Foam
{

class TwoEquationSoot
{

    // true constants
    const scalar rho_s = 2000.; // soot density [kg/m^3]
    const scalar Na = 6.232e26; //Avagadros number [particles/kMol]
    const scalar sigma = 1.38e-23; // boltzman [J/K]
    const scalar Ca = 9.; //agglomeration constant [-]
    const scalar n_c = 60.; // minimum particle number for agglomeration [-]
    const scalar pi = 3.1415;

    // Const ref to the thermo
    const psiReactionThermo& thermo_;
    // Const ref to the composition
    const basicSpecieMixture& composition_;
    // Const ref to mesh
    const fvMesh& mesh_;
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
    
    // Source for enthalpy equation
    scalarField Qdot_;

    // DEBUG
    Switch storeReactionRates_;
    volScalarField& rGro_;
    volScalarField& rOxO2_;
    volScalarField& rOxOH_;
    volScalarField& rGasH2O_;
    volScalarField& rGasCO2_;
    volScalarField& rAgg_;
    volScalarField& rKronOH_;
    HashTable<scalar> reactionRateTable_;

    // ...........................................................//
    // // Member functions

    // Given the current thermophysical/chemical state
    // in a cell determine the reaction rates of the 
    // reactions in the soot mass fraction equation
    // rates from Kronenburg(2000) and Josephson(2017)
    void ratesOfChange
    (
        scalar& r_growth,
        scalar& r_oxidation_O2,
        scalar& r_oxidation_OH,
        scalar& r_gasification_H2O,
        scalar& r_gasification_CO2,
        scalar& r_agglomeration
    );

    // For a cell perform a single exlicit step
    // with the provided sub time step to calculate the
    // mass fractions at the end of subdt (Y_final). 
    void explicitStep
    (
        const scalarField& Y_initial,
        scalarField& Y_final,
        const label cellNumber,
        const scalar subdt
    );

    // In the likely event that the required explicit timestep
    // leads to negative mass fractions just advance the problem specie to 0.0.
    // Interpolate that change to the other species and the soot
    // variables then find the corresponding time as well (smallSubTime)
    // This function is called by the function exhaustLowSpecie
    void advanceToZero
    (
        const scalar totalTime,
        scalar& smallSubTime,
        const scalarField& Y_initial,
        scalarField& Y_current
    );

    // If there are negative mass fractions after explicit integration
    // This will call advanceToZero (multiple times if need be) to correct that
    // It returns in the references:
    // - the corresponding smallSubtime step needed to make the worst 
    //   problem specie reach zero but not go negative.
    // - the corresponding other species in Y_final. As though they had been 
    //   evolved with timestep smallSubTime.
    void exhaustLowSpecie
    (
        const scalar totalTime,
        scalar& smallSubTime,
        const scalarField& Y_initial,
        scalarField& Y_final
    );

    // Correct the soot model Qdot based on the sources for the
    // species equations and thier enthalpies of formation
    // Very much how the main combustion chemistry model does it.
    void correctQdot();

    // Calculate and update specieSources for cell given time step 
    // and for the specified cell
    void calcSpecieSources
    (
        const scalar dt,
        const label nSubSteps,
        const List<word> frozenSpecieNames,
        const label cellNumber
    );
public: 

    // Constructor
    TwoEquationSoot
    (
        const psiReactionThermo& thermo, 
        const basicSpecieMixture& composition,
        const volScalarField& Ns,
        const fvMesh& mesh,
        volScalarField& rGro,
        volScalarField& rOxO2,
        volScalarField& rOxOH,
        volScalarField& rGasH2O,
        volScalarField& rGasCO2,
        volScalarField& rAgg,
        volScalarField& rKronOH
    );


    // Main function for calculating the sources.
    // Manages substeps of the main time step dt
    // And decides if the simulation is LTS or not
    void updateSources
    (
        const label nSubSteps = 3
    );

    // **************** Access Functions ************************//
    
    // Return the source for the N_s equation
    tmp<volScalarField> sourceN();
    
    // This will be called in all species equations
    // it returns the appropriate sources depending
    // on the specie mass fraction argument.
    tmp<volScalarField> sourceY
    (
        const volScalarField& Y_field
    );

    // Return the enthalpy source from the reactions
    tmp<volScalarField> Qdot();

}; // end of class EulerImlicitSystem

} // End namepace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
