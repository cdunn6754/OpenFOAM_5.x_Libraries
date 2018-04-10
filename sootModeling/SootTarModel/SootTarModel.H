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

#ifndef SootTarModel_H
#define SootTarModel_H

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
#include "Dictionary.H"
#include "CellState.H"
#include "localEulerDdtScheme.H"


namespace Foam
{

class SootTarModel
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

    // table for relevant specie molecular weights [kg/kmol]
    HashTable<scalar> MW_;

    // table to hold specie mass fractions
    // For each species (i.e. table entry)
    // an initial value and a final value is stored
    // This enables us to check for negative mass fractions
    // following the explicit steps before recording the values.
    HashTable<Pair<scalar>> speciesY_;

    // List of relevant species for tar breakdown
    Dictionary<word> relevantSpecies_;

    // Sources for the other related species
    HashTable<scalarField> speciesSources_;
    
    // Source for enthalpy equation
    scalarField Qdot_;

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
    // species equations and their enthalpies of formation
    // Very similar how the main combustion chemistry model does it.
    void correctQdot();

    // Calculate and update specieSources for cell given time step 
    // and for the specified cell
    void calcSpecieSources
    (
        const scalar dt,
        const label nSubSteps,
        const label cellNumber
    );

public: 

    // Constructor
    SootTarModel
    (
        const psiReactionThermo& thermo, 
        const basicSpecieMixture& composition,
        const volScalarField& Ns,
        const fvMesh& mesh
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
