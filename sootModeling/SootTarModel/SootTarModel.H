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

    // Convenience typedef
    typedef FixedList<scalar,3> threeList;

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
    // for each specie store a final, current and initial
    // mass fraction. All three are needed to do explicit stepping
    // Organized in FixedList as [initial, current, final]
    HashTable<threeList> speciesY_;

    // Sources for the other related species
    HashTable<scalarField> speciesSources_;
    
    // Source for enthalpy equation
    scalarField Qdot_;

    // List of scalar fields to hold rate references
    scalarField& OxidationRate_;
    scalarField& CrackingRate_;
    scalarField& SootFormationRate_;

    // List to temporarily store rates over sub time steps.
    // Both rate list have the same order
    // Oxidation, cracking and then soot formation
    List<scalarField> tempRateList_;

    // Number of sub steps taken. If integration is going well
    // then these should equal the default, 3, if not then 
    // they will be higher.
    scalarField& subStepsTaken_;
    

    // ...........................................................//
    // // Member functions

    // Given the current thermophysical/chemical state
    // in a cell determine the reaction rates of the 
    // reactions in the soot mass fraction equation
    // rates from Fletcher(1998) and Xu(2016)
    void ratesOfChange
    (
        scalar& r_oxidation,
        scalar& r_gasification,
        scalar& r_sootFormation
    );

    // For a cell perform a single exlicit step
    // with the provided sub time step to calculate the
    // mass fractions at the end of subdt. 
    void explicitStep
    (
        const label cellNumber,
        const scalar subdt
    );

    // Correct the soot model Qdot based on the sources for the
    // species equations and their enthalpies of formation
    // Very similar how the main combustion chemistry model does it.
    void correctQdot();

    // To be used after explicitStep() to ensure that 
    // the reactant mass fractions don't go negative.
    Switch negativeFinalMassFraction();

    // This is called when negativeFinalMassFraction is true. That means
    // that our explicit sub-step has produced a negative mass fraction.
    // This will call findWorstSpecies() in the case where there are two
    // negative mass fractions to determine which became negative first.
    // 
    // Then it will determine the smallTimeStep which would result in the
    // species going to zero only (not negative). It then uses smallTimeStep
    // to evolve the rest of the species and returns the smallTimeStep 
    // to the main integrator loop.
    void exhaustWorstSpecies
    (
        const scalar& subdt, 
        scalar& smallTimeStep
    );

    // Called when there is a negative final mass fraction in speciesY_
    // folowing an excessively large explicit step. This will return the 
    // name of the species that first reaches zero while being integrated.
    word findWorstSpecies();

    // Following the successful comletion of an explicit step
    // update the current mass fraction to equal the final mass fraction
    void updateSpeciesMassFractions();

    // Following the successful completion of an explicit step
    // the speciesY_ table will be updated this function uses
    // it to update the speciesSources table.
    void setSpeciesSources
    (
        const scalar& dt
    );

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
        const fvMesh& mesh,
        scalarField& OxidationRate,
        scalarField& CrackingRate,
        scalarField& SootFormationRate,
        scalarField& subStepsTaken
    );


    // Main function for calculating the sources.
    // Manages substeps of the main time step dt
    // And decides if the simulation is LTS or not
    void updateSources
    (
        const label nSubSteps = 3
    );

    // **************** Access Functions ************************//
    
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
