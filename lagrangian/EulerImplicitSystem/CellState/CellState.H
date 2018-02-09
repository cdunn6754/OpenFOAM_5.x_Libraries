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

Class to hold cell thermophysical state for the EulerImplicit integration.

\*---------------------------------------------------------------------------*/

#ifndef CellState_H
#define CellState_H

#include "fvCFD.H"
#include "psiReactionThermo.H"

#include <time.h>
#include <printf.h>

namespace Foam
{

class CellState
{

    const psiReactionThermo& thermo_;
    const basicSpecieMixture& composition_;
    const fvMesh& mesh_;
    const scalarField& Ns_;
    // mixture molar mass field[kg/mol]
    // only changes when the main cfd time step changes
    scalarField molarMassField_;
    HashTable<scalar> frozenSpecieMassFractions_;
    // Partial pressures of species
    HashTable<scalar> frozenSpeciePPressures_;
    HashTable<scalar> thermoProperties_;
    scalar Ysoot_;
    scalar Nsoot_;
    // Mixture molecular weight [kg/kmol]
    scalar W_;
    // Current cell volume [m^3]
    scalar cellVolume_;
    // Cell index
    scalar cellNumber_;


public:

    CellState
    (
        const psiReactionThermo& thermo,
        const basicSpecieMixture& composition,
        const scalarField& Ns,
        const fvMesh& mesh
    );

    // When a new cfd time step is reached we need
    // to update some of the fields stored in this class
    void updateCellStateFields();

    // Update the current state based on the cell number
    // This pulls the data from the main volFields.
    // - Optional updating of the thermo variables
    void updateCellState
    (
        const scalar cellNumber,
        const bool updateThermo = true
    );
    
    // In the sub time steps we need to update just the 
    // - mass fractions of the frozen species
    // - mass fraction of the soot
    // - particle number density
    // In the cell while leaving the cell thermo state as it is since
    // it doesnt evolve during the sub dt's anyway.
    void updateCellState
    (
        const scalarField Y_current,
        const List<word> frozenSpecieNames,
        const scalar Ysoot_current,
        const scalar Nsoot_current
    );

    // access functions
    HashTable<scalar> frozenSpecieMassFractions();
    HashTable<scalar> frozenSpeciePPressures();
    HashTable<scalar> thermoProperties();
    scalar Ysoot();
    scalar Nsoot();
    scalar cellVolume();
    scalar cellNumber();


}; // end class CellState

} // end namespace Foam

#endif
