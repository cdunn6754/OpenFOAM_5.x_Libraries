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

namespace Foam
{

class CellState
{

    const psiReactionThermo& thermo_;
    const basicSpecieMixture& composition_;
    const fvMesh& mesh_;
    const scalarField& Ns_;
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


public:

    CellState
    (
        const psiReactionThermo& thermo,
        const basicSpecieMixture& composition,
        const scalarField& Ns,
        const fvMesh& mesh
    );

    // Update the current state based on the cell number
    void updateCellState
    (
        const scalar cellNumber
    );

    // access functions
    HashTable<scalar> frozenSpecieMassFractions();
    HashTable<scalar> frozenSpeciePPressures();
    HashTable<scalar> thermoProperties();
    scalar Ysoot();
    scalar Nsoot();
    scalar cellVolume();


}; // end class CellState

} // end namespace Foam

#endif
