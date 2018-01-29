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

#include "CellState.H"

// Constructor
Foam::CellState::CellState
(
    const psiReactionThermo& thermo,
    const basicSpecieMixture& composition,
    const scalarField& Ns
)
    :
    thermo_(thermo),
    composition_(composition),
    Ns_(Ns),
    frozenSpecieMassFractions_(6),
    thermoProperties_(3),
    Ysoot_(0.0),
    Nsoot_(0.0)
{
    frozenSpecieMassFractions_.insert("C2H2", 0.0);
    frozenSpecieMassFractions_.insert("O2", 0.0);
    frozenSpecieMassFractions_.insert("OH", 0.0);
    frozenSpecieMassFractions_.insert("H2", 0.0);
    frozenSpecieMassFractions_.insert("CO", 0.0);
    frozenSpecieMassFractions_.insert("H", 0.0);

    thermoProperties_.insert("rho", 0.0);
    thermoProperties_.insert("T", 0.0);
    thermoProperties_.insert("p",0.0);     
}

void Foam::CellState::updateCellState
(
    const scalar cellNumber
)
{
    // get the mass fractions 
    forAllIter(HashTable<scalar>, this->frozenSpecieMassFractions_, iter)
    {
        iter() = this->composition_.Y(iter.key()).primitiveField()[cellNumber];
    }
    Ysoot_ = this->composition_.Y("SOOT").primitiveField()[cellNumber];
    Nsoot_ = this->Ns_[cellNumber];

    // thermo properties
    thermoProperties_.set("T",this->thermo_.T().primitiveField()[cellNumber]);
    thermoProperties_.set("p",this->thermo_.p().primitiveField()[cellNumber]);
    thermoProperties_.set("rho", 
    this->thermo_.rho().ref().primitiveField()[cellNumber]);
}

Foam::HashTable<scalar> Foam::CellState::frozenSpecieMassFractions()
{
    return frozenSpecieMassFractions_;
}

Foam::HashTable<scalar> Foam::CellState::thermoProperties()
{
    return thermoProperties_;
}
Foam::scalar Foam::CellState::Ysoot()
{
    return Ysoot_;
}
Foam::scalar Foam::CellState::Nsoot()
{
    return Nsoot_;
}
			
