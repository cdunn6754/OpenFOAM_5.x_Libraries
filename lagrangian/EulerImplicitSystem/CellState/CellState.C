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
    const scalarField& Ns,
    const fvMesh& mesh
)
    :
    thermo_(thermo),
    composition_(composition),
    mesh_(mesh),
    Ns_(Ns),
    frozenSpecieMassFractions_(8),
    frozenSpeciePPressures_(4),
    thermoProperties_(3),
    Ysoot_(0.0),
    Nsoot_(0.0),
    W_(0.0),
    cellVolume_(0.0)
{
    frozenSpecieMassFractions_.insert("C2H2", 0.0);
    frozenSpecieMassFractions_.insert("O2", 0.0);
    frozenSpecieMassFractions_.insert("OH", 0.0);
    frozenSpecieMassFractions_.insert("H2", 0.0);
    frozenSpecieMassFractions_.insert("CO", 0.0);
    frozenSpecieMassFractions_.insert("H", 0.0);
    frozenSpecieMassFractions_.insert("CO2", 0.0);
    frozenSpecieMassFractions_.insert("H20", 0.0);

    frozenSpeciePPressures_.insert("O2", 0.0);
    frozenSpeciePPressures_.insert("OH", 0.0);
    frozenSpeciePPressures_.insert("H2O", 0.0);
    frozenSpeciePPressures_.insert("CO2", 0.0);

    thermoProperties_.insert("rho", 0.0);
    thermoProperties_.insert("T", 0.0);
    thermoProperties_.insert("p",0.0);     
}



void Foam::CellState::updateCellState
(
    const scalar cellNumber
)
{

    // thermo properties
    this->thermoProperties_.set("T",this->thermo_.T().primitiveField()[cellNumber]);
    this->thermoProperties_.set("p",this->thermo_.p().primitiveField()[cellNumber]);
    this->thermoProperties_.set("rho", 
    thermo_.rho().ref().primitiveField()[cellNumber]);
    this->W_ = this->composition_.W().ref()[cellNumber];

    // get the mass fractions and partial pressures
    // where P_i = X_i * P and X_i = (W * Y_i)/W_i
    forAllIter(HashTable<scalar>, this->frozenSpecieMassFractions_, iter)
    {
        word key_ = iter.key();
        scalar& massFraction = iter();

        // set the mass fraction for this specie
        massFraction = this->composition_.Y(key_).primitiveField()[cellNumber];

        // if the specie is also one we need a pp for do that now
        if (this->frozenSpeciePPressures_.found(key_))
        {
            scalar molFraction = (this->W_ * massFraction) /
                this->composition_.W(this->composition_.species()[key_]);

            scalar PPressure = molFraction * this->thermoProperties_["p"];

            frozenSpeciePPressures_.insert(key_, PPressure);
        }
    }
    this->Ysoot_ = this->composition_.Y("SOOT").primitiveField()[cellNumber];
    this->Nsoot_ = this->Ns_[cellNumber];
    
    this->cellVolume_ = this->mesh_.V()[cellNumber];
     
}

Foam::HashTable<scalar> Foam::CellState::frozenSpecieMassFractions()
{
    return frozenSpecieMassFractions_;
}

Foam::HashTable<scalar> Foam::CellState::frozenSpeciePPressures()
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
Foam::scalar Foam::CellState::cellVolume()
{
    return cellVolume_;
}
			
