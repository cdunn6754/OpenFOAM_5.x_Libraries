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

// Constructors
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
    molarMassField_(Ns.size(),0.0),
    densityField_(Ns.size(),0.0),
    frozenSpecieMassFractions_(8),
    frozenSpeciePPressures_(4),
    thermoProperties_(3),
    Ysoot_(0.0),
    Nsoot_(0.0),
    W_(0.0),
    cellVolume_(0.0),
    cellNumber_(0.0),
    cellCenter_(0.0,0.0,0.0)

{

    frozenSpecieMassFractions_.insert("C2H2", 0.0);
    frozenSpecieMassFractions_.insert("O2", 0.0);
    frozenSpecieMassFractions_.insert("OH", 0.0);
    frozenSpecieMassFractions_.insert("H2", 0.0);
    frozenSpecieMassFractions_.insert("CO", 0.0);
    frozenSpecieMassFractions_.insert("H", 0.0);
    frozenSpecieMassFractions_.insert("CO2", 0.0);
    frozenSpecieMassFractions_.insert("H2O", 0.0);

    frozenSpeciePPressures_.insert("O2", 0.0);
    frozenSpeciePPressures_.insert("OH", 0.0);
    frozenSpeciePPressures_.insert("H2O", 0.0);
    frozenSpeciePPressures_.insert("CO2", 0.0);

    thermoProperties_.insert("rho", 0.0);
    thermoProperties_.insert("T", 0.0);
    thermoProperties_.insert("p",0.0);

    // Calculate the density and mixture molecular weight fields
    updateCellStateFields();
}


Foam::CellState::CellState
(
    const psiReactionThermo& thermo,
    const basicSpecieMixture& composition,
    const scalarField& Ns,
    const fvMesh& mesh,
    const wordList speciesMassFractions,
    const wordList speciesPartialPressures
)
    :
    thermo_(thermo),
    composition_(composition),
    mesh_(mesh),
    Ns_(Ns),
    molarMassField_(Ns.size(),0.0),
    densityField_(Ns.size(),0.0),
    frozenSpecieMassFractions_(8),
    frozenSpeciePPressures_(4),
    thermoProperties_(3),
    Ysoot_(0.0),
    Nsoot_(0.0),
    W_(0.0),
    cellVolume_(0.0),
    cellNumber_(0.0),
    cellCenter_(0.0,0.0,0.0)

{

    // Instantiate mass fraction and partial pressure table entries
    forAll(speciesMassFractions, specieIdx)
    {
        frozenSpecieMassFractions_.insert
            (
                speciesMassFractions[specieIdx],
                0.0
            );
    }

    forAll(speciesPartialPressures, specieIdx)
    {
        frozenSpeciePPressures_.insert
            (
                speciesPartialPressures[specieIdx],
                0.0
            );
    }

    // These thermo properties are always included
    thermoProperties_.insert("rho", 0.0);
    thermoProperties_.insert("T", 0.0);
    thermoProperties_.insert("p",0.0);

    // Calculate the density and mixture molecular weight fields
    updateCellStateFields();
}


void Foam::CellState::updateCellStateFields()
{
    molarMassField_ = composition_.W().ref().primitiveField();
    densityField_ = thermo_.rho().ref().primitiveField();
}

void Foam::CellState::updateCellState
(
    const scalar cellNumber,
    const bool updateThermo
)
{

    if (updateThermo)
    {
        // thermo properties
        this->thermoProperties_.set("T",this->thermo_.T().primitiveField()[cellNumber]);
        this->thermoProperties_.set("p",this->thermo_.p().primitiveField()[cellNumber]);
        this->thermoProperties_.set("rho", this->densityField_[cellNumber]);
        //thermo_.rho().ref().primitiveField()[cellNumber]);
        this->W_ = this->molarMassField_[cellNumber];
        // not really thermo but this shouldn't change either.
        // not going to work with dynamic mesh I guess.
        this->cellVolume_ = this->mesh_.V()[cellNumber];
        this->cellNumber_ = cellNumber;
        this->cellCenter_ = this->mesh_.C().internalField()[cellNumber];
    }
    
    // Get the mass fractions and partial pressures
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

            frozenSpeciePPressures_.set(key_, PPressure);

        }

    }
    this->Ysoot_ = this->composition_.Y("SOOT").primitiveField()[cellNumber];
    this->Nsoot_ = this->Ns_[cellNumber];
}

void Foam::CellState::updateCellState
(
    const scalarField Y_current,
    const List<word> frozenSpecieNames,
    const scalar Ysoot_current,
    const scalar Nsoot_current
)
{
 
    // loop through the frozen species
    forAll(frozenSpecieNames, specieIndex)
    {
        word specieName = frozenSpecieNames[specieIndex];
        scalar Y_specie = Y_current[specieIndex];
        // maybe add some error handling just in case the specieName
        // isnt int the mass fractions hash table.
        scalar& Y_cellState = this->frozenSpecieMassFractions_[specieName];
        
        // update mass fraction
        Y_cellState = Y_specie;

        // if we track pp on this specie update that now too.
        if (this->frozenSpeciePPressures_.found(specieName))
        {
            scalar molFraction = (this->W_ * Y_specie) /
                this->composition_.W(this->composition_.species()[specieName]);
            
            scalar PPressure = molFraction * this->thermoProperties_["p"];

            frozenSpeciePPressures_.set(specieName, PPressure);
        }
    }
    
    // Finally do the soot species too.
    this->Ysoot_ = Ysoot_current;
    this->Nsoot_ = Nsoot_current;
    
}

void Foam::CellState::updateCellState
(
    const HashTable<threeList>& specieY_
)
{
    forAllConstIter(HashTable<threeList>, specieY_, iter)
    {
        const threeList mf = iter();
        const word name = iter.key();

        frozenSpecieMassFractions_[name] = mf[2];

        // Not implementing the partial pressure and mol fraction stuff
        // because we dont use that in SootTarModel.
    }
}

Foam::HashTable<scalar> Foam::CellState::frozenSpecieMassFractions()
{
    return frozenSpecieMassFractions_;
}

Foam::HashTable<scalar> Foam::CellState::frozenSpeciePPressures()
{
    return frozenSpeciePPressures_;
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
Foam::scalar Foam::CellState::cellNumber()
{
    return cellNumber_;
}
Foam::vector Foam::CellState::cellCenter()
{
    return cellCenter_;
}
			
