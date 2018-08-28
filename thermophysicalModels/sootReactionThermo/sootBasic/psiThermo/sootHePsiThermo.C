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

\*---------------------------------------------------------------------------*/

#include "sootHePsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
typename MixtureType::thermoType
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::sootCellMixture
(
    const label celli
) const
{

    // Get list of mass fractions
    const PtrList<volScalarField>& Y_ = MixtureType::Y();

    const scalar Ysoot = MixtureType::Y("SOOT")[celli];
    const label sootIdx = MixtureType::species()["SOOT"];

    // find a good index to start from, not at the soot idx
    label startIdx = 0;
    if (sootIdx == 0)
    {
        startIdx = 1;
    }

    // Start the mixture averaging with first species, not at soot
    typename MixtureType::thermoType mixture = 
        (Y_[startIdx][celli]/(1.0 - Ysoot))*MixtureType::speciesData()[startIdx];

    // Loop through all other species and get the average 
    // properties, excluding soot and whatever startIdx is
    for (label n=(startIdx + 1); n<Y_.size(); n++)
    {
        if (n != sootIdx && n != startIdx)
        {
            // adjusted Y to renomalize mass fractions with soot excluded
            const scalar adjY = Y_[n][celli]/(1.0 - Ysoot);

            mixture += adjY * MixtureType::speciesData()[n];
        }
    }
    
    return mixture;
}

template<class BasicPsiThermo, class MixtureType>
typename MixtureType::thermoType
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::sootPatchFaceMixture
(
    const label patchi,
    const label facei
) const
{

    // Get those mass fractions
    const PtrList<volScalarField>& Y_ = MixtureType::Y();

    const scalar Ysoot = MixtureType::Y("SOOT").boundaryField()[patchi][facei];
    const label sootIdx = MixtureType::species()["SOOT"];

    // find a good index to start from, not at the soot idx
    label startIdx = 0;
    if (sootIdx == 0)
    {
        startIdx = 1;
    }

    // Start the mixture averaging with first species, not at soot
    typename MixtureType::thermoType mixture = 
        (Y_[startIdx].boundaryField()[patchi][facei]/(1.0 - Ysoot))
        * MixtureType::speciesData()[startIdx];

    // Loop through all other species and get the average 
    // properties, excluding soot and whatever startIdx is
    for (label n=(startIdx + 1); n<Y_.size(); n++)
    {
        if (n != sootIdx && n != startIdx)
        {
            // adjusted Y to renomalize mass fractions with soot excluded
            const scalar adjY = Y_[n].boundaryField()[patchi][facei]/(1.0 - Ysoot);

            mixture += adjY * MixtureType::speciesData()[n];
        }
    }

    return mixture;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        // Use the sootCellMixture function to exclude 
        // soot from the calculation
        const typename MixtureType::thermoType sootMixture_ =
            this->sootCellMixture(celli);

        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
        
        psiCells[celli] = sootMixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& sootMixture_ =
                    this->sootPatchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = sootMixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& sootMixture_ =
                    this->sootPatchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = sootMixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::sootHePsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    sootVolume_
    (
        IOobject
        (
            "ThermoSootVolumeFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    sootDensity_
    (
        "sootDensity", 
        dimensionSet(1,-3,0,0,0),
        1000.0
    )
{
    calculate();
    updateSootVolume();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::~sootHePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }
    
    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();
    updateSootVolume();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::updateSootVolume()
{
    
    // Hardcoded soot density
    scalar sootDensity(2000.0); // kg\m^3 from dasgupta thesis

    // To be the sum overall species of [m^3_specie / kg_total]
    scalarField specificVolumeSum = scalarField(this->sootVolume_.size(), 0.0);
    
    // As we iterate we will grab the SOOT specie index
    label sootIdx(-1);

    // Pointer to the mixture for this thermo
    basicSpecieMixture& mixture_ = this->composition();

    forAll(this->Y(), specieI)
    {
        const scalarField& Yi = mixture_.Y()[specieI];
        const word specieName = mixture_.Y()[specieI].name();
          
        if (specieName == "SOOT")
        {
            sootIdx = specieI;
            specificVolumeSum += Yi/sootDensity;
        }
        else
        {
            // loop through cells for non-constant density
            forAll(specificVolumeSum, celli)
            {
                specificVolumeSum[celli] += Yi[celli]/
                    mixture_.rho(specieI, this->p_[celli], this->T_[celli]);
            }
        }
            
    }// end loop through species
    
    // now find and set the soot volume fraction as
    // [V_soot/kg_total] / [V_total/kg_total]
    this->sootVolume_.primitiveFieldRef() = 
        (this->Y()[sootIdx]/sootDensity) / (specificVolumeSum);

    // Set this to 0 so that psi_* P_ is used for the boundary density field.
    this->sootVolume_.boundaryFieldRef() = 0.0;
}

template<class BasicPsiThermo, class MixtureType> 
const Foam::scalarField& 
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::sootVolume() const
{
    return this->sootVolume_;
}


template<class BasicPsiThermo, class MixtureType> 
Foam::tmp<Foam::volScalarField>
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::rho() const
{

    //return this->p_ * this->psi_;

    return (1.0 - this->sootVolume_)*(this->p_*this->psi_) + 
        (this->sootVolume_) * this->sootDensity_;
}

// ************************************************************************* //
