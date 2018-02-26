/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "ternaryAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(ternaryAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            ternaryAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::ternaryAbsorptionEmission::ternaryAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    model1_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model1"), mesh)
    ),
    model2_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model2"), mesh)
    ),
    model3_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model3"), mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::ternaryAbsorptionEmission::~ternaryAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::ternaryAbsorptionEmission::aCont(const label bandI) const
{
    return model1_->aCont(bandI) + model2_->aCont(bandI)
        + model3_->aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::ternaryAbsorptionEmission::aDisp(const label bandI) const
{
    return model1_->aDisp(bandI) + model2_->aDisp(bandI)
        + model3_->aDisp(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::ternaryAbsorptionEmission::eCont(const label bandI) const
{
    return model1_->eCont(bandI) + model2_->eCont(bandI)
        + model3_->eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::ternaryAbsorptionEmission::eDisp(const label bandI) const
{
    return model1_->eDisp(bandI) + model2_->eDisp(bandI)
        + model3_->eDisp(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::ternaryAbsorptionEmission::ECont(const label bandI) const
{
    return model1_->ECont(bandI) + model2_->ECont(bandI)
        + model3_->ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::ternaryAbsorptionEmission::EDisp(const label bandI) const
{
    return model1_->EDisp(bandI) + model2_->EDisp(bandI)
        + model3_->EDisp(bandI);
}


// ************************************************************************* //
