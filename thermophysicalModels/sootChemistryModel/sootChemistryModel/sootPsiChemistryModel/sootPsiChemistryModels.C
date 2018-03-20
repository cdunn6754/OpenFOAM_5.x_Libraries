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

InClass
    Foam::sootPsiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "sootPsiChemistryModel.H"
#include "chemistryModel.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        icoPoly8HThermoPhysics
    );


    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        icoPoly8HThermoPhysics
    );


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        sootPsiChemistryModel,
        icoPoly8EThermoPhysics
    );


    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        sootPsiChemistryModel,
        icoPoly8EThermoPhysics
    );
}

// ************************************************************************* //
