/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "basicPcReactingMultiphaseCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelPcDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// Stuff added 09-18-17 to get the turbulence stuff too
#include "makeParcelTurbulenceDispersionModels.H"
// More added to incorporate the smithIntrinsic surface reaction model we wrote
#include "makeCoalParcelSmithIntrinsicSurfaceReactionModel.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicPcReactingMultiphaseCloud);

// Kinematic sub-models
makeThermoParcelForces(basicPcReactingMultiphaseCloud);
makeParcelDispersionModels(basicPcReactingMultiphaseCloud);
makeReactingMultiphaseParcelInjectionModels(basicPcReactingMultiphaseCloud);
makeParcelPatchInteractionModels(basicPcReactingMultiphaseCloud);
makeReactingMultiphaseParcelStochasticCollisionModels
(
    basicPcReactingMultiphaseCloud
);
makeParcelSurfaceFilmModels(basicPcReactingMultiphaseCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicPcReactingMultiphaseCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels
(
    basicPcReactingMultiphaseCloud
);
makeReactingParcelPhaseChangeModels(basicPcReactingMultiphaseCloud);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelPcDevolatilisationModels
(
    basicPcReactingMultiphaseCloud
);
makeReactingMultiphaseParcelSurfaceReactionModels
(
    basicPcReactingMultiphaseCloud
);

// Pc dispersion
makeParcelTurbulenceDispersionModels(basicPcReactingMultiphaseCloud);

// smithIntrinsic
makeCoalParcelSmithIntrinsicSurfaceReactionModel(basicPcReactingMultiphaseCloud);


// ************************************************************************* //
