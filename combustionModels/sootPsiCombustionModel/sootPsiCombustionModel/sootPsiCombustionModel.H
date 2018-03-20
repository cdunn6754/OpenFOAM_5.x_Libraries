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

Class
    Foam::sootPsiCombustionModel

Description
    Combustion models for compressibility-based thermodynamics

SourceFiles
    sootPsiCombustionModelI.H
    sootPsiCombustionModel.C
    sootPsiCombustionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef sootPsiCombustionModel_H
#define sootPsiCombustionModel_H

#include "combustionModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "sootPsiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     class sootPsiCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class sootPsiCombustionModel
:
    public combustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        sootPsiCombustionModel(const sootPsiCombustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const sootPsiCombustionModel&);


public:

    typedef sootPsiReactionThermo ReactionThermo;


    //- Runtime type information
    TypeName("sootPsiCombustionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        sootPsiCombustionModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& combustionProperties,
            const word& phaseName
        ),
        (modelType, mesh, combustionProperties, phaseName)
    );


    // Constructors

        //- Construct from components
        sootPsiCombustionModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& combustionProperties,
            const word& phaseName
        );


        //- Selector
        static autoPtr<sootPsiCombustionModel> New
        (
            const fvMesh& mesh,
            const word& combustionProperties=combustionPropertiesName,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~sootPsiCombustionModel();


    // Member Functions

        //- Return access to the thermo package
        virtual sootPsiReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const sootPsiReactionThermo& thermo() const = 0;


    // IO

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
