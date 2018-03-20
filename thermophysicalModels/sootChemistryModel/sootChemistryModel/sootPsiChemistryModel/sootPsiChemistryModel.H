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

Class
    Foam::sootPsiChemistryModel

Description
    Chemistry model for compressibility-based thermodynamics

SourceFiles
    sootPsiChemistryModelI.H
    sootPsiChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef sootPsiChemistryModel_H
#define sootPsiChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "sootPsiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class sootPsiChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class sootPsiChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        sootPsiChemistryModel(const sootPsiChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const sootPsiChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<sootPsiReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("psi");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        sootPsiChemistryModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        sootPsiChemistryModel(const fvMesh& mesh, const word& phaseName);


    //- Selector
    static autoPtr<sootPsiChemistryModel> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~sootPsiChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline sootPsiReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const sootPsiReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sootPsiChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
