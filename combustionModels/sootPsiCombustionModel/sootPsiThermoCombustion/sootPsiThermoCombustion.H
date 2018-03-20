/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::sootPsiThermoCombustion

Description
    Compressibility-based thermo model wrapper for combustion models

SourceFiles
    sootPsiThermoCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef sootPsiThermoCombustion_H
#define sootPsiThermoCombustion_H

#include "autoPtr.H"
#include "sootPsiCombustionModel.H"
#include "sootPsiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class sootPsiThermoCombustion Declaration
\*---------------------------------------------------------------------------*/

class sootPsiThermoCombustion
:
    public sootPsiCombustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        sootPsiThermoCombustion(const sootPsiThermoCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const sootPsiThermoCombustion&);


protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<sootPsiReactionThermo> thermoPtr_;


public:

    // Constructors

        //- Construct from components
        sootPsiThermoCombustion
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~sootPsiThermoCombustion();


    // Member Functions

        //- Return access to the thermo package
        virtual sootPsiReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const sootPsiReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
