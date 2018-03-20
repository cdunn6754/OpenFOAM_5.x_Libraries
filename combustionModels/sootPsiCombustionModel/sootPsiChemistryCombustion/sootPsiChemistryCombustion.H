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
    Foam::sootPsiChemistryCombustion

Description
    Compressibility-based chemistry model wrapper for combustion models

SourceFiles
    sootPsiChemistryCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef sootPsiChemistryCombustion_H
#define sootPsiChemistryCombustion_H

#include "autoPtr.H"
#include "sootPsiCombustionModel.H"
#include "sootPsiChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class sootPsiChemistryCombustion Declaration
\*---------------------------------------------------------------------------*/

class sootPsiChemistryCombustion
:
    public sootPsiCombustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        sootPsiChemistryCombustion(const sootPsiChemistryCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const sootPsiChemistryCombustion&);


protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<sootPsiChemistryModel> chemistryPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        sootPsiChemistryCombustion
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& combustionProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~sootPsiChemistryCombustion();


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
