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
    Foam::sootReactingMixture

Description
    Foam::sootReactingMixture

SourceFiles
    sootReactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef sootReactingMixture_H
#define sootReactingMixture_H

#include "speciesTable.H"
#include "chemistryReader.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sootReactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class sootReactingMixture
:
    public reactingMixture<ThermoType>
{
    // Private member data

        //- Table of species composition
        speciesCompositionTable speciesComposition_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sootReactingMixture(const sootReactingMixture&);

        //- Disallow default bitwise assignment
        void operator=(const sootReactingMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    //typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        sootReactingMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~sootReactingMixture()
    {}


    // Member functions

        // //- Read dictionary
        // void read(const dictionary&);

        // label size() const
        // {
        //     return PtrList<Reaction<ThermoType>>::size();
        // }

        // Reaction<ThermoType>& operator[](const label i)
        // {
        //     return PtrList<Reaction<ThermoType>>::operator[](i);
        // }

        // const Reaction<ThermoType>& operator[](const label i) const
        // {
        //     return PtrList<Reaction<ThermoType>>::operator[](i);
        // }

        // //- Table of species composition
        // const speciesCompositionTable& specieComposition() const
        // {
        //     return speciesComposition_;
        // }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sootReactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
