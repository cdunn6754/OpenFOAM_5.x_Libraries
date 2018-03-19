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

Class
    Foam::sootPureMixture

Description
    Foam::sootPureMixture

SourceFiles
    sootPureMixture.C

\*---------------------------------------------------------------------------*/

#ifndef sootPureMixture_H
#define sootPureMixture_H

#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sootPureMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class sootPureMixture
:
    public basicMixture
{
    // Private data

        ThermoType mixture_;

        // - Species data (just a dummy so that we can 
        //  build call speciesData() so that we can in turn
        //  build this lower level library.
        PtrList<ThermoType> speciesData_;

        // - (just a dummy so that we can 
        //  build call speciesData() so that we can in turn
        //  build this lower level library.
        PtrList<volScalarField> Y_;

        //- Construct as copy (not implemented)
        sootPureMixture(const sootPureMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        sootPureMixture(const dictionary&, const fvMesh&, const word&);


    // Member functions

        const ThermoType& mixture() const
        {
            return mixture_;
        }

        const ThermoType& cellMixture(const label) const
        {
            return mixture_;
        }

        const ThermoType& patchFaceMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        const ThermoType& cellVolMixture
        (
            const scalar,
            const scalar,
            const label
        ) const
        {
            return mixture_;
        }

        const ThermoType& patchFaceVolMixture
        (
            const scalar,
            const scalar,
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        //- Return the raw specie thermodynamic data
        // dummy function so that we can build this library
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Return the raw specie thermodynamic data
        // dummy function so that we can build this library
        const PtrList<volScalarField>& Y() const
        {
            return Y_;
        }

        //- Dummy function to make the basic library compile
        scalar rho
        (
            const label speciei,
            const scalar p,
            const scalar T
        )
        {
            return 0;
        }

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "sootPureMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
