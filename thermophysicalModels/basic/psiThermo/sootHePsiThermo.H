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
    Foam::sootHePsiThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    sootHePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef sootHePsiThermo_H
#define sootHePsiThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sootHePsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class sootHePsiThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{

    // Private Data
    
        //- Soot volume fraction, used to determine multiphase 
        // density.
        scalarField sootVolume_;
    
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        // - Specialized version of cellMixture adapted from multiComponentMixture.
        // This will ignore the soot contribution when making the mixture
        typename MixtureType::thermoType sootCellMixture(const label celli) const;

        //- Construct as copy (not implemented)
        sootHePsiThermo(const sootHePsiThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("sootHePsiThermo");


    // Constructors

        //- Construct from mesh and phase name
        sootHePsiThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~sootHePsiThermo();


    // Member functions

        //- Update properties
        virtual void correct();


    // New 03-2018
    //- Update soot volume fraction
    void updateSootVolume();

    //- Const access to soot volume fraction
    const scalarField& sootVolume() const;

    //- Density calculation that overides virtual psiThermo function.
    //  Treats soot fraction in proper multiphase fashion.
    tmp<volScalarField> rho();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sootHePsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
