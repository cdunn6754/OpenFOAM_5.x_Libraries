/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::smithIntrinsicRate

Description
    Adapted from openFoam COxidationIntrinsicRate model. It now matches the 
    intrinsic model "from smith" as detailed in 
    
    Michael Stollinger, Bertrand Naud, Dirk Roekaerts, Nijso Beishizen, and
    Stefan Heinz. Pdf modeling and simulations of pulverized coal combustion
    - part 2 applcation. Combustion and Flame, 160:396–410, 2012.

    There are two issues currently (07-02-17)
    - Stoellinger paper calls for mole fraction in the main equation, here we
        we use the mass fraction of oxygen
    - I had trouble tracking down some of the constant's units. So currently 
        it doenst work out. But that could be fixed if the consants's
	units are correct. In any case I used the numbers from Stoellinger's
	fortran implementation so it is a least the same as that.

\*---------------------------------------------------------------------------*/

#ifndef smithIntrinsicRate_H
#define smithIntrinsicRate_H

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward class declarations
template<class CloudType>
class smithIntrinsicRate;

/*---------------------------------------------------------------------------*\
                   Class smithIntrinsicRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class smithIntrinsicRate
:
    public SurfaceReactionModel<CloudType>
{
    // Private data

        // Model constants

            //- Stoichiometry of reaction []
            const scalar Sb_;

            //- Mass diffusion limited rate constant [kg/m2sPa]
            const scalar Cdiff_;

            //- Mean pore radius [m]
            const scalar rMean_;

            //- Char porosity [] = 1 - rho_apparent/rho_true
            const scalar theta_;

            //- Pre-exponential factor [kg/m2sPa]
            const scalar Ai_;

            //- Activation energy [kJ/mol]
            const scalar Ei_;

            //- Char specific internal area [m2/g]
            const scalar Ag_;

            //- Pore tortuosity []; default to sqrt(2)
            const scalar zeta_;

  // Clint 06-21-17
  // Oxygen diffusion constant parameter [m2/s]
  const scalar D0_;

  // Knudsen diffusion constant parameter [kg*m/K0.5*kmol0.5*s]
  const scalar Ckn_;


        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- O2 position in global list
            label O2GlobalId_;

            //- CO2 positions in global list
            label CO2GlobalId_;


        // Local copies of thermo properties

            //- Molecular weight of C [kg/kmol]
            scalar WC_;

            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;

            //- Formation enthalpy for CO2 [J/kg]
            scalar HcCO2_;


public:

    //- Runtime type information
    TypeName("smithIntrinsicRate");


    // Constructors

        //- Construct from dictionary
        smithIntrinsicRate
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        smithIntrinsicRate
        (
            const smithIntrinsicRate<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<SurfaceReactionModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceReactionModel<CloudType>>
            (
                new smithIntrinsicRate<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~smithIntrinsicRate();


    // Member Functions

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "smithIntrinsicRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
