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

#include "AnalyticalSFORDevolatilization.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AnalyticalSFORDevolatilization<CloudType>::
AnalyticalSFORDevolatilization
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    YVolatile0_(volatileData_.size()),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff"))),
    actUnits_(word(this->coeffDict().lookup("ActivationEnergyUnits"))),
    R_(-1.0)
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "Devolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and cloud gas components
        const label idGas = owner.composition().idGas();
        const scalar YGasTot = owner.composition().YMixture0()[idGas];
        const scalarField& YGas = owner.composition().Y0(idGas);
        forAll(volatileData_, i)
        {
            const word& specieName = volatileData_[i].name();
            const label id = owner.composition().localId(idGas, specieName);
            volatileToGasMap_[i] = id;
            YVolatile0_[i] = YGasTot*YGas[id];

            Info<< "    " << specieName << ": particle mass fraction = "
                << YVolatile0_[i] << endl;
        }
        
        // Check the activation energy units
        if (actUnits_ == "PCCL")
        {
            Info << "Using " << actUnits_ << " -> [kcal/(mol K)],"
                << " units for devolatilization rates." 
                << endl;
            // convert from [J/(kmol K)] to [kcal/(mol K)]
            // factors 4184 J = 1 kcal, 1000 mol = 1 kmol
            R_ = RR/(1000 * 4184);

        }
        else if (actUnits_ == "OF")
        {
            Info << "Using " << actUnits_ << " -> [J/(kmol K)],"
                << " units for devolatilization rates." 
                << endl;

            // set local R_ to built_in RR in [J/kmol]
            R_ = RR;
        }
        else // If it was neither option throw and error
        {
            FatalErrorInFunction 
                << "Please select either 'PCCL' -> [kcal/mol] or 'OF' -> [J/kmol] for "
                    << "the ActivationEnergyUnits input." 
                    << "\nYou selected: " << actUnits_ << abort(FatalError);
        }
    }
}


template<class CloudType>
Foam::AnalyticalSFORDevolatilization<CloudType>::
AnalyticalSFORDevolatilization
(
    const AnalyticalSFORDevolatilization<CloudType>& dm
)
:
    DevolatilisationModel<CloudType>(dm),
    volatileData_(dm.volatileData_),
    YVolatile0_(dm.YVolatile0_),
    volatileToGasMap_(dm.volatileToGasMap_),
    residualCoeff_(dm.residualCoeff_),
    actUnits_(dm.actUnits_),
    R_(dm.R_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AnalyticalSFORDevolatilization<CloudType>::
~AnalyticalSFORDevolatilization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::AnalyticalSFORDevolatilization<CloudType>::calculate
(
    const scalar dt,
    const scalar age,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV
) const
{
    bool done = true;

    // Loop through species
    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];
        const scalar massVolatile0 = mass0*YVolatile0_[i];
        const scalar massVolatile = mass*YGasEff[id];

        // Combustion allowed once all volatile components evolved
        done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients
        const scalar A1 = volatileData_[i].A1();
        const scalar E = volatileData_[i].E();

        // Kinetic rate
        const scalar kappa = A1*exp(-E/(R_*T));

        // Analytical integration
        // dm/dt = kappa * m => m(t) = c * exp(kappa * t)
        // To find c we know that m(0) = massVolatile, for this time step
        const scalar c = massVolatile;
        
        const scalar dm = c * Foam::exp(kappa * dt) - massVolatile;

        // Mass transferred from particle to carrier gas phase
        dMassDV[id] = min(dm, massVolatile);
    }

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
