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

#include "nthOrderDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::nthOrderDevolatilisation<CloudType>::
nthOrderDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    YVolatile0_(volatileData_.size()),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff")))
{
    if (volatileData_.empty())
    {
        WarningIn
        (
            "Foam::nthOrderDevolatilisation<CloudType>::"
            "nthOrderDevolatilisation"
            "("
                "const dictionary& dict, "
                "CloudType& owner"
            ")"
        )   << "Devolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
      Info<< "This is the nth order model! ...Participating volatile species:" <<endl;

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
    }
}


template<class CloudType>
Foam::nthOrderDevolatilisation<CloudType>::
nthOrderDevolatilisation
(
    const nthOrderDevolatilisation<CloudType>& dm
)
:
    DevolatilisationModel<CloudType>(dm),
    volatileData_(dm.volatileData_),
    YVolatile0_(dm.YVolatile0_),
    volatileToGasMap_(dm.volatileToGasMap_),
    residualCoeff_(dm.residualCoeff_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::nthOrderDevolatilisation<CloudType>::
~nthOrderDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::nthOrderDevolatilisation<CloudType>::calculate
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
    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];

	// Fraction of mass0 of this particular species
	// that has the potential to be devoled
	// note in the beggining the particle mass == mass0
        const scalar massFracVolatile0 = YVolatile0_[i];

	// Fraction of current particle mass that has the potential to be devoled
	// for this particular species
	const scalar massFracVolatile = YGasEff[id];

	// New conversion stuff
	
	// This approach assumes that the mass of the char/ash is 
	// unchanged while devolatilisation takes place

	
	// mass of gas initialy in the coal
	const scalar initial_mass_gas = massFracVolatile0 * mass0;

	// Current mass of solid (ash/char) presumed to be constant 
	// in time
	const scalar mass_solid = mass * (YSolidEff[0] + YSolidEff[1]);

	// The mass fraction that remains when all of this species is 
	// devolatilised
	const scalar mf_end_solid = mass_solid/mass0;

	const scalar initial_mass_sg = mass_solid + initial_mass_gas;

	const scalar current_mass_sg = mass_solid + massFracVolatile*mass;

	const scalar mf_still_solid = current_mass_sg/initial_mass_sg;

	// This conversion formula was taken from C3M, its important to base it on the total
	// mass of solid + gas, not just the mass remaining in the solid particle
	const scalar conversion = 1.0 - (mf_still_solid - mf_end_solid)/(1.0 - mf_end_solid);
	

        // Combustion allowed once all volatile components evolved
        done = done && (massFracVolatile <= residualCoeff_*massFracVolatile0);

        // Model coefficients
        const scalar A1 = volatileData_[i].A1(); //[1/s]
        const scalar E = volatileData_[i].E(); //[J/mol]
	const scalar n = volatileData_[i].n(); //[-]
	const scalar R = RR/1000.0; //convert R to [J/mol*k]

        // Kinetic rate
        const scalar kappa = A1*exp(-E/(R*T));

	// Rate
	const scalar rd = kappa * pow( 1.0 - conversion,n);  //rate [1/s]

	// Fraction of mass0 to be transfered (i.e. devolatilised) to carrier gas
	const scalar Y_devol = rd * dt;

        // Mass transferred from particle to carrier gas phase
	dMassDV[id] = min( Y_devol * mass0, massFracVolatile*mass);
    }

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }


}


// ************************************************************************* //
