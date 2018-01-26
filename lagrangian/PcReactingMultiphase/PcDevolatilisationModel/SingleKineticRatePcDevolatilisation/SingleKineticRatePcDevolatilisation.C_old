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

#include "SingleKineticRatePcDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SingleKineticRatePcDevolatilisation<CloudType>::
SingleKineticRatePcDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PcDevolatilisationModel<CloudType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    YVolatile0_(volatileData_.size()),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff"))),
    Ydaf0_(1.0),
    YdafInfTar_(0.24),
    Atar_(382.0),
    Etar_(7.88),
    AtarDy_(3.41e+8),
    EtarDy_(44.0),
    owner_(owner)
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "PcDevolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and cloud gas components
        const label idGas = owner.composition().idGas();
        const scalar YGasTot = owner.composition().YMixture0()[idGas];
        const scalarField& YGas0 = owner.composition().Y0(idGas);
        forAll(volatileData_, i)
        {
            const word& specieName = volatileData_[i].name();
            const label id = owner.composition().localId(idGas, specieName);
            volatileToGasMap_[i] = id;
            YVolatile0_[i] = YGasTot*YGas0[id];

            Info<< "    " << specieName << ": particle mass fraction = "
                << YVolatile0_[i] << endl;
        }

	// Setting the proper value for Ydaf0
	const label idSolid = owner.composition().idSolid();
	const label idLiquid = owner.composition().idLiquid();
	const scalarField& YSolid0 = owner.composition().Y0(idSolid);
	const scalarField& YLiquid0 = owner.composition().Y0(idLiquid);
	const label ashId = owner.composition().localId(idSolid, "ash");
	const label waterId = owner.composition().localId(idLiquid, "H2O");
	const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
	const scalar YLiquidTot = owner.composition().YMixture0()[idLiquid];
	const scalar Yash = YSolidTot * YSolid0[ashId];
	const scalar Ywater = YLiquidTot * YLiquid0[waterId];
	Ydaf0_ = 1.0 - Yash - Ywater;

	ashId_ = ashId;
	waterId_ = waterId;
    }


}


template<class CloudType>
Foam::SingleKineticRatePcDevolatilisation<CloudType>::
SingleKineticRatePcDevolatilisation
(
    const SingleKineticRatePcDevolatilisation<CloudType>& dm
)
:
    PcDevolatilisationModel<CloudType>(dm),
    volatileData_(dm.volatileData_),
    YVolatile0_(dm.YVolatile0_),
    volatileToGasMap_(dm.volatileToGasMap_),
    residualCoeff_(dm.residualCoeff_),
    owner_(dm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SingleKineticRatePcDevolatilisation<CloudType>::
~SingleKineticRatePcDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SingleKineticRatePcDevolatilisation<CloudType>::calculate
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
    scalarField& dMassDV,
    scalarField& tarFields
) const
{
    bool done = true;

    // Initial daf mass of particle
    // needs to incorporate the tar mass 
    // that isnt included
    // so we are basically tacking on a fake 
    // mass here
    const scalar dafMass0 = Ydaf0_ * mass0 / (1.0 - YdafInfTar_);

    // Convert the builtin RR from [J/kmol K] 
    // to [kcal/mol K]
    const scalar pcR = ((RR/1000.)/4184);

    // Find the current daf mass fraction of tar within the particle
    const scalar Yash = YSolidEff[ashId_];
    const scalar Ywater = YLiquidEff[waterId_];
    const scalar YdafTar = (YdafInfTar_ - tarFields[0]);
    const scalar dafMass = (1.0 - Ywater - Yash) * mass /(1.0 - YdafTar);

    


    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];

	// Find the mass of this specie within the particle
	// but we need to account for the mass of tar that
	// is not included in the composition list
	// hence the (1.0 - YdafInfTar_) factor
        const scalar massVolatile0 = dafMass0*YVolatile0_[i];
	// get the current mass of this specie volatiles in the particle
        const scalar massVolatile = dafMass*YGasEff[id];

	// For the PC coal lab devol rate laws we need daf based
	// YdafVolatile0, as opposed to the YVolatile0 we already have.
	// it should also be adjusted to account for the tar fraction
	// in the particle that isnt included in the composition list
	const scalar YdafVolatile0 = (YVolatile0_[i]/Ydaf0_) * (1.0 - YdafInfTar_);

	// Find the percentage of the intial daf mass 
	// that has been devolatilized
	// YDevoled starts at 0.0 and approaches YdafVolatile0
	const scalar massDevoled = massVolatile0 - massVolatile;
	// again the total mass needs to account for the tar fraction
	// that isn't included through the normal mechanisms
	scalar YdafDevoled = massDevoled/(dafMass0);	
 
	Info << "YdafDevoled " <<  volatileData_[i].name() << ": "
	     << YdafDevoled << endl;

	Info << "YdafDevole0 " <<  volatileData_[i].name() << ": "
	     << YdafVolatile0 << endl;
        // Combustion allowed once all volatile components evolved
        done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients
        const scalar Ap = volatileData_[i].Ap();
        const scalar Ep = volatileData_[i].Ep();

	// make sure we dont exceed the amount of volatiles in the coal
	if (YdafDevoled >= YdafVolatile0)
	  {
	    YdafDevoled = YdafVolatile0;
	  }

        // Kinetic rate
        const scalar kappa = Ap*exp(-Ep/(pcR*T));

        // Mass transferred from particle to carrier gas phase
	// The rates are also based on the daf mass
	const scalar YdafTransfered = dt * kappa * (YdafVolatile0 - YdafDevoled);

	dMassDV[id] = min(YdafTransfered*dafMass0, massVolatile);

	

    } // end per specie primary devolatilization



    // --- Secondary Pyrolysis

    // - Get the particle data from tarFields
    // Primary released by particle daf mass fraction
    scalar& Yp = tarFields[0];
    // difference between Primary tar released (Yp)
    // and actual tar remaining (i.e. secondary tar (Ys))
    scalar& dY = tarFields[1];

    // - Devolatilize the primary tar (just like any other species)
    const scalar kappaTar = Atar_ * exp(- Etar_/(pcR*T));
    const scalar YpNew = dt * kappaTar * (YdafInfTar_ - Yp);
     
    // Add to primary tar and decompose tar mass fraction
    Yp += YpNew;


    // -  Now evaluate the dY equation (Y here based on daf)

    // rate constant for dY equation
    const scalar kappaDy = AtarDy_ * exp(- EtarDy_/(pcR*T) );
    const scalar dYnew = dt * kappaDy * (YdafInfTar_ - dY);
    dY += dYnew;
    const scalar dYmass = dYnew * dafMass0;

    // Separate the tar breakdown mass into its contituents
    forAll(volatileData_, j)
      {
	const label ids = volatileToGasMap_[j];
	dMassDV[ids] += volatileData_[j].Yinfs() * (dYmass);
      }

    Info << "TAR: " << Yp << endl;
    Info << "dY: " << dY << endl;

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
