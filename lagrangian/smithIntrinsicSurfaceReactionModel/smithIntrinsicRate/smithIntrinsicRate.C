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

\*---------------------------------------------------------------------------*/

#include "smithIntrinsicRate.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::smithIntrinsicRate<CloudType>::smithIntrinsicRate
(
    const dictionary& dict,
    CloudType& owner
)
:
  SurfaceReactionModel<CloudType>(dict, owner, typeName),
  Sb_(this->coeffDict().lookupOrDefault("Sb", 1.33)),
  Cdiff_(this->coeffDict().lookupOrDefault("Cdiff", 5.e-12)),
  rMean_(this->coeffDict().lookupOrDefault("rMean", 6.e-8)),
  theta_(this->coeffDict().lookupOrDefault("theta",0.1)),
  Ai_(this->coeffDict().lookupOrDefault("Ai", 0.052)),
  Ei_(this->coeffDict().lookupOrDefault("Ei", 161.5)),
  Ag_(this->coeffDict().lookupOrDefault("Ag", 104.)),
  zeta_(this->coeffDict().lookupOrDefault("zeta", sqrt(2.0))),
  D0_(this->coeffDict().lookupOrDefault("D0", 3.13e-4)),
  Ckn_(this->coeffDict().lookupOrDefault("Ckn", 97.0)),
  CsLocalId_(-1),
  O2GlobalId_(owner.composition().carrierId("O2")),
  CO2GlobalId_(owner.composition().carrierId("CO2")),
  WC_(0.0),
  WO2_(0.0),
  HcCO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().W(O2GlobalId_);
    const scalar WCO2 = owner.thermo().carrier().W(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);

    if (Sb_ < 0)
    {
        FatalErrorInFunction
            << "Stoichiometry of reaction, Sb, must be greater than zero" << nl
            << exit(FatalError);
    }

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::smithIntrinsicRate<CloudType>::smithIntrinsicRate
(
    const smithIntrinsicRate<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    Sb_(srm.Sb_),
    Cdiff_(srm.Cdiff_),
    rMean_(srm.rMean_),
    theta_(srm.theta_),
    Ai_(srm.Ai_),
    Ei_(srm.Ei_),
    Ag_(srm.Ag_),
    zeta_(srm.zeta_),
    D0_(srm.D0_),
    Ckn_(srm.Ckn_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::smithIntrinsicRate<CloudType>::
~smithIntrinsicRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::smithIntrinsicRate<CloudType>::calculate
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
) const
{
    // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar Ychar = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion until combustible fraction is consumed
    if (Ychar < SMALL)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();

    // Local mass fraction of O2 in the carrier phase []
    const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[celli];

    // Quick exit if oxidant not present
    if (YO2 < ROOTVSMALL)
    {
        return 0.0;
    }

    // Average temperature between particle and gas 
    const scalar Tm = 0.5*(T + Tc);

    // Diffusion rate coefficient [m2/s]
    const scalar Rox = Cdiff_/d*pow(Tm, 0.75);

    // Apparent density of pyrolysis char [kg/m3]
    const scalar rhop = 6.0*mass/(constant::mathematical::pi*pow3(d));

    // Knusden diffusion coefficient [m2/s]
    const scalar Dkn = Ckn_*rMean_*sqrt(T/WO2_);

    // Reference temperature for the oxygen diffusion [K]
    const scalar T0 = 1500 ;

    // Oxygen diffusion coefficient [m2/s]
    const scalar Dox = D0_*pow(Tm/T0, 1.75);

    // Effective diffusion [m2/s]
    const scalar De = theta_/sqr(zeta_)/(1.0/Dkn + 1/Dox);

    // Intrinsic reactivity [1/s]
    const scalar ki = Ai_*exp(-(Ei_*1000)/(RR/1000.0)/T);

    // Thiele modulus []
    const scalar phi =
      max(0.5*d*sqrt(Sb_*rhop*(Ag_*1000.)*ki*pc/(De*rhoc)), ROOTVSMALL);

    // Effectiveness factor []
    const scalar eta = max(3.0/sqr(phi)*(phi/tanh(phi) - 1.0), 0.0);

    // Chemical rate [kmol/m2/s]
    const scalar Rch = eta*d/6.0*rhop*(Ag_*1000.)*ki;

    // Particle surface area [m2]
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Change in C mass [kg]
    scalar dmC = Ap*pc*YO2*Rox*Rch/(Rox + Rch)*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass*Ychar, dmC);

    // Molar consumption [kmol]
    const scalar dOmega = dmC/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega*Sb_*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = dOmega*(WC_ + Sb_*WO2_);

    // Update local particle C mass
    dMassSolid[CsLocalId_] += dOmega*WC_;

    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;

    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    return dmC*HsC - dmCO2*HcCO2_;
}


// ************************************************************************* //
