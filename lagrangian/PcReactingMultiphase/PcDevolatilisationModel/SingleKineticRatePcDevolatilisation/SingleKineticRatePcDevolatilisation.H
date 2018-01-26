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
    Foam::SingleKineticRatePcDevolatilisation

Description
    Single kinetic rate devolatisation model for the 
    primary particle pyrolysis (i.e. devolatilization) and
    secondary pyrolysis (i.e. the breakdown of tar into
    smaller molecules that are actually in the chemical 
    mechanism).

    Primary Pyrolysis (PP):
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A*exp(- E/R*T)

    Where:
        kappa = rate constant
        Ap    = pre-exponential factor (user input) [1/s]
        Ep     = activation energy (user input) [kcal/mole]
        R     = universal gas constant
        T     = temperature


    Secondary Pyrolysis (SP):
    - Through PP TAR is treated as a normal species (therefore
      there needs to be an entry named 'TAR' in the chemical mechanism).
      In PP the tar is evolved and returned appropriately in dMassDV, still
      as if it is going to be released into the gas phase.
    - In SP (within this class) we calculate, given the mass of tar being
      evolved, the mass of the actual gas species that should be released
      in its place and store them in dMassSP (which is ordered just 
      as dMassDV so dMassDV + dMassSP makes sense). dMassSP is returned by 
      the main calculate function in this class as well.
    - In the outer function (calc() in PcReactingMultiphaseParcel.C) the tar
      in dMassDV is left untill the parcel mass fractions have updated =>
      the tar fraction of the particle is update appropriately. BUT ! the tar 
      is then removed and replaced with its breakdown products stored in dMassSP
      before the source terms to the gas phase are calculated => No TAR is actually
      release to the gas phase.

      


    Usage:

        SingleKineticRateDevolatilisationCoeffs
        {
            volatileData
            (
	        (TAR     382    7.88   3.41e8   44    .24)
                (CH4     12     0.5    0.0      0.0   .23) 
                (CO2     12     0.5    0.0      0.0   .35)
		(SOOT    0      0      0        0     .42)
            );

            volatileResidualCoeff 1e-6;
        }

	The parameters are passed in in this form
	(speciesName    Ap    Ep    As    Es    YInfs)
	
	The forms of the entry are dependant on whether the specie is
	TAR or not:

	For non-TAR:
	 -Ap and Ep are the standard primary pyrolysis rates, they have nothing 
	  to do with the breakdown of tar. From PCCL V. 5.0 this information
	  is obtained in the file 'SFORC1T1.rpt'. Depending on whether it is a
	  hydrocabon or some other PCCL specie class. Enter here in units 
          provided in PCCL i.e. A=[1/s] and E=[kcal/mol]

	 -In these species Yinfs stands for the 
	  mass fraction of the species from tar breakdown. In other words
	  for species X with Yinfs = 0.5 there will be 0.5 grams of X produced
	  when 1 gram of Tar is broken down (through secondary pyrolysis). This
	  fraction information is obtained from the PCCL file 'SFTRC1T1.rpt'.
	  Note that excluding tar sum(Yinfs) = 1.0, i.e. you need to fully account
	  for where the tar mass will go.

	 -As and Es should be set to 0.

	For TAR:
	 -Ap and Ep are just the normal speciated PP rate parameters. In this 
	  sense TAR is treated as any other volatile specie. These parameters
	  are also obtained from the 'SFORC1T1.rpt' file (taking the rates for 
	  tar in that file).

	 -Yinfs in the tar line is the daf mass fraction of tar initially in the 
	  particle. This is also the ultimate daf mass fraction  yield of Tar.
	  Also can be obtained from 'SFORC1T1.rpt' file. Note that it has no 
	  relationship to the Yinfs from the other species.

	 -As and Es are the rate parameters for the Arrhenius reation that
	  describes the breakdown of tar into its constituents through SP.
	  This corresponds to Eqn. 9 from page 13-5 in the PCCL manual (i.e.
	  the 'dY' equation).

	

	 /// ********************************************************
	IMPORTANT BUG NOTE:
	 Right now there is one known problem in the code, TAR must be the 
	 first specie listed in the cloudProperties when specifying the gas
	 mass fractions within the singleMixtureFractionCoeffs dict.
	 /// ********************************************************
	 
	

\*---------------------------------------------------------------------------*/

#ifndef SingleKineticRatePcDevolatilisation_H
#define SingleKineticRatePcDevolatilisation_H

#include "PcDevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class SingleKineticRatePcDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SingleKineticRatePcDevolatilisation
:
    public PcDevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
        // Private data

            //- Specie name
            word name_;

            //- Activation energy Primary
            scalar Ap_;

            //- Pre-exponential factor Primary
            scalar Ep_;

            //- Activation energy Secondary
            scalar As_;

            //- Pre-exponential factor Secondary
            scalar Es_;

            //- Ultimate species mass fraction (within tar)
            scalar Yinfs_;

    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
                Ap_(0.0),
                Ep_(0.0),
                As_(0.0),
                Es_(0.0),
		Yinfs_(0.0)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
                Ap_(readScalar(is)),
                Ep_(readScalar(is)),
                As_(readScalar(is)),
                Es_(readScalar(is)),
                Yinfs_(readScalar(is))
            {}

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
                Ap_(vd.Ap_),
                Ep_(vd.Ep_),
                As_(vd.As_),
                Es_(vd.Es_),
                Yinfs_(vd.Yinfs_)
            {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the activation energy
                scalar Ap() const
                {
                    return Ap_;
                }

                //- Return const access to the pre-exponential factor
                scalar Ep() const
                {
                    return Ep_;
                }
                //- Return const access to the activation energy
                scalar As() const
                {
                    return As_;
                }

                //- Return const access to the pre-exponential factor
                scalar Es() const
                {
                    return Es_;
                }
                //- What do you suppose this function does?
                scalar Yinfs() const
                {
                    return Yinfs_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.Ap_  >> vd.Ep_ >> 
		  vd.As_ >> vd.Es_ >> vd.Yinfs_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.Ap_ << token::SPACE
                    << vd.Ep_ << token::SPACE
                    << vd.As_ << token::SPACE
                    << vd.Es_ << token::SPACE
		    << vd.Yinfs_
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants

            //- List of volatile data - (name A1 E)
            List<volatileData> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;

  // In order to work with daf basis in main funtion
  // we need the mass fraction of ash and water of
  // the particles so Ydaf0 is the initial 
  // mass fraction of the particle that is daf
  // it is assumed that these fractions are the same for all particles
  scalar Ydaf0_;


public:

    //- Runtime type information
    TypeName("singleKineticRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        SingleKineticRatePcDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        SingleKineticRatePcDevolatilisation
        (
            const SingleKineticRatePcDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<PcDevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<PcDevolatilisationModel<CloudType>>
            (
                new SingleKineticRatePcDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SingleKineticRatePcDevolatilisation();


    // Member Functions

        //- Update model
        virtual void calculate
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
	    scalarField& tarProps,
	    scalarField& dMassSP
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SingleKineticRatePcDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
