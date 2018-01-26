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
    Foam::PcDevolatilisationModel

Description
    Templated devolatilisation model class

SourceFiles
    PcDevolatilisationModel.C
    PcDevolatilisationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PcDevolatilisationModel_H
#define PcDevolatilisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PcDevolatilisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PcDevolatilisationModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //- Runtime type information
    TypeName("devolatilisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PcDevolatilisationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PcDevolatilisationModel(CloudType& owner);

        //- Construct from dictionary
        PcDevolatilisationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        PcDevolatilisationModel(const PcDevolatilisationModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<PcDevolatilisationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~PcDevolatilisationModel();


    //- Selector
    static autoPtr<PcDevolatilisationModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


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
        ) const = 0;

        //- Add to devolatilisation mass
        void addToDevolatilisationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePcDevolatilisationModel(CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::pcReactingMultiphaseCloudType                   \
        pcReactingMultiphaseCloudType;                                       \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::PcDevolatilisationModel<pcReactingMultiphaseCloudType>,      \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            PcDevolatilisationModel<pcReactingMultiphaseCloudType>,        \
            dictionary                                                         \
        );                                                                     \
    }


#define makePcDevolatilisationModelType(SS, CloudType)                       \
                                                                               \
    typedef Foam::CloudType::pcReactingMultiphaseCloudType                   \
        pcReactingMultiphaseCloudType;                                       \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<pcReactingMultiphaseCloudType>, 0);                        \
                                                                               \
    Foam::PcDevolatilisationModel<pcReactingMultiphaseCloudType>::         \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<pcReactingMultiphaseCloudType>>                            \
        add##SS##CloudType##pcreactingMultiphaseCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PcDevolatilisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
