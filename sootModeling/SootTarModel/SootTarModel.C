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

Description

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SootTarModel.H"

// Constructor
Foam::SootTarModel::SootTarModel
(
    const psiReactionThermo& thermo, 
    const basicSpecieMixture& composition,
    const volScalarField& Ns,
    const fvMesh& mesh
)
    :
    thermo_(thermo),
    composition_(composition),
    mesh_(mesh),
    cellState_
    (
        thermo,
        composition,
        Ns,
        mesh,
        wordList({"TAR", "SOOT", "CO", "CO2", "O2", "H2"}),
        wordList()
    ),
    MW_(6),
    speciesY_(6),
    speciesSources_(6),
    Qdot_(Ns.size(),0.0)
{
    Info<< "Creating Tar Breakdown  Model Solver \n\n" << endl;

    // // Get the dictionary
    // IOdictionary coalPropertiesDict
    // (
    //     IOobject
    //     (
    //         "coalCloud1Properties",    // dictionary name
    //         mesh.time().constant(),     // dict is found in "constant"
    //         mesh,                   // registry for the dict
    //         IOobject::MUST_READ,    // must exist, otherwise failure
    //         IOobject::NO_WRITE      // dict is only read by the solver
    //     )
    // );

    // // Get the actual tar breakdown data
    // dictionary tarBreakdownDict(coalPropertiesDict.subDict("subModels").subDict("tarBreakdown"));

    // ITstream tokenStream(tarBreakdownDict.lookup("species"));

    // label idx =tokenStream.tokenIndex();
    
    // Info <<  idx 
    //     << endl;

    // char temp(tokenStream.readBeginList("List"));

    // idx =tokenStream.tokenIndex();
    
    // Info <<  idx 
    //     << endl;

    // word tempNum(tokenStream);

    // idx =tokenStream.tokenIndex();

    // Info <<  idx
    //     << endl;

    // FatalErrorInFunction << "END" << abort(FatalError);
    //relevantSpecies_ = temp;
    
    // list of relevant species names
    const wordList relevantSpecies_({"TAR", "SOOT", "CO", "CO2", "O2", "H2"});

    // Initialize all hash table entries
    forAll(relevantSpecies_, s)
    {
        const word name =relevantSpecies_[s];

        if (! speciesSources_.found(name))
        {
            speciesSources_.insert(name, scalarField(Ns.size(),0.0));
            MW_.insert(name, composition.W(composition.species()[name]));
            speciesY_.insert(name, threeList(0.0));
        }
        else // Pretty much deprecated, don't worry about it.
        {
            FatalErrorInFunction << "Species " << name
                << " listed twice in tarBreakdown species dictionary." 
                << relevantSpecies_ << abort(FatalError);
        }
    }
}

// Member functions
void Foam::SootTarModel::ratesOfChange
(
    scalar& r_oxidation,
    scalar& r_gasification,
    scalar& r_sootFormation
)
{
    // Get the necessary cell reactant value from CellState
    // I realize this is redundant it just seems safer this way
    const scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    const scalar cell_T = this->cellState_.thermoProperties()["T"];
    const scalar cell_tar = this->cellState_.frozenSpecieMassFractions()["TAR"];
    const scalar cell_o2 = this->cellState_.frozenSpecieMassFractions()["O2"];
    // Universal R in [kJ/(mol K)]
    const scalar R(8.314e-3);

    // All from Fletcher Revision paper 1998
    const scalar A_ox(6.77e5);   //[m^3/(kg s)]
    const scalar A_gas(9.77e10); //[1/s]
    const scalar A_soot(5.02e8);   //[1/s]
    // Activation Energies all in [kJ/mol]
    const scalar E_ox(52.3);
    const scalar E_gas(286.9);
    const scalar E_soot(198.9);

    // Rates based on current cell state
    // modified from [kg_tar /(m^3*s)] to [kmol/(m^3*s)] 
    // by division with the molar mass of tar
    r_oxidation = cell_rho * (cell_tar * cell_o2) * A_ox * 
        Foam::exp(-E_ox/(R*cell_T)) / MW_["TAR"];
    r_gasification = cell_rho * cell_tar * A_gas * 
        Foam::exp(-E_gas/(R*cell_T)) / MW_["TAR"];
    r_sootFormation = cell_rho * cell_tar * A_soot *
        Foam::exp(-E_soot/(R*cell_T)) / MW_["TAR"];

} //end ratesOfChange

void Foam::SootTarModel::explicitStep
(
    const label cellNumber,
    const scalar subdt
)
{
    // density in this cell
    const scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    
    // rates for soot mass fraction reactions [kmol/(m^3 * s)]
    scalar r_oxidation(0.0);
    scalar r_gasification(0.0);
    scalar r_sootFormation(0.0);
    
    // calculate the rates
    ratesOfChange(r_oxidation, r_gasification, r_sootFormation);

    HashTable<threeList>& s(speciesY_);

    // factor to change rate from [kg/m^3s] (once multilpied with MW)
    // to an increment in mass fraction => [-] no unit
    const scalar factor(subdt / cell_rho);

    // There are three reactions in this model and TAR is assumed to be Pyrene (C16H10)
    // 1. Soot Formation
    //     1 TAR -> 29.3 SOOT (where SOOT == C, here only mass is conserved)
    // 2. Combustion
    //    1 TAR -> 16 CO2 + 5 H2 (here both mass and elements are conserved)
    // 3. Cracking
    //    1 TAR -> 7.065 CO + 2.208 H2 (here both mass and C/H ratio are conserved)
    
    // Reactants
    s["TAR"][2] = s["TAR"][1] - 
        MW_["TAR"] * (r_oxidation + r_gasification + r_sootFormation) * factor;
    s["O2"][2] = s["O2"][1] - MW_["O2"] * 16.0 * r_oxidation * factor;
    // Products
    s["SOOT"][2] = s["SOOT"][1] + MW_["SOOT"] * 29.3 * r_sootFormation * factor;
    s["CO"][2] = s["CO"][1] + MW_["CO"] * 7.065 * r_gasification * factor;
    s["CO2"][2] = s["CO2"][1] + MW_["CO2"] * 16.0 * r_oxidation * factor;
    s["H2"][2] = s["H2"][1] + MW_["H2"] * 
        (2.208 * r_gasification + 5 * r_oxidation) * factor;

}// end explicitStep

void Foam::SootTarModel::correctQdot()
{

    // Reset from last time
    this->Qdot_ = 0.0;

    // loop through the species sources that have
    // been updated

    forAllIter(HashTable<scalarField>, speciesSources_, iter)
    {
        word specieName = iter.key();
        label specieIndex = composition_.species()[specieName];
        scalarField& sourceField = iter(); //[kg/(m^3 * s)]

        // get the enthalpy of formation for this specie
        scalar hF = composition_.Hc(specieIndex);
        
        Qdot_ -= hF * sourceField;
    }    
}

Switch Foam::SootTarModel::negativeFinalMassFraction()
{
    // Would be best to make this just check the reactants probably.
    // The products shouldn't be a problem in this regard
    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        // species name
        const word name(iter.key());
        // species mass fraction threeList
        const threeList& mf(iter());

        if (mf[2] < 0.0) // if the final mf is less than 0.
        {
            return true;
        }
    }

    return false;
}

void Foam::SootTarModel::exhaustWorstSpecies
(
    const scalar& subdt,
    scalar& smallTimeStep
)
{

    // Get the species that goes to zero first (in time sense)
    const word name = findWorstSpecies();
    
    // Interpolate species mf to zero and apply that to the subdt
    // to find smallTimeStep. Assumes positive current value [1]
    // and negative final value [2]
    const scalar fraction = speciesY_[name][1] / 
        (speciesY_[name][1] - speciesY_[name][2]);

    smallTimeStep = subdt * fraction;    
    
    // now reset the final species values to only be integrated
    // through smallTimeStep instead of subdt as they are now
    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        // get the mf list for this species
        threeList& mf = iter();

        // Apply interpolation to the species
        mf[2] = mf[1] + fraction * (mf[2] - mf[1]);
    }
     
    if 
    (
        (Foam::cmptMag(speciesY_[name][2]) >= Foam::SMALL)
    )
    {
        FatalErrorInFunction << "Expected " 
            <<  name << " mass fraction of very near 0.0 but found it to be " 
            << speciesY_[name][2] << "." << abort(FatalError);
    }

    // As long as the interpolation got it near zero, just set it to 
    // zero here
    speciesY_[name][2] = 0.0;
}

word Foam::SootTarModel::findWorstSpecies()
{
    
    
    // Assume that there are more than one negative species
    DynamicList<word> negativeSpecies;

    // I realize that we looped through already in negativeFinalMassFraction()
    // Maybe they should be combined in the future.
    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        // species name
        const word name(iter.key());
        // species mass fraction threeList
        const threeList& mf(iter());

        if (mf[2] < 0.0)
        {
            negativeSpecies.append(name);
        }
    }
    
    if (negativeSpecies.size() == 1)
    {
        return negativeSpecies[0];
    }
    
    else // Now we need to find the worst one
    {

        // The ratio is question is Y_current/(Y_current - Y_final)
        // where Y_current is assumed to be positve and Y_final negative.
        // Of the species that have gone negative 
        // that which has the lowest ratio will be the one 
        // that goes negative first and is the name that should be returned.
        scalar smallRatio(1.0);
        word worstSpeciesName;

        forAll(negativeSpecies, s)
        {
            word name = negativeSpecies[s];
            
            scalar ratio = speciesY_[name][1]/ 
                (speciesY_[name][1] - speciesY_[name][2]);

            if (ratio < smallRatio)
            {
                smallRatio = ratio;
                worstSpeciesName = name;
            }
        } // end loop through negativeSpecies

        return worstSpeciesName;
    }
}

void Foam::SootTarModel::updateSpeciesMassFractions()
{

    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        threeList& mf = iter();
        mf[1] = mf[2];
    }
    
}

void Foam::SootTarModel::setSpeciesSources
(
    const scalar& dt
)
{
    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        const word name = iter.key();
        const threeList mf = iter();
        const label cellNumber = cellState_.cellNumber();
        const scalar cellRho = cellState_.thermoProperties()["rho"];

        // Note: density multiplication because the transport
        // equations are d(rho Y)/dt and rho is constant
        // in this context.
        speciesSources_[name][cellNumber] = (mf[2] - mf[0]) * (cellRho/dt);
    }
}

void Foam::SootTarModel::calcSpecieSources
(
    const scalar dt,
    const label nSubSteps,
    const label cellNumber
)
{

    // Local storage for mass fractions of species and soot variables
    scalarField Y_current(10,0.0);
    // Create another field for storage of Mass Fractions/Ns
    // at the end of the sub time step
    scalarField Y_final(Y_current.size(), 0.0);

    // set default initial sub time step for this cell
    scalar subdt = dt/nSubSteps;

    // Grab the current main field data for this cell
    // this will update the cellState thermo data too
    this->cellState_.updateCellState(cellNumber);

    // Update the local speciesY_ mass fractions from CellState
    // this pulls the state from the main volScalarFields
    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        const word name = iter.key();
        threeList& massFractions = iter();
        
        scalar cellMassFraction = cellState_.frozenSpecieMassFractions()[name];

        // Set the initial, current and final mass fractions
        massFractions[0] = cellMassFraction; // initial
        massFractions[1] = cellMassFraction; // current
        massFractions[2] = cellMassFraction; // final
    }

    // Integrate
    scalar integratedTime(0.0);
    scalar smallTimeStep(0.0);
    while(integratedTime < dt)
    {
        // Advance, explicitly, local cell copy of species and Ns
        this->explicitStep(cellNumber, subdt);

        // It's possible that the explicit step may be too long 
        // and we exhaust the supply of a reactant specie
        if (negativeFinalMassFraction())
        {
            // Integrates only until the first specie to reach zero does 
            // so. Then calculates how much time that would take, that is
            // 'smallTimeStep'. It resets the speciesY_ final fields appropriately
            // as integrated only for smallTimeStep, and not subdt.
            exhaustWorstSpecies(subdt, smallTimeStep);
                
            integratedTime += smallTimeStep;
            subdt = subdt - smallTimeStep;
        }
        else
        {
            // If there was not problem with negative mass fractions
            // record the time step time taken
            // and then make sure subdt is set back to the default.
            integratedTime += subdt;
            subdt = dt/nSubSteps;
        }

        // set current mass fractions to final
        updateSpeciesMassFractions();
        
        // Update cellState mass fractions
        cellState_.updateCellState(speciesY_);

    } // end while (integratedTime < dt)

    // If mass fractions are negative here then there is a 
    // problem.
    if (negativeFinalMassFraction()) 
    {
        FatalErrorInFunction 
            << "Negative species mass fraction detected."
                << abort(FatalError);
    }

    setSpeciesSources(dt);

}// End calcSpecieSources


//****************** Public member functions ********************//


Foam::tmp<Foam::volScalarField> Foam::SootTarModel::Qdot()
{
        tmp<volScalarField> tQdot
        (
            new volScalarField
            (
                IOobject
                (
                    "sootQdot",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
            )
        );

        tQdot.ref().primitiveFieldRef() = this->Qdot_;

        return tQdot;
}

Foam::tmp<Foam::volScalarField> Foam::SootTarModel::sourceY
(
    const volScalarField& Y_field
)
{
    word specieName = Y_field.name();

    tmp<volScalarField> tSourceY
        (
            new volScalarField
            (
                IOobject
                (
                    "sourceY",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );

    // get the reference to volScalarField that tmp is wrapped around
    volScalarField& sourceY = tSourceY.ref();

    if (this->speciesSources_.found(specieName))
    {
        sourceY.primitiveFieldRef() = this->speciesSources_[specieName];
        
        return tSourceY;
    }
    else
    {
        // return empty source because this specie is unrelated to this model
        return tSourceY;
    }
}
 

// Main function for calculating the sources.
// Used in transient non-LTS solvers with same dt for 
// all cells
void Foam::SootTarModel::updateSources
(
    const label nSubSteps
)
{  
    // Update the fields stored in cellState to new time
    // Currently stores mixture MW field and gas density field.
    this->cellState_.updateCellStateFields(); 
    
    // Calculate and set the species sources into this->speciesSources_ table.
    // If statement to decide if LTS or transient.
    // In either case get the time step/steps and iterate through cells
    if (fv::localEulerDdt::enabled(this->mesh_))
    {
        Info << "Operating tar model in LTS mode" << endl;
        // get the LTS timsteps inverses (if there is a function
        // that just gets the time steps then use that
        const scalarField& localDtInv = 
            fv::localEulerDdt::localRDeltaT(this->mesh_);
        forAll(this->mesh_.C(), cellNumber)
        {
            this->calcSpecieSources(1.0/localDtInv[cellNumber], nSubSteps, 
             cellNumber);
        }
    }
    // if single time step (i.e. transient simulation)
    else if (! fv::localEulerDdt::enabled(this->mesh_))
    {
        Info << "Operating tar model in transient mode" << endl;
        // get global time step
        const scalar& dt = this->mesh_.time().deltaTValue();
        forAll(this->mesh_.C(), cellNumber)
        {
            this->calcSpecieSources(dt, nSubSteps, 
            cellNumber);
        }
    }
    else
    {
        // I'm not really sure what the localEulerDdt::enabled() is capable of 
        // returning.
        FatalErrorInFunction << "Not LTS or Transient" << abort(FatalError); 
    }

    
    // based on the updated species sources calculate the enthalpy source
    this->correctQdot();
}


// **************** Access Functions ************************//
