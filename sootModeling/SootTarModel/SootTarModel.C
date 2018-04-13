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
    cellState_(thermo, composition, Ns, mesh),
    MW_(),
    speciesY_(),
    relevantSpecies_({"TAR", "SOOT", "CO", "CO2", "O2", "H2"}),
    speciesSources_(),
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

    // append CO, O2, H2 and TAR for combustion if not already present
    // List<Switch> combSpecies(4, false); // [CO present, O2, H2, TAR]
    
    // forAll(relevantSpecies_, speciesIdx)
    // {
    //     const word name = relevantSpecies_[speciesIdx];

    //     if (name == "CO")
    //     { 
    //         combSpecies[0] = true;
    //     }
    //     if (name == "O2") 
    //     {
    //         combSpecies[1] = true;
    //     }
    //     if (name == "H2") 
    //     {
    //         combSpecies[2] = true;
    //     }
    //     if (name == "TAR") 
    //     {
    //         combSpecies[3] = true;
    //     }
    // }
    
    // forAll(combSpecies, s)
    // {
    //     if (!combSpecies[s] && s==0) 
    //     {
    //         relevantSpecies_.append("CO");
    //     }
    //     if (!combSpecies[s] && s==1) 
    //     {
    //         relevantSpecies_.append("O2");
    //     }
    //     if (!combSpecies[s] && s==2) 
    //     {
    //         relevantSpecies_.append("H2");
    //     }
    //     if (!combSpecies[s] && s==3) 
    //     {
    //         relevantSpecies_.append("TAR");
    //     }
    // }

    const label speciesNumber(relevantSpecies_.size());

    // make sure the default tables are big enough
    if (speciesNumber > MW_.capacity())
    {
        MW_.resize(speciesNumber);
        speciesY_.resize(speciesNumber);
        speciesSources_.resize(speciesNumber);
        
    }

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
        else 
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
   
    // ############################################################### ///
    // // get relevant state variables for the cell
    // scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    // scalar cell_T = this->cellState_.thermoProperties()["T"];
    // scalar cell_c2h2 = this->cellState_.frozenSpecieMassFractions()["C2H2"];
    // scalar cell_Ys = this->cellState_.Ysoot();
    // scalar cell_Ns = this->cellState_.Nsoot();
    // // Partial pressures
    // scalar P_co2 = this->cellState_.frozenSpeciePPressures()["CO2"];
    // scalar P_o2 = this->cellState_.frozenSpeciePPressures()["O2"];
    // scalar P_h2o = this->cellState_.frozenSpeciePPressures()["H2O"];
    // scalar P_oh = this->cellState_.frozenSpeciePPressures()["OH"];
    // // Specie concentrations
    // const scalar C_c2h2(cell_rho * cell_c2h2 * (1/this->MW_["C2H2"]));
    // // Molar Mass of soot
    // const scalar M_s = this->MW_["SOOT"];

    // // Soot particle diameter and surface area for cell
    // // (these are considered uniform in the cell)
    // scalar cell_dp(0.0);
    // scalar cell_As(0.0);
    // if (cell_Ns > Foam::SMALL)
    // {
    //     cell_dp = 
    //         Foam::pow
    //         (
    //             (6.*cell_Ys)/(this->pi * this->rho_s * cell_Ns),
    //             (1./3.)
    //         );

    //     cell_As = this->pi * Foam::sqr(cell_dp) *
    //         cell_rho * cell_Ns;
    // }
      
    // // Set constants from papers
    // // Kronenburg(2000)
    // scalar A_growth(7.5e2);
    // scalar T_growth(12100.0);
    // // Josephson (2017)
    // scalar A_o2(7.98e-1);
    // scalar A_oh(1.89e-3);
    // scalar A_co2(3.06e-17);
    // scalar A_h2o(6.27e4);
    // scalar E_o2(1.77e5);
    // scalar E_co2(5.56e3);
    // scalar E_h20(2.95e5);
    // scalar n_(0.13);
    // scalar R_(8.314); //[kJ/(kmol*K)]

    // // Now calculate the rates
    // // Soot Mass fraction [kmol/(m^3*s)]
    // r_growth = A_growth * Foam::exp(-T_growth/cell_T) * cell_As * C_c2h2;
    // r_oxidation_O2 = Foam::pow(cell_T, -0.5) * 
    //     A_o2 * P_o2 * Foam::exp(-E_o2/(R_*cell_T)) *
    //     (cell_As/M_s);
    // r_oxidation_OH = Foam::pow(cell_T, -0.5) * A_oh * P_oh * (cell_As/M_s);
        
    // r_gasification_CO2 = A_co2 * Foam::pow(P_co2,0.5) * Foam::sqr(cell_T) *
    //     Foam::exp(-E_co2/(R_*cell_T)) * (cell_As/M_s);
       
    // r_gasification_H2O = A_h2o * Foam::pow(P_h2o,n_) * Foam::pow(cell_T, -0.5) *
    //     Foam::exp(-E_h20/(R_ * cell_T)) * (cell_As/M_s);
    
    // r_agglomeration = 2 * this->Ca * Foam::sqrt(cell_dp) * 
    //     Foam::sqrt((6.0 * this->sigma * cell_T)/(this->rho_s)) * 
    //     Foam::sqr(cell_rho * cell_Ns);

    // scalar cell_oh = this->cellState_.frozenSpecieMassFractions()["OH"];
    // const scalar C_oh(cell_rho * cell_oh * (1/this->MW_["OH"]));

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

    // reference to make it easier to type this
    HashTable<threeList>& s(speciesY_);

    // factor to change rate from [kg/m^3s] (once multilpied with MW)
    // to an increment in mass fraction => [-] no unit
    const scalar factor(subdt / cell_rho);

    // Take explicit step
    s["TAR"][2] = s["TAR"][1] - 
        MW_["TAR"] * (r_oxidation + r_gasification + r_sootFormation) * factor;
    // Soot formation
    s["SOOT"][2] = s["SOOT"][1] + MW_["SOOT"] * r_sootFormation * factor;
    // Tar cracking/gasification
    s["CO"][2] = s["CO"][1] + MW_["CO"] * r_gasification * factor;
    // Tar combustion
    s["CO2"][2] = s["CO2"][1] + MW_["CO2"] * r_oxidation * factor;
    s["O2"][2] = s["O2"][1] - MW_["O2"] * r_oxidation * factor;
    // Both Tar combustion and gasification
    s["H2"][2] = s["H2"][1] + MW_["H2"] * (r_gasification + r_oxidation) * factor;

}// end explicitStep

void Foam::SootTarModel::advanceToZero
(
    const scalar totalTime,
    scalar& smallSubTime,
    const scalarField& Y_initial,
    scalarField& Y_current
)
{
    // change in each field over totalTime
    scalarField Y_change = Y_initial - Y_current;

    
    label worstIndex(-1);
    scalar worstCurrent(0.0);
    // Find the most negative mass fraction and it's index
    forAll(Y_current, specieIndex)
    {
        if (Y_current[specieIndex] < worstCurrent)
        {
            worstCurrent = Y_current[specieIndex];
            worstIndex = specieIndex;
        }
    }
    
    if (worstIndex == -1)
      {
	FatalErrorInFunction << "No negative mass fractions." << abort(FatalError);
      }

    // grab the value of the worst specie and the beginning of the cfd ts
    scalar worstInitial = Y_initial[worstIndex];

    // Interpolate the problem specie to zero. Apply that interpolation
    // to the evolution of the other species as well.
     
    // assume it's decreasing and worstInitial is positive. Find the total change
    scalar totalChange =  worstInitial - worstCurrent;
    scalar fraction = worstInitial/totalChange;

    // Now apply that interpolation to the entire system
    Y_current = Y_initial -  Y_change*fraction;
    smallSubTime= totalTime*fraction;

    if (smallSubTime < 0.0)
      {
	FatalErrorInFunction << "Negative Time Step" << abort(FatalError);
      }

    // Just make sure that the problem specie is pretty much zero
    Y_current[worstIndex] = 0.0;//Foam::VSMALL;
}

void Foam::SootTarModel::exhaustLowSpecie
(
    const scalar totalTime, 
    scalar& smallSubTime,
    const scalarField& Y_initial,
    scalarField& Y_final
)
{
    smallSubTime = 0.0;
    // TODO: this whole function can be done better with a rewrite of 
    // advanceToZero(). It just fixes the most negative value. That leaves
    // open the possibility that another negative number wont be fixed 
    // depending on the initial values. What should be done is to fix
    // the specie that goes to zero first, not the one that goes the
    // most negative. They don't have the same rates.
    // For now just do an initial correction and then correct
    // more over the new timestep if there are still negative values
    // after the first advancetoZero call.
    this->advanceToZero(totalTime, smallSubTime, Y_initial, Y_final);

    // check to see that we got it otherwise take care of other
    // negative mass fractions
    bool allPositive = (min(Y_final) >= 0.0);
    scalar newDt(0.0);
    while (!allPositive)
    {	
        this->advanceToZero(smallSubTime, newDt, Y_initial, Y_final);
        smallSubTime = newDt;
        allPositive = (min(Y_final) >= 0.0);
    }
}

void Foam::SootTarModel::correctQdot()
{

    // Reset from last time
    this->Qdot_ = 0.0;

    // loop through the species sources that have
    // been updated

    forAllIter(HashTable<scalarField>, this->speciesSources_, iter)
    {
        word specieName = iter.key();
        label specieIndex = this->composition_.species()[specieName];
        scalarField& sourceField = iter(); //[kg/(m^3 * s)]

        // get the enthalpy of formation for this specie
        scalar hF = this->composition_.Hc(specieIndex);
        
        this->Qdot_ -= hF * sourceField;
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
    scalar subdt, 
    scalar smallTimeStep
)
{

    // Get the species that goes to zero first (in time sense)
    const word name = findWorstSpecies();
    
    // Interpolate species mf to zero and apply that to the subdt
    // to find smallTimeStep. Assumes positive current value [1]
    // and negative final value [2]
    scalar fraction = speciesY_[name][1] / 
        (speciesY_[name][1] - speciesY_[name][2]);

    smallTimeStep = subdt * fraction;    
    
    // now reset the final species values to only be integrated
    // through smallTimeStep instead of subdt as they are now
    forAllIter(HashTable<threeList>, speciesY_, iter)
    {
        // get the mf list for this species
        threeList mf = iter();

        // Apply interpolation to the species
        mf[2] = mf[1] + fraction * (mf[2] - mf[1]);
    }
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
                
            // Record the small time integration
            integratedTime += smallTimeStep;
            // Set subdt for the next time step
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

        // Update both the soot variables and the frozen species for
        // the next iteration
        Y_current = Y_final;

        // before the rates are recalculated on the next iteration update the
        // cell state but only frozen species and soot mass fractions,
        // and soot particle number density. 
        //The thermo/cell volume shouldn't be changing anyway.
        cellState_.updateCellState(Y_current, Y_current[8], Y_current[9]);

    } // end while (integratedTime < dt)

    // density in this cell (assumed constant over these steps)
    const scalar cellRho = this->cellState_.thermoProperties()["rho"];

    if  // There is mass fraction exceeding one
        (
            Y_final[0] > 1.0 ||
            Y_final[1] > 1.0 ||
            Y_final[2] > 1.0 ||
            Y_final[3] > 1.0 ||
            Y_final[4] > 1.0 ||
            Y_final[5] > 1.0 ||
            Y_final[6] > 1.0 ||
            Y_final[7] > 1.0 ||
            Y_final[8] > 1.0
        )
    {
        Info << "WARNING: Mass Fraction Greater than one" << endl;
        Info << "After sources:\n"  << endl;
        Info << "\nCell number: " << cellNumber << endl;
        // Info << "Y_initial: \n " << Y_initial << endl;
        // Info << "Y_final: \n " << Y_final << "\n\n" <<  endl;
    }

    if (min(Y_final) < 0.0)
    {
        FatalErrorInFunction 
            << "Negative mass fraction after soot model integration"
                << abort(FatalError);
    }

    // Now set the source terms given the changes in soot mass fraction,
    // frozen specie mass fraction  and soot number density. These 
    // source terms will be used when the respective transport
    // equations are solved
    // NOTE: Multiply by density because the transport
    //  equations are d(rho*Y)/dt not just d(Y)/dt, 
    //  we take rho to be constant here so just pull it out.

    // Soot number density
    // this->N_source[cellNumber] = 
    //     (Y_final[9] - Y_initial[9])*cellRho/dt;
    // Soot mass fraction
    // this->speciesSources_["SOOT"][cellNumber] = 
    //     (Y_final[8] - Y_initial[8])*cellRho/dt;

    // //  Other species
    // NOTE: Multiply by density because the transport
    //  equations are d(rho*Y)/dt not just d(Y)/dt
    // and rho is constant over this step

    // // Reactant species 
    // this->speciesSources_["C2H2"][cellNumber] = 
    //     (Y_final[0] - Y_initial[0])*cellRho/dt;
    // this->speciesSources_["O2"][cellNumber] = 
    //     (Y_final[1] - Y_initial[1])*cellRho/dt;
    // this->speciesSources_["OH"][cellNumber] = 
    //     (Y_final[2] - Y_initial[2])*cellRho/dt;
    // this->speciesSources_["H2O"][cellNumber] = 
    //     (Y_final[7] - Y_initial[7])*cellRho/dt;
    // // Product species
    // this->speciesSources_["H2"][cellNumber] = 
    //     (Y_final[3] - Y_initial[3])*cellRho/dt;
    // this->speciesSources_["CO"][cellNumber] = 
    //     (Y_final[4] - Y_initial[4])*cellRho/dt;
    // this->speciesSources_["H"][cellNumber] = 
    //     (Y_final[5] - Y_initial[5])*cellRho/dt;
    // // Both a product and a reactant
    // this->speciesSources_["CO2"][cellNumber] = 
    //     (Y_final[6] - Y_initial[6])*cellRho/dt;

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

Foam::tmp<Foam::volScalarField> Foam::SootTarModel::sourceN()
{
    //tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(N_soot, N_source_dims));

    tmp<volScalarField> tSourceN
        (
            new volScalarField
            (
                IOobject
                (
                    "sourceN",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", dimless/dimVolume/dimTime, 0.0)
            )
        );

    // // get reference to the fvmatrix inside the tmp
    // volScalarField& Nsource = tSourceN.ref();

    // Nsource.primitiveFieldRef() = this->N_source;

    // // reset field to zero in preparation for next time step.
    // this->N_source = 0.0;

    return tSourceN;
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
        
        
        // reset field to zero in preparation for next time step.
        this->speciesSources_[specieName] = 0.0;
        return tSourceY;
    }
    else
    {
        // return empty source because this specie is unrelated to the soot model
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
        Info << "Operating soot model in LTS mode" << endl;
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
        Info << "Operating soot model in transient mode" << endl;
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
