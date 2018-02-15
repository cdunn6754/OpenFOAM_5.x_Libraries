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

#include "EulerImplicitSystem.H"


// Constructor
Foam::EulerImplicitSystem::EulerImplicitSystem
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
    Y_source_dims(1,0,-1,0,0),// Source dimensions (already multipied by volume)
    N_source_dims(0,0,-1,0,0),// Source dimensions (already multipied by volume)
    MW_(9),
    N_source(Ns.size(),0.0),
    speciesSources(9),
    Qdot_(Ns.size(),0.0)
{
    Info<< "Creating Soot Model Solver \n\n" << endl;
    
    // Initialize all specie sources to zero
    speciesSources.insert("SOOT", scalarField(Ns.size(),0.0));
    speciesSources.insert("C2H2", scalarField(Ns.size(),0.0));
    speciesSources.insert("O2", scalarField(Ns.size(),0.0));
    speciesSources.insert("OH", scalarField(Ns.size(),0.0));
    speciesSources.insert("H2", scalarField(Ns.size(),0.0));
    speciesSources.insert("CO", scalarField(Ns.size(),0.0));
    speciesSources.insert("H", scalarField(Ns.size(),0.0));
    speciesSources.insert("CO2", scalarField(Ns.size(),0.0));
    speciesSources.insert("H2O", scalarField(Ns.size(),0.0));

    // set up MW table
    MW_.insert("SOOT", composition.W(composition.species()["SOOT"]));
    MW_.insert("C2H2", composition.W(composition.species()["C2H2"]));
    MW_.insert("O2", composition.W(composition.species()["O2"]));
    MW_.insert("OH", composition.W(composition.species()["OH"]));
    MW_.insert("H2", composition.W(composition.species()["H2"]));
    MW_.insert("CO", composition.W(composition.species()["CO"]));
    MW_.insert("H", composition.W(composition.species()["H"]));
    MW_.insert("CO2", composition.W(composition.species()["CO2"]));
    MW_.insert("H2O", composition.W(composition.species()["H2O"]));
}

// Member functions
void Foam::EulerImplicitSystem::ratesOfChange
(
    scalar& r_growth,
    scalar& r_oxidation_O2,
    scalar& r_oxidation_OH,
    scalar& r_gasification_H2O,
    scalar& r_gasification_CO2,
    scalar& r_agglomeration
)
{
    // get relevant state variables for the cell
    scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    scalar cell_T = this->cellState_.thermoProperties()["T"];
    scalar cell_c2h2 = this->cellState_.frozenSpecieMassFractions()["C2H2"];
    scalar cell_Ys = this->cellState_.Ysoot();
    scalar cell_Ns = this->cellState_.Nsoot();
    // Partial pressures
    scalar P_co2 = this->cellState_.frozenSpeciePPressures()["CO2"];
    scalar P_o2 = this->cellState_.frozenSpeciePPressures()["O2"];
    scalar P_h2o = this->cellState_.frozenSpeciePPressures()["H2O"];
    scalar P_oh = this->cellState_.frozenSpeciePPressures()["OH"];
    // Specie concentrations
    const scalar C_c2h2(cell_rho * cell_c2h2 * (1/this->MW_["C2H2"]));
    // Molar Mass of soot
    const scalar M_s = this->MW_["SOOT"];

    // Soot particle diameter and surface area for cell
    // (these are considered uniform in the cell)
    scalar cell_dp(0.0);
    scalar cell_As(0.0);
    if (cell_Ns > Foam::SMALL)
    {
        cell_dp = 
            Foam::pow
            (
                (6.*cell_Ys)/(this->pi * this->rho_s * cell_Ns),
                (1./3.)
            );

        cell_As = this->pi * Foam::sqr(cell_dp) *
            cell_rho * cell_Ns;
    }
      
    // Set constants from papers
    // Kronenburg(2000)
    scalar A_growth(7.5e2);
    scalar T_growth(12100.0);
    // Fletcher(2017)
    scalar A_o2(7.98e-1);
    scalar A_oh(1.89e-3);
    scalar A_co2(3.06e-17);
    scalar A_h2o(6.27e4);
    scalar E_o2(1.77e5);
    scalar E_co2(5.56e3);
    scalar E_h20(2.95e5);
    scalar n_(0.13);
    scalar R_(8.314); //[J/(mol*K)]

    // Now calculate the rates
    // Soot Mass fraction [kmol/(m^3*s)]
    r_growth = A_growth * Foam::exp(-T_growth/cell_T) * cell_As * C_c2h2;
    r_oxidation_O2 = Foam::pow(cell_T, -0.5) * 
        A_o2 * P_o2 * Foam::exp(-E_o2/(R_*cell_T)) *
        (cell_As/M_s);
    r_oxidation_OH = Foam::pow(cell_T, -0.5) * A_oh * P_oh * (cell_As/M_s);
        
    r_gasification_CO2 = A_co2 * Foam::pow(P_co2,0.5) * Foam::sqr(cell_T) *
        Foam::exp(-E_co2/(R_*cell_T)) * (cell_As/M_s);
       
    r_gasification_H2O = A_h2o * Foam::pow(P_h2o,n_) * Foam::pow(cell_T, -0.5) *
        Foam::exp(-E_h20/(R_ * cell_T)) * (cell_As/M_s);
    
    r_agglomeration = 2 * this->Ca * Foam::sqrt(cell_dp) * 
        Foam::sqrt((6.0 * this->sigma * cell_T)/(this->rho_s)) * 
        Foam::sqr(cell_rho * cell_Ns);

    // DEBUG
    r_gasification_CO2 = 0.0;
    r_gasification_H2O = 0.0;
    
} //end ratesOfChange

void Foam::EulerImplicitSystem::explicitStep
(
    const scalarField& Y_initial,
    scalarField& Y_final,
    const label cellNumber,
    const scalar subdt,
    const bool advanceProductSpecies
)
{
    // density in this cell
    const scalar cell_rho = this->cellState_.thermoProperties()["rho"];
    // get relevant state variables for the cell
    // start with per volume mass in cell [kg/m^3]
    const scalar massC2H2(Y_initial[0]*cell_rho);
    const scalar massO2(Y_initial[1]*cell_rho);
    const scalar massOH(Y_initial[2]*cell_rho);
    const scalar massCO2(Y_initial[6]*cell_rho);
    const scalar massH2O(Y_initial[7]*cell_rho);
    
    // rates for soot mass fraction reactions [kmol/(m^3 * s)]
    scalar r_growth(0.0); // growth of soot particle from c2h2 not nucleation
    scalar r_oxidation_O2(0.0);
    scalar r_oxidation_OH(0.0);
    scalar r_gasification_H2O(0.0);    
    scalar r_gasification_CO2(0.0);
    // rate for soot particle number density reaction[1/(m^3 * s)] (not used here)
    scalar r_agglomeration(0.0);
    // calculate the rates
    this->ratesOfChange(r_growth, r_oxidation_O2,r_oxidation_OH, 
    r_gasification_H2O, r_gasification_CO2, r_agglomeration);
    
    // Explicit time step
    // Reactant species
    scalar massC2H2Final = massC2H2 - this->MW_["C2H2"]*r_growth*subdt;
    scalar massO2Final = massO2 - this->MW_["O2"]*r_oxidation_O2*subdt;
    scalar massOHFinal = massOH - this->MW_["OH"]*r_oxidation_OH*subdt;
    scalar massCO2Final = massCO2 - this->MW_["CO2"]*r_gasification_CO2*subdt;
    scalar massH2OFinal = massH2O - this->MW_["H2O"]*r_gasification_H2O*subdt;


    // Assume the overall cell density has not changed 
    // convert back to mass fractions
    Y_final[0]= massC2H2Final / cell_rho;
    Y_final[1]= massO2Final / cell_rho;
    Y_final[2]= massOHFinal / cell_rho;

    // Product species
    if (advanceProductSpecies)
    {
        // get the mass /volume of the product species
        const scalar massH2(Y_initial[3]*cell_rho);
        const scalar massCO(Y_initial[4]*cell_rho);
        const scalar massH(Y_initial[5]*cell_rho);
        
        
        // Add the mass/volume from the current explict step
        // (it's always a source for theses species)
        scalar massH2Final = massH2 + this->MW_["H2"]*(r_growth + r_gasification_H2O)
            *subdt;
        scalar massCOFinal = massCO + this->MW_["CO"]*(r_oxidation_OH +
            r_gasification_H2O + 2*r_gasification_CO2) *subdt;
        scalar massHFinal = massH + this->MW_["H"]*(r_oxidation_OH)*subdt;
        // CO2 is also a reactant in gasification
        massCO2Final += this->MW_["CO2"]*(r_oxidation_O2)*subdt;
        

        // Convert back to mass fraction
        Y_final[3] = massH2Final/cell_rho;
        Y_final[4] = massCOFinal/cell_rho;
        Y_final[5] = massHFinal/cell_rho;
        Y_final[6] = massCO2Final/cell_rho;
        Y_final[7] = massH2OFinal/cell_rho;
    }


    // Soot Stuff

    const scalar massSOOT(Y_initial[8]*cell_rho);

    // total rate of soot consumption for conveneince
    scalar r_soot_consumption = r_oxidation_O2 + r_oxidation_OH + r_gasification_CO2
        + r_gasification_H2O;

    // Take explicit Euler step
    scalar massSootFinal = massSOOT + this->MW_["SOOT"] *
        (2*r_growth  - r_soot_consumption) * subdt;

    // convert back to mass fraction
    Y_final[8] = massSootFinal/cell_rho;
    
    // Need to divide by rho because this rate is for d(N*rho)/dt
    // and we only want to evolve N
    Y_final[9] = Y_initial[9] - (1.0/cell_rho) * r_agglomeration * subdt;


}// end AdvanceFrozenspecies

void Foam::EulerImplicitSystem::advanceToZero
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

void Foam::EulerImplicitSystem::exhaustLowSpecie
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



void Foam::EulerImplicitSystem::correctQdot()
{

    // loop through the species sources that have 
    // been updated 

    forAllIter(HashTable<scalarField>, this->speciesSources, iter)
    {
        word specieName = iter.key();
        label specieIndex = this->composition_.species()[specieName];
        scalarField& sourceField = iter(); //[kg/(m^3 * s)]

        // get the enthalpy of formation for this specie
        scalar hF = this->composition_.Hc(specieIndex);
        
        this->Qdot_ -= hF * sourceField;
    }    
}


//****************** Public member functions ********************//


Foam::tmp<Foam::volScalarField> Foam::EulerImplicitSystem::Qdot()
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
            
        // Reset for next time
        this->Qdot_ = 0.0;

        return tQdot;
}

Foam::tmp<Foam::fvScalarMatrix> Foam::EulerImplicitSystem::sourceN
(
    const volScalarField& N_soot
)
{
    tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(N_soot, N_source_dims));

    // get reference to the fvmatrix inside the tmp
    fvScalarMatrix& fvm = tfvm.ref();

    fvm.source() = -this->N_source;

    // reset field to zero in preparation for next time step.
    this->N_source = 0.0;

    return tfvm;
}


Foam::tmp<Foam::volScalarField> Foam::EulerImplicitSystem::sourceY
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

    // get the reference to fvMatrix that tmp is wrapped around
    volScalarField& sourceY = tSourceY.ref();

    if (this->speciesSources.found(specieName))
    {
        sourceY.primitiveFieldRef() = this->speciesSources[specieName];
        
        // reset field to zero in preparation for next time step.
        this->speciesSources[specieName] = 0.0;
        return tSourceY;
    }
    else
    {
        // return empty source because this specie is unrelated to the soot model
        return tSourceY;
    }
}
 

// Main function for calculating the sources.
// Manages substeps of the main time step dt
void Foam::EulerImplicitSystem::updateSources
(
    const scalar dt,
    const scalarField meshVolumes,
    const label nSubSteps,
    const scalar relTol
)
{  

    // Update the fields stored in cellState to new time
    this->cellState_.updateCellStateFields(); 

    // Local storage for mass fractions of species and soot variables
    scalarField Y_current(10,0.0);
    // Create another field for storage of Mass Fractions/Ns
    // at the end of the sub time step
    scalarField Y_final(Y_current.size(), 0.0);
    
    // Index the species list with the names of the species
    List<word> frozenSpecieNames(8,"none");
    frozenSpecieNames[0] = word("C2H2");
    frozenSpecieNames[1] = word("O2");
    frozenSpecieNames[2] = word("OH");
    frozenSpecieNames[3] = word("H2");
    frozenSpecieNames[4] = word("CO");
    frozenSpecieNames[5] = word("H");
    frozenSpecieNames[6] = word("CO2");
    frozenSpecieNames[7] = word("H2O");

    // Loop through cells
    forAll(this->thermo_.rho().ref(), cellNumber)
      {   
	// Info << "cellNumber: " << cellNumber << endl;
       
        // set default initial sub time step for this cell
        scalar subdt = dt/nSubSteps;

        // Grab the current field data for this cell
        // this will update the cellState thermo data too
        this->cellState_.updateCellState(cellNumber);
        
        // Local storage of frozen specie mass fractions.
        Y_current[0] = this->cellState_.frozenSpecieMassFractions()["C2H2"];
        Y_current[1] = this->cellState_.frozenSpecieMassFractions()["O2"];
        Y_current[2] = this->cellState_.frozenSpecieMassFractions()["OH"];
        Y_current[3] = this->cellState_.frozenSpecieMassFractions()["H2"];
        Y_current[4] = this->cellState_.frozenSpecieMassFractions()["CO"];
        Y_current[5] = this->cellState_.frozenSpecieMassFractions()["H"];
        Y_current[6] = this->cellState_.frozenSpecieMassFractions()["CO2"];
        Y_current[7] = this->cellState_.frozenSpecieMassFractions()["H2O"];
        Y_current[8] = this->cellState_.Ysoot();
        Y_current[9] = this->cellState_.Nsoot();

        // if (cellNumber == 966)
        // {
        //     Info << "Thermo Properties: \n" << 
        //         cellState_.thermoProperties() << endl;
        //     Info << "Mass fractions: \n" << 
        //         cellState_.frozenSpecieMassFractions() <<
        //         "\nSOOT:" << cellState_.Ysoot() <<
        //     "\nNs:" << cellState_.Nsoot() << endl;
        //     Info << "Mixture MM: " << composition_.W().ref()[cellNumber] << endl;
        // }

        // Store a constant version of these initial mass fractions
        const scalarField Y_initial(Y_current);
        Y_final = Y_initial;

        // Time that integration has progressed in this cell
        scalar integratedTime(0.0);
        scalar smallTimeStep(0.0);
        while(integratedTime < dt)
        {
            // Advance, explicitly, local cell copy of species an Ns
            this->explicitStep(Y_current, Y_final, cellNumber, subdt, true);

            // It's possible that the explicit step may be too long 
            // and we exhaust the supply of a reactant specie
            if (min(Y_final) < 0.0)
            {
                // Integrates only until the first specie to reach zero does 
                // so. Then calculates how much time that would take, that is
                // 'smallTimeStep'. Y_final is the whole specie list evolved
                // from Y_current state over time step smalltimeStep
                this->exhaustLowSpecie(subdt, smallTimeStep,
                Y_current, Y_final);
                
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
            cellState_.updateCellState(Y_current, frozenSpecieNames,
            Y_current[8], Y_current[9]);

        } // end while (integratedTime < dt)

        // density in this cell (assumed constant over these steps)
        const scalar cellRho = this->cellState_.thermoProperties()["rho"];

	if 
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
	    Info << "\nCell number: " << cellNumber << endl;
	    Info << "Y_initial: \n " << Y_initial << endl;
	    Info << "Y_final: \n " << Y_final << endl;
			
	    FatalErrorInFunction << "Mass Fraction Above one" << abort(FatalError);
	  }

        if (min(Y_final) < 0.0)
        {
            FatalErrorInFunction << "Negatives" << abort(FatalError);
        }

        // Now set the source terms given the changes in soot mass fraction,
        // frozen specie mass fraction  and soot number density. These 
        // source terms will be used when the respective transport
        // equations are solved
        // NOTE: Multiply by density because the transport
        //  equations are d(rho*Y)/dt not just d(Y)/dt, 
        //  we take rho to be constant here so just pull it out.
        const scalar cellVolume(meshVolumes[cellNumber]);
        // Soot number density
        this->N_source[cellNumber] = 
	  (Y_final[9] - Y_initial[9])*cellVolume*cellRho/dt;
        // Soot mass fraction
        this->speciesSources["SOOT"][cellNumber] = 
            (Y_final[8] - Y_initial[8])*cellRho/dt;

        // //  Other species
        // NOTE: Multiply by density because the transport
        //  equations are d(rho*Y)/dt not just d(Y)/dt
        // and rho is constant over this step


        // Reactant species 
        this->speciesSources["C2H2"][cellNumber] = 
            (Y_final[0] - Y_initial[0])*cellRho/dt;
        this->speciesSources["O2"][cellNumber] = 
            (Y_final[1] - Y_initial[1])*cellRho/dt;
        this->speciesSources["OH"][cellNumber] = 
            (Y_final[2] - Y_initial[2])*cellRho/dt;
        this->speciesSources["H2O"][cellNumber] = 
            (Y_final[7] - Y_initial[7])*cellRho/dt;
        // Product species
        this->speciesSources["H2"][cellNumber] = 
            (Y_final[3] - Y_initial[3])*cellRho/dt;
        this->speciesSources["CO"][cellNumber] = 
            (Y_final[4] - Y_initial[4])*cellRho/dt;
        this->speciesSources["H"][cellNumber] = 
            (Y_final[5] - Y_initial[5])*cellRho/dt;
	// Both a product and a reactant
        this->speciesSources["CO2"][cellNumber] = 
            (Y_final[6] - Y_initial[6])*cellRho/dt;

    }// end loop through cells
    
    // based on the updated species sources calculate the enthalpy source
    this->correctQdot();
}

// **************** Access Functions ************************//
