## Libraries for use in OpenFOAM 5.x

### Description of Contents

Many of these libraries requre specialized solvers to be used, see the corresponding repository for 
[OpenFOAM 5.x Solvers](https://github.com/cdunn6754/OpenFOAM_5.x_Applications) for examples. For instance to use
the Two Equation soot model library, TwoEquationSoot, you need modify your solver to create the soot mass fraction and
particle number density fields with certain names. That is all taken care of in the 
[SootTarFoam solver](https://github.com/cdunn6754/OpenFOAM_5.x_Applications/tree/master/lagrangian/SootTarFoam).

#### Lagrangian:

Collection of libraries that relate to modifications within the Lagrangian classes of OpenFOAM

  * AnalyticalSFORDevolatilization
    A new devolatilization model. Almost identical to the OF 
	  singleKineticRateDevolatilization model. Here we use an analytical
	  solution to the mass loss rate equation rather than the explicit integration
	  used in singleKineticRateDevolatilization.
    
  * PCReactingMultiphase
    This is an alternative to the OF ReactingMultiphase parcel and cloud classes. It includes the new 
	  devolatilization model that permits tar to be retained within a 
	  particle, the new PCReactingMultiphaseParcel class and the new
	  PCReactingMultiphaseCloud class. The parcel and cloud classes 
	  were developed because we need the parcels to track additional 
	  information about primary and secondary tar fractions. Other than
	  that they are very similar to the OF ReactingMultiphase classes.

  * smithIntrinsicSurfaceReactionModel
    New surface reaction model for coal particles. Very similar to the 
	  OF model 'COxidationIntrinsicRate' model. This version is implemented
	  as described in 
    [1](https://www.sciencedirect.com/science/article/pii/S0010218012002994?via%3Dihub).
    
   * nthOrderDevolatilization
     Another new devolatlilization model. Similar to the single first order
	   rate equation models but has an additional exponential parameter. This
	   rate equation is used by C3M for their CPD predictions.
    
#### sootModeling:
  	
This directory collects all of the libraries that implement the 
Two Equation soot model and the tracked tar model. Since 
the rate constants are hard-coded 
I ended up making different versions of each for different cases.
    
  * TwoEquationSoot
    Used with SootTarFoam solvers to model soot mass fraction and particle number 
    density with two transport equations (see Brown 1998). Includes a soot surface
    growth reaction to account for gas phase acetylene addition. Fully accounts for
    chemical and energetic coupling between this model and the main gas phase chemistry/enthalpy.
    
   * NoGrowthTwoEquationSoot
     This is the same as TwoEquationSoot but the soot growth reaction rate
	   is set to 0.
     
   * SootTarModel
     Similar to the TwoEquationSoot model but implementing the tar transport equation (Ma 1996)
     rather than the two soot equations. Again fully accounts for chemical and energetic coupling.
     
   * PCFSootTarModel
     Same as the SootTarModel but it uses rates for soot formation and 
	   tar cracking that were determined from the [PCF](https://github.com/cdunn6754/PCCLConversion)
     rather than Ma's (Ma 1996) experiment/paper.
    
#### thermoPhysicalModels:

Collection of classes that relate to the thermo-physical models in OpenFOAM
  * radiation
    Composed of two models. First is the ternaryAbsorptionEmission model. 
	  It is an extension of the OF binaryAbsorptionEmission model. Instead
	  of just capturing two sources for A/E we now can specify three.

	  The other model here is the greyMeanSootAbsorptionEmission model. This 
	  is the way to calculate the soot contribution which can thenbe included as 
	  the third A/E model in the ternaryModel above. It calculates the A/E based
	  on the soot volume fraction per Xu (2017).
    
  * sootReactionThermo
     This model corrects the calculation of density for the 
	   presence of soot by applying a multiphase density approach. It 
	   is pretty complicated, there is a new version of a few of the 
	   lower level mixture classes. The entire purpose is to exclude 
	   soot from the calculation of gas phase density through the IGL and then
	   sum the soot density and gas density scaled by thier respective volume 
	   fractions, i.e. proper multiphase density calculation. 
    
    
### Environment Set-up for Building:
Libraries for use with OpenFOAM 5.x (OF). You need to have a few important OF specific enviroment variables set
before building this project. Most importantly you will need to ensure the 
`FOAM_USER_LIBBIN`  variable is set. It should point to something like
`~/OpenFOAM/<username>/platforms/<your platform>/lib` where `<username>` is your Linux username and 
`<your platform>` is a description of your system, mine is `linux64Gcc48DPint64Opt`, the OF installation should
create that for you. This is where the library
files will live after you build them. 

### Building the Libraries
With the enviroment set and assuming you have OF-5.x installed properly you should
be able to clone this repo into your source code folder, `~/OpenFOAM/<username>/src/`. Then you can
build it with the OF command `wmake` from the root directory or any of the subdirectories that have a `Make/`
folder. But note that some of the libraries depend on eachother so it may be best to just build the whole thing.
