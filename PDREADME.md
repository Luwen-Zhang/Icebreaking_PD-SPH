# Preface

Before getting to know the following modifications on establishing NOSB-PD method based on SPHinXsys, you must have a preliminary knowledge of SPHinXsys (https://www.sphinxsys.org/).

# Influence Function and its Radius

..\SPHinXsys-master-source\SPHINXsys\src\shared\kernels\kernel_wenland_c2.h  
Derived `class Winfunc` from `class Kernel` to introduce the *influence function* of PD, variable `kernel_size_` maintain 2.0.  
Note: Calculating the derivative of the influence function.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\bodies\solid_body.h  
Derived `class PDBody` from `class RealBody`, implement the unique horizon size of PD by function `defineAdaptation<SPHAdaptation>(1.5075)`.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\adaptations\adaptation.h  
Add variable `std::string body_name_`, use the first 6 characters *"PDBody"* as the unique identifier for building `class Winfunc`.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\adaptations\adaptation.cpp  
Add `body_name_ = sph_body.getName()` in the constructor of `class SPHAdaptation`.  
Constructing `class Winfunc` based on  the unique identifier *"PDBody"*, 
```cpp
if (body_name_ == "PDBody") {
			kernel_ptr_ = makeUnique<Winfunc>(h_ref_);
		}
```

# Data Containers and Data Structures

..\SPHinXsys-master-source\SPHINXsys\src\shared\particles\solid_particles.h  
Declare `class PDParticles : public ElasticSolidParticles`, common member variables `particleLive，bond，damage` are added for PD  
Declare `class NosbPDParticles : public PDParticles`, exclusive variables for NOSB-PD added here.  
Declare `class NosbPDPlasticParticles : public NosbPDParticles`, exclusive variables for NOSB-PD with elastoplastic constitutive model added here.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\particles\solid_particles.cpp  
Define the constructor of the above `class`, initialize the member variables, define the necessary member functions.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_neighborhood\neighborhood.h  
Add member variable `StdLargeVec<bool> bondLive_` to introduce the concept *bond*.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_neighborhood\neighborhood.cpp  
Initialize variable `bondLive_` in  member function `class NeighborBuilder::createNeighbor` and `NeighborBuilder::initializeNeighbor`  

# Constitutive Model

..\SPHinXsys-master-source\SPHINXsys\src\shared\materials\elastic_solid.h  
Declare `class HughesWingetSolid : public LinearElasticSolid` to implement an incremental constitutive with finite strain and finite rotation.  
The constructor and member function are defined in the corresponding *.cpp* file.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\materials\inelastic_solid.h  
Declare `class PlasticSolidforPD : public HughesWingetSolid`, which is the base class of elastoplastic constitutive model.  
The *Yield strength, Isotropic hardening modulus and its ISVs, Kinematic hardening modulus and its ISVs* are introduced by member variables.  
Virtual functions are prepared for *Yield function, Isotropic hardening, Kinematic hardening, elastoplastic constitutive relationship*.  

Declare `class J2PlasticityforPD : public PlasticSolidforPD`. Based on the decomposition of strain and incremental objectivity, the J2 plasticity with associated flow rules is realized concerning isotropic hardening, kinematic hardening. This model mainly applys to the macroscopic simulations of metal (Note: Only available for 3D).  

Declare `class DruckerPragerPlasticityforPD : public PlasticSolidforPD`. Based on the decomposition of strain and incremental objectivity, the Drucker-Prager plasticity with associated flow rules is realized concerning isotropic hardening. This model mainly applys to the macroscopic simulations of ice (Note: Only available for 3D).  

All the constructors and member functions of above classes have been defined in the corresponding *.cpp* files.  

# Calculation Methods of Governing Equations

## Data Containers

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h  
Define `typedef DataDelegateSimple<NosbPDParticles> NosbPDSolidDataSimple`.  
Define `typedef DataDelegateInner<NosbPDParticles> NosbPDSolidDataInner`  
Define `typedef DataDelegateInner<NosbPDPlasticParticles> NosbPDPlasticSolidDataInner` 

## Spatial Integration

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h  
Declare `class PDTimeStepInitialization`, the member function is overridden in the corresponding *.cpp* file to exert gravity.  
Declare `class NosbPDShapeMatrix`, the member function is overridden in the corresponding *.cpp* file to calculate *shape tensor K*.  
Declare `class NosbPDSecondStep`, the member function is overridden in the corresponding *.cpp* file to calculate *deformation gradient F, Euler velocity gradient G and intermediate variables*. Then the *Cauchy stress, PK1, force density T* are calculated.  
Declare `class NosbPDSecondStepPlastic`, the member function is overridden in the corresponding *.cpp* file. The plastic version if `class NosbPDSecondStep` in which the elastoplastic constitutive model is called.  
Declare `class NosbPDThirdStep`, the member function is overridden in the corresponding *.cpp* file. Calculating the first term of dynamic equation of NOSB-PD.  
Declare `class LittleWoodHourGlassControl`, the member function is overridden in the corresponding *.cpp* file. Hourglass control based on a background bond force.  
Declare `class PairNumericalDampingforPD`, the member function is overridden in the corresponding *.cpp* file. The Kelvin-Voigt damper with artificial viscosity in TL-SPH is introduced in NOSB-PD.  

## Criteria for Bond Breaking

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h  
Declare `class NosbPDCheckBondLive`, the member function is overridden in the corresponding *.cpp* file. Base class for bond breaking.  
Declare `class BondBreakByPrinStress`, the member function is overridden in the corresponding *.cpp* file. Criterion based on MAX principal stress.  
Declare `class BondBreakBySigma1andSigma3`, the member function is overridden in the corresponding *.cpp* file. Criterion based on MAX principal stress and MAX shear stress.  
Declare `class BondBreakByPlasticStrain`, the member function is overridden in the corresponding *.cpp* file. Criterion based on MAX principal plastic strain.  

## Time Integration

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\elastic_dynamics.h  
Declare `class NosbPDFirstStep`, the member function is overridden in the corresponding *.cpp* file. Position-Verlet scheme for NOSB-PD, the first half step.  
Declare `class NosbPDFourthStep`, the member function is overridden in the corresponding *.cpp* file. Position-Verlet scheme for NOSB-PD, the second half step.  
Declare `class ADRFirstStep`, the member function is overridden in the corresponding *.cpp* file. Adaptive Dynamic Relaxation (ADR) for quasi-static problems: the first step.  
Declare `class ADRSecondStep`, the member function is overridden in the corresponding *.cpp* file. Adaptive Dynamic Relaxation (ADR) for quasi-static problems: the second step.  
Declare `class NosbPDFourthStepWithADR`, the member function is overridden in the corresponding *.cpp* file. Position-Verlet scheme for NOSB-PD, the second half step with ADR.  

## FSI Algorithm

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_neighborhood\neighborhood.cpp  
Kernel selection strategy for the interaction between *PDBody* and *WaterBody* is added in the constructor of `class NeighborBuilderContact`.  
```cpp
std::string source_name = body.getName().substr(0, 6);
		std::string target_name = contact_body.getName().substr(0, 6);
		if (source_name == "PDBody" && target_name == "WaterB") {
			kernel_ = target_kernel;
		}
		else if (source_name == "WaterB" && target_name == "PDBody") {
			kernel_ = source_kernel;
		}
		else {
			kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;
		}		
```
Reason of modification: Based on the principle of multi-resolution FSI, when the two phases are in different resolutions, simply search for fluid-solid particle pairs based on the smooth length *h* of the fluid domain to achieve multi-resolution fluid-structure coupling. Assuming the average particle distance of solid *ds* is smaller than the fluid *df*, in the original coupling of WCSPH (*h=1.3*) and TL-SPH (*h=1.1*), one can simply find the searching radius by comparing *1.3df > 1.1ds*. However, in the coupling of WCSPH (*h=1.3*) and NOSB-PD (*h=1.5015*), the relationship *1.3df > 1.5015ds* may not always hold true. Therefore, here, we definitely determine the fluid domain by the name of the classes.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\fluid_dynamics\fluid_dynamics_complex.h  
Declare `class BaseIntegration1stHalfWithWallforPD`, the member function is overridden in the corresponding *.cpp* file. FSI for PD. source term: pressure gradient, density diffusion.  
Declare `class BaseIntegration2ndHalfWithWallforPD`, the member function is overridden in the corresponding *.cpp* file. FSI for PD. source term: pressure diffusion, velocity divergence.  

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\solid_dynamics\fluid_structure_interaction.h  
Declare `class BasePressureForceAccelerationFromFluidforPD`, the member function is overridden in the corresponding *.cpp* file. Pressure on FSI interface is monitored for output.  

## Pre-Processing

..\SPHinXsys-master-source\SPHINXsys\src\shared\particle_dynamics\base_local_dynamics.h  
Declare `class BaseSPHBodyRotation`, Rigid body rotation in a coordinate plane for SPHBody particles, by the *orthogonal tensor Q*

## Post-Processing

..\SPHinXsys-master-source\SPHINXsys\src\shared\particles\solid_particles_variable.h  
Define `typedef DataDelegateSimple<NosbPDParticles> PDSolidDataSimple`.  
Define `typedef DataDelegateSimple<NosbPDPlasticParticles> PDPlasticSolidDataSimple`.  

Declare `class VonMisesStressforPD` to output the equivalent Mises stress.  
Declare `class VonMisesPlasticStrainforPD` to output the equivalent Mises strain (Note: Only available for 3D).  

..\SPHinXsys-master-source\SPHINXsys\src\for_2D_build\particles\solid_particles_supplementary.cpp  
..\SPHinXsys-master-source\SPHINXsys\src\for_3D_build\particles\solid_particles_supplementary.cpp  
Override the member functions in the corresponding *.cpp* files of 2D and 3D, respectively.  




































