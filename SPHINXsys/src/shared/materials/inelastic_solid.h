/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	inelastic_solid.h
 * @brief 	These are classes for define properties of elastic solid materials.
 *			These classes are based on isotropic linear elastic solid.
 * 			Several more complex materials, including neo-hookean, FENE noe-hookean
 *			and anisotropic muscle, are derived from the basic elastic solid class.
 * @author	Xiangyu Hu and Chi Zhang
 */
#pragma once

#include "elastic_solid.h"

namespace SPH
{
	/**
	* @class PlasticSolid
	* @brief Abstract class for a generalized plastic solid
	*/
	class PlasticSolid : public NeoHookeanSolid
	{
	protected:
		Real yield_stress_;

		virtual void initializePlasticParameters() = 0;

	public:
		/** Constructor */
		explicit PlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress)
			: NeoHookeanSolid(rho0, youngs_modulus, poisson_ratio), yield_stress_(yield_stress)
		{
			material_type_name_ = "PlasticSolid";
		};
		virtual ~PlasticSolid(){};

		Real YieldStress() { return yield_stress_; };
		/** compute the stress through deformation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd &deformation, size_t index_i, Real dt = 0.0) = 0;

		virtual PlasticSolid *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class HardeningPlasticSolid
	 * @brief Class for plastic solid with hardening
	 */
	class HardeningPlasticSolid : public PlasticSolid
	{
	protected:
		Real hardening_modulus_;
		const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
		StdLargeVec<Matd> inverse_plastic_strain_; /**< inverse of plastic right cauchy green strain tensor */
		StdLargeVec<Real> hardening_parameter_;	   /**< hardening parameter */

		virtual void initializePlasticParameters() override;

	public:
		/** Constructor */
		explicit HardeningPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real hardening_modulus)
			: PlasticSolid(rho0, youngs_modulus, poisson_ratio, yield_stress), hardening_modulus_(hardening_modulus)
		{
			material_type_name_ = "HardeningPlasticSolid";
		};
		virtual ~HardeningPlasticSolid(){};

		Real HardeningModulus() { return hardening_modulus_; };
		/** assign particles to this material */
		virtual void assignBaseParticles(BaseParticles *base_particles) override;;
		/** compute the stress through deformation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd &deformation, size_t index_i, Real dt = 0.0) override;

		virtual HardeningPlasticSolid *ThisObjectPtr() override { return this; };
	};

	//from here the Plasticity for NosbPD is constructed 

	/**
	* @Created by Haotian Shi from SJTU
	* @class PlasticSolidforPD
	* @brief Abstract class for a generalized plastic solid for PD combining the isotropic-kinematic hardening
	*/
	class PlasticSolidforPD : public HughesWingetSolid
	{
	protected:
		const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
		Real yield_stress_;
		Real isotropic_hardening_modulus_;
		Real kinematic_hardening_modulus_; //default

		//StdLargeVec<Matd>& plastic_strain_;
		StdLargeVec<Real> isotropic_hardening_ISV1_;
		StdLargeVec<Matd> kinematic_hardening_ISV2_; //default

		virtual void initializePlasticParameters();

	public:
		/** Constructor */
		explicit PlasticSolidforPD(Real rho0, Real youngs_modulus, Real poisson_ratio, 
			Real yield_stress, Real isotropic_hardening_modulus = 0.0, Real kinematic_hardening_modulus = 0.0)
			: HughesWingetSolid(rho0, youngs_modulus, poisson_ratio), yield_stress_(yield_stress), 
			isotropic_hardening_modulus_(isotropic_hardening_modulus), kinematic_hardening_modulus_(kinematic_hardening_modulus)
		{
			material_type_name_ = "PlasticSolidforPD";			
		};
		virtual ~PlasticSolidforPD() {};

		Real YieldStress() { return yield_stress_; };
		Real IsotropicHardeningModulus() { return isotropic_hardening_modulus_; };
		Real KinematicHardeningModulus() { return kinematic_hardening_modulus_; };

		/** assign particles to this material */
		virtual void assignBaseParticles(BaseParticles* base_particles) override;;
		/* compute yield function */
		virtual Real YieldFunc(const Matd & stress, const Real& iso_q = 0.0, const Matd& kin_q = Matd::Zero()) = 0;
		/** compute the stress through deformation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd& vel_grad , const Matd& stress_old, Matd& epsilon_p, size_t index_i, Real dt = 0.0) = 0;

		virtual PlasticSolidforPD* ThisObjectPtr() override { return this; };
	};
	/**
	* @Created by Haotian Shi from SJTU
	* @class J2PlasticityforPD
	* @brief Associative plasticity for PD under J2 theory, combine isotropic-kinematic hardening, rate-independent
	*/
	class J2PlasticityforPD : public PlasticSolidforPD
	{
	protected:
		//Parameter for Saturation-typed hardening function
		Real sigmaY_infinite_;
		Real sigmaY_0_;
		Real index_kesi_;

	public:
		/** Constructor */
		explicit J2PlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress,
			Real isotropic_hardening_modulus = 0.0, Real kinematic_hardening_modulus = 0.0)
			: PlasticSolidforPD(rho0, youngs_modulus, poisson_ratio, yield_stress, 
				isotropic_hardening_modulus, kinematic_hardening_modulus), 
			sigmaY_infinite_(1.589 * yield_stress), sigmaY_0_(yield_stress), index_kesi_(16.93)
		{
			material_type_name_ = "J2PlasticityforPD";
		};
		virtual ~J2PlasticityforPD() {};

		/* compute yield function */
		virtual Real YieldFunc(const Matd& stress, const Real& iso_q = 0.0, const Matd& kin_q = Matd::Zero()) override;
		/* compute hardening varibles by hardening constitution */
		inline Real IsoHardeningFunc(const Real& ISV1) 
		{
			Real exponent = exp(-index_kesi_ * ISV1);
			return isotropic_hardening_modulus_ * ISV1 + (sigmaY_infinite_ - sigmaY_0_) * (1 - exponent);
		};
		inline Matd KinHardeningFunc(const Matd& ISV2) { return kinematic_hardening_modulus_ * ISV2; };
		/** compute the stress through deformation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd& vel_grad, const Matd& stress_old, Matd& epsilon_p, size_t index_i, Real dt = 0.0) override;

	};
	/**
	* @Created by Haotian Shi from SJTU
	* @class DruckerPragerPlasticityforPD
	* @brief Associative plasticity for PD under Drucker-Prager theory, only isotropic hardening "k" being concerned, rate-independent
	*/
	class DruckerPragerPlasticityforPD : public PlasticSolidforPD
	{
	protected:
		//
		const Real sqrt2 = sqrt(2);
		Real alpha_phi_;
		Real k_phi_;		
		//Parameter for Saturation-typed hardening function
		Real sigmaY_infinite_;
		Real sigmaY_0_;
		Real index_kesi_;

	public:
		/** Constructor */
		explicit DruckerPragerPlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real alpha_phi,
			Real isotropic_hardening_modulus = 0.0, Real kinematic_hardening_modulus = 0.0)
			: PlasticSolidforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
				isotropic_hardening_modulus, kinematic_hardening_modulus), k_phi_(yield_stress), alpha_phi_(alpha_phi),
			sigmaY_infinite_(1.589 * yield_stress), sigmaY_0_(yield_stress), index_kesi_(16.93)
		{
			material_type_name_ = "DruckerPragerPlasticityforPD";
			
		};
		virtual ~DruckerPragerPlasticityforPD() {};

		/* compute yield function */
		virtual Real YieldFunc(const Matd& stress, const Real& iso_q = 0.0, const Matd& kin_q = Matd::Zero()) override;
		/* compute hardening varibles by hardening constitution */
		inline Real IsoHardeningFunc(const Real& ISV1)
		{
			Real exponent = exp(-index_kesi_ * ISV1);
			return isotropic_hardening_modulus_ * ISV1 + (sigmaY_infinite_ - sigmaY_0_) * (1 - exponent);
		};
		inline Matd KinHardeningFunc(const Matd& ISV2) { return kinematic_hardening_modulus_ * ISV2; };
		/** compute the stress through deformation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd& vel_grad, const Matd& stress_old, Matd& epsilon_p, size_t index_i, Real dt = 0.0) override;

	};
}
