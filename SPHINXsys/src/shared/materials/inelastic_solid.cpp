#include "inelastic_solid.h"
#include "base_particles.hpp"

namespace SPH
{
	//=================================================================================================//
	void HardeningPlasticSolid::initializePlasticParameters()
	{
		base_particles_->registerVariable(inverse_plastic_strain_, "InversePlasticRightCauchyStrain", 
											[&](size_t i) -> Matd { return Matd::Identity(); });
		base_particles_->registerVariable(hardening_parameter_, "HardeningParameter");
		base_particles_->addVariableToRestart<Matd>("InversePlasticRightCauchyStrain");
		base_particles_->addVariableToRestart<Real>("HardeningParameter");
	}
	//=================================================================================================//
	void HardeningPlasticSolid::assignBaseParticles(BaseParticles *base_particles)
	{
		ElasticSolid::assignBaseParticles(base_particles);
		initializePlasticParameters();
	}
	//=================================================================================================//
	Matd HardeningPlasticSolid::PlasticConstitutiveRelation(const Matd &F, size_t index_i, Real dt)
	{
		Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
		Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
		Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
		Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
		Real deviatoric_PK_norm = deviatoric_PK.norm();
		Real trial_function = deviatoric_PK_norm -
							  sqrt_2_over_3_ * (hardening_modulus_ * hardening_parameter_[index_i] + yield_stress_);
		if (trial_function > 0.0)
		{
			Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
			Real relax_increment = 0.5 * trial_function / (renormalized_shear_modulus + hardening_modulus_ / 3.0);
			hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
			deviatoric_PK -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_PK / deviatoric_PK_norm;
			Matd relaxed_be = deviatoric_PK / G0_ + normalized_be_isentropic * Matd::Identity();
			normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
		}
		Matd inverse_F = F.inverse();
		Matd inverse_F_T = inverse_F.transpose();
		inverse_plastic_strain_[index_i] = inverse_F * normalized_be * inverse_F_T;

		return (deviatoric_PK + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_F_T;
	}
	//=================================================================================================//
	void PlasticSolidforPD::initializePlasticParameters()
	{
		base_particles_->registerVariable(isotropic_hardening_ISV1_, "IsotropicHardeningParam");
		base_particles_->registerVariable(kinematic_hardening_ISV2_, "KinematicHardeningParam");
	}
	//=================================================================================================//
	void PlasticSolidforPD::assignBaseParticles(BaseParticles* base_particles)
	{
		ElasticSolid::assignBaseParticles(base_particles);
		initializePlasticParameters();
	}
	//=================================================================================================//
	Real J2PlasticityforPD::YieldFunc(const Matd& stress, const Real& iso_q, const Matd& kin_q)
	{
		return 0;
	}
	//=================================================================================================//
	Matd J2PlasticityforPD::PlasticConstitutiveRelation(const Matd& G, const Matd& stress_old, Matd& epsilon_p, size_t index_i, Real dt)
	{
		/** Hughes-Winget incremental objectivity **/
		//Symmetric part of G: rate of deformation tensor / rate of strain tensor
		Matd Gsymm = (G + G.transpose()) * 0.5;
		//Antisymmetric(skew) part of G: rate of rotation tensor / spin tensor / vorticity tensor
		Matd Gskew = (G - G.transpose()) * 0.5;
		//Coordinate transformation tensor Q: rotation of principal direction
		Matd R = Matd::Identity() - Gskew * 0.5;
		Matd Q = Matd::Identity() + R.inverse() * Gskew;

		/** Elastic Predictor **/
		Matd delta_sigma = lambda0_ * Gsymm.trace() * Matd::Identity()
			+ 2.0 * G0_ * Gsymm;//regard Gsymm as shear strain tensor
		//Matd delta_sigma = 0.5 * lambda0_ * Gsymm.trace() * Matd::Identity()
		//	+ G0_ * Gsymm;//regard Gsymm as engineering shear strain!
		//delta_sigma.diagonal() = delta_sigma.diagonal() * 2.0;//shear strain tensor

		Matd trial_sigma = Q * stress_old * Q.transpose() + delta_sigma;
		Matd trial_epsilon_p = Q * epsilon_p * Q.transpose();
		//Predict the ISVs
		Real trial_isoHardening_ISV1 = isotropic_hardening_ISV1_[index_i];
		Matd trial_kinHardening_ISV2 = Q * kinematic_hardening_ISV2_[index_i] * Q.transpose();
		//Calculate Hardening variables 
		Real trial_isoHardening_q = IsoHardeningFunc(trial_isoHardening_ISV1);
		Matd trial_kinHardening_q = KinHardeningFunc(trial_kinHardening_ISV2);

		//intermediate variables
		Matd trial_relative_stress = trial_sigma - trial_kinHardening_q;
		Real norm1 = trial_relative_stress.trace() * OneOverDimensions;
		Matd dev_eta = trial_relative_stress - norm1 * Matd::Identity();
		Real dev_eta_norm = dev_eta.norm();

		/** Yield Function **/
		//Real trial_function = YieldFunc(trial_sigma, trial_isoHardening_q);
		Real trial_function = dev_eta_norm - 
			sqrt_2_over_3_ * (trial_isoHardening_q + yield_stress_);
				
		if (trial_function > TinyReal)
		{
			/** Plastic Corrector **/
			//delta gamma: plastic consistency increment / parameter / multiplier
			Real relax_increment = 0.5 * trial_function / 
				(G0_ + (isotropic_hardening_modulus_ + kinematic_hardening_modulus_) / 3.0);
			//plastic potential
			Matd flow_direction = dev_eta / dev_eta_norm;
			//update cauchy stress
			trial_sigma -= 2.0 * G0_ * relax_increment * flow_direction;
			//update plastic strain
			trial_epsilon_p += relax_increment * flow_direction;
			//update the ISVs
			trial_isoHardening_ISV1 += sqrt_2_over_3_ * relax_increment;
			trial_kinHardening_ISV2 += relax_increment * flow_direction;			
		}
		
		epsilon_p = trial_epsilon_p;
		isotropic_hardening_ISV1_[index_i] = trial_isoHardening_ISV1;
		kinematic_hardening_ISV2_[index_i] = trial_kinHardening_ISV2;

		return trial_sigma;
	}
}
