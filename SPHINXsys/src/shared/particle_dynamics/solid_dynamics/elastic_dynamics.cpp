#include "elastic_dynamics.h"
#include "general_dynamics.h"

#include <numeric>


namespace SPH
{
//=========================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real CFL)
			: LocalDynamicsReduce<Real, ReduceMin>(sph_body, Real(MaxRealNumber)),
			  ElasticSolidDataSimple(sph_body), CFL_(CFL),
			  vel_(particles_->vel_), acc_(particles_->acc_), acc_prior_(particles_->acc_prior_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()),
			  c0_(particles_->elastic_solid_.ReferenceSoundSpeed()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return CFL_ * SMIN(sqrt(smoothing_length_ / ((acc_[index_i] + acc_prior_[index_i]).norm() + TinyReal)),
							   smoothing_length_ / (c0_ + vel_[index_i].norm()));
		}
		//=================================================================================================//
		ElasticDynamicsInitialCondition::ElasticDynamicsInitialCondition(SPHBody &sph_body)
			: LocalDynamics(sph_body),
			  ElasticSolidDataSimple(sph_body),
			  pos_(particles_->pos_), vel_(particles_->vel_) {}
		//=================================================================================================//
		UpdateElasticNormalDirection::UpdateElasticNormalDirection(SPHBody &sph_body)
			: LocalDynamics(sph_body),
			  ElasticSolidDataSimple(sph_body),
			  n_(particles_->n_), n0_(particles_->n0_), F_(particles_->F_) {}
		//=================================================================================================//
		DeformationGradientBySummation::
			DeformationGradientBySummation(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), ElasticSolidDataInner(inner_relation),
			  pos_(particles_->pos_), B_(particles_->B_), F_(particles_->F_) {}
		//=================================================================================================//
		void DeformationGradientBySummation::interaction(size_t index_i, Real dt)
		{
			Vecd &pos_n_i = pos_[index_i];

			Matd deformation = Matd::Identity();
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				deformation -= (pos_n_i - pos_[index_j]) * gradW_ijV_j.transpose();
			}

			F_[index_i] = deformation * B_[index_i];
		}
		//=================================================================================================//
		BaseElasticIntegration::
			BaseElasticIntegration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), ElasticSolidDataInner(inner_relation),
			  rho_(particles_->rho_), mass_(particles_->mass_),
			  pos_(particles_->pos_), vel_(particles_->vel_), acc_(particles_->acc_),
			  B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_) {}
		//=================================================================================================//
		BaseIntegration1stHalf::
			BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
			: BaseElasticIntegration(inner_relation),
			  elastic_solid_(particles_->elastic_solid_),
			  acc_prior_(particles_->acc_prior_)
		{
			rho0_ = particles_->elastic_solid_.ReferenceDensity();
			inv_rho0_ = 1.0 / rho0_;
			smoothing_length_ = sph_body_.sph_adaptation_->ReferenceSmoothingLength();
		}
		//=================================================================================================//
		void BaseIntegration1stHalf::update(size_t index_i, Real dt)
		{
			vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
		}
		//=================================================================================================//
		Integration1stHalf::
			Integration1stHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration1stHalf(inner_relation)
		{
			particles_->registerVariable(stress_PK1_B_, "CorrectedStressPK1");
			numerical_dissipation_factor_ = 0.25;
		}
		//=================================================================================================//
		void Integration1stHalf::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_[index_i] = rho0_ / F_[index_i].determinant();
			// obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress
			// it seems using reproducing correction here increases convergence rate near the free surface
			stress_PK1_B_[index_i] = F_[index_i] * elastic_solid_.StressPK2(F_[index_i], index_i) * B_[index_i];
		}
		//=================================================================================================//
		void Integration1stHalf::interaction(size_t index_i, Real dt)
		{
			// including gravity and force from fluid
			Vecd acceleration = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];
				Real dim_r_ij_1 = Dimensions / r_ij;
				Vecd pos_jump = pos_[index_i] - pos_[index_j];
				Vecd vel_jump = vel_[index_i] - vel_[index_j];
				Real strain_rate = dim_r_ij_1 * dim_r_ij_1 * pos_jump.dot(vel_jump);
				Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
				Matd numerical_stress_ij =
					0.5 * (F_[index_i] + F_[index_j]) * elastic_solid_.PairNumericalDamping(strain_rate, smoothing_length_);
				acceleration += inv_rho0_ *  inner_neighborhood.dW_ijV_j_[n] *
								(stress_PK1_B_[index_i] + stress_PK1_B_[index_j] +
								 	numerical_dissipation_factor_ * weight * numerical_stress_ij) * e_ij;
			}

			acc_[index_i] = acceleration;
		}
		//=================================================================================================//
		KirchhoffParticleIntegration1stHalf::
			KirchhoffParticleIntegration1stHalf(BaseInnerRelation &inner_relation)
			: Integration1stHalf(inner_relation){};
		//=================================================================================================//
		void KirchhoffParticleIntegration1stHalf::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_[index_i] = rho0_ / F_[index_i].determinant();
			Real J = F_[index_i].determinant();
			Real one_over_J = 1.0 / J;
			rho_[index_i] = rho0_ * one_over_J;
			Real J_to_minus_2_over_dimension = pow(one_over_J, 2.0 * OneOverDimensions);
			Matd normalized_b = (F_[index_i] * F_[index_i].transpose()) * J_to_minus_2_over_dimension;
			Matd deviatoric_b = normalized_b - Matd::Identity() * normalized_b.trace() * OneOverDimensions;
			Matd inverse_F_T = F_[index_i].inverse().transpose();
			// obtain the first Piola-Kirchhoff stress from the Kirchhoff stress
			// it seems using reproducing correction here increases convergence rate
			// near the free surface however, this correction is not used for the numerical dissipation
			stress_PK1_B_[index_i] = ( Matd::Identity() * elastic_solid_.VolumetricKirchhoff(J) +
									  	elastic_solid_.DeviatoricKirchhoff(deviatoric_b) ) * inverse_F_T * B_[index_i];
		}
		//=================================================================================================//
		KirchhoffIntegration1stHalf::
			KirchhoffIntegration1stHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration1stHalf(inner_relation)
		{
			particles_->registerVariable(J_to_minus_2_over_dimension_, "DeterminantTerm");
			particles_->registerVariable(stress_on_particle_, "StressOnParticle");
			particles_->registerVariable(inverse_F_T_, "InverseTransposedDeformation");
		};
		//=================================================================================================//
		void KirchhoffIntegration1stHalf::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			Real J = F_[index_i].determinant();
			Real one_over_J = 1.0 / J;
			rho_[index_i] = rho0_ * one_over_J;
			J_to_minus_2_over_dimension_[index_i] = pow(one_over_J * one_over_J, OneOverDimensions);
			
			inverse_F_T_[index_i] = F_[index_i].inverse().transpose();
			stress_on_particle_[index_i] = inverse_F_T_[index_i] * 
				(elastic_solid_.VolumetricKirchhoff(J) - correction_factor_ * elastic_solid_.ShearModulus() *
				 J_to_minus_2_over_dimension_[index_i] * (F_[index_i] * F_[index_i].transpose()).trace() * OneOverDimensions) 
				+ elastic_solid_.NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_, index_i) * inverse_F_T_[index_i];
		}
		//=================================================================================================//
		void KirchhoffIntegration1stHalf::interaction(size_t index_i, Real dt)
		{
			// including gravity and force from fluid
			Vecd acceleration = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd shear_force_ij = correction_factor_ * elastic_solid_.ShearModulus() *
									  (J_to_minus_2_over_dimension_[index_i] + J_to_minus_2_over_dimension_[index_j]) *
									  (pos_[index_i] - pos_[index_j]) / inner_neighborhood.r_ij_[n];
				acceleration += ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) * inner_neighborhood.e_ij_[n] + shear_force_ij) *
								inner_neighborhood.dW_ijV_j_[n] * inv_rho0_;
			}
			acc_[index_i] = acceleration;
		}
		//=================================================================================================//
		void Integration2ndHalf::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void Integration2ndHalf::interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_n_i = vel_[index_i];

			Matd deformation_gradient_change_rate = Matd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				deformation_gradient_change_rate -= (vel_n_i - vel_[index_j]) * gradW_ij.transpose();
			}

			dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
		}
		//=================================================================================================//
		void Integration2ndHalf::update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		//=================================================================================================//
		NosbPDShapeMatrix::
			NosbPDShapeMatrix(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), NosbPDSolidDataInner(inner_relation),
			Vol_(particles_->Vol_), shape_K_(particles_->shape_K_),
			shape_K_1_(particles_->shape_K_1_) {}
		//=================================================================================================//
		void NosbPDShapeMatrix::interaction(size_t index_i, Real dt)
		{
			Matd local_configuration = Eps * Matd::Identity(); // a small number added to diagonal to avoid divide zero
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration += inner_neighborhood.W_ij_[n] * Vol_[index_j] * (r_ij * r_ij.transpose());
			}
			shape_K_[index_i] += local_configuration;
		}
		//=================================================================================================//
		void NosbPDShapeMatrix::update(size_t index_i, Real dt)
		{
			shape_K_1_[index_i] = shape_K_[index_i].inverse();
		}

		//=================================================================================================//
		NosbPDFirstStep::NosbPDFirstStep(SPHBody& sph_body)
			: LocalDynamics(sph_body),
			NosbPDSolidDataSimple(sph_body), 
			elastic_solid_(particles_->elastic_solid_),	rho_(particles_->rho_), 
			pos_(particles_->pos_), vel_(particles_->vel_), 
			acc_(particles_->acc_), acc_prior_(particles_->acc_prior_), F_(particles_->F_)
		{
			rho0_ = particles_->elastic_solid_.ReferenceDensity();		
		}
		//=================================================================================================//
		void NosbPDFirstStep::update(size_t index_i, Real dt)
		{
			acc_[index_i] = acc_prior_[index_i];
			rho_[index_i] = rho0_ / F_[index_i].determinant();
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}

		//=================================================================================================//
		NosbPDSecondStep::
			NosbPDSecondStep(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), NosbPDSolidDataInner(inner_relation),
			particleLive_(particles_->particleLive_), Vol_(particles_->Vol_),
			pos_(particles_->pos_), vel_(particles_->vel_), F_(particles_->F_), shape_K_1_(particles_->shape_K_1_),
			N_(particles_->N_), N_deltaU_(particles_->N_deltaU_), N_half_(particles_->N_half_),
			F_half_(particles_->F_half_), F_delta_(particles_->F_delta_), F_1_(particles_->F_1_), F_1_half_(particles_->F_1_half_),
			PK1_(particles_->PK1_), T0_(particles_->T0_) {}
		//=================================================================================================//
		void NosbPDSecondStep::interaction(size_t index_i, Real dt)
		{
			if (particleLive_[index_i] == 1) {

				Matd N = Matd::Zero();
				Matd N_deltaU = Matd::Zero();
				Matd N_half = Matd::Zero();

				Vecd eta = Vecd::Zero();
				Vecd eta_U = Vecd::Zero();
				Vecd eta_half = Vecd::Zero();

				const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) 
				{
					size_t index_j = inner_neighborhood.j_[n];

					Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
					if (inner_neighborhood.bondLive_[n]) {
						eta = pos_[index_j] - pos_[index_i];
						eta_U = (vel_[index_j] - vel_[index_i]) * dt;
						eta_half = eta - eta_U * 0.5;
					}
					else {
						eta = r_ij;
						eta_U = Vecd::Zero();
						eta_half = r_ij;
					}
					Matd dN = eta * r_ij.transpose();
					Matd dN_deltaU = eta_U * r_ij.transpose();
					Matd dN_half = eta_half * r_ij.transpose();
					N += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN;
					N_deltaU += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN_deltaU;
					N_half += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN_half;
				}
				N_[index_i] = N;
				N_deltaU_[index_i] = N_deltaU;
				N_half_[index_i] = N_half;
			}
		}
		//=================================================================================================//
		void NosbPDSecondStep::update(size_t index_i, Real dt)
		{
			if (particleLive_[index_i] == 1) {
				F_[index_i] = N_[index_i] * shape_K_1_[index_i];
				F_delta_[index_i] = N_deltaU_[index_i] * shape_K_1_[index_i];
				F_half_[index_i] = N_half_[index_i] * shape_K_1_[index_i];

				F_1_[index_i] = F_[index_i].inverse();
				F_1_half_[index_i] = F_half_[index_i].inverse();

				Real detF = F_[index_i].determinant();
				Real detF_half = F_half_[index_i].determinant();
				if (detF < 0.0 || detF_half < 0.0) {
					std::cout << "Particle_index = " << index_i << "has been disabled since det F < 0"
						<< GlobalStaticVariables::physical_time_ << "dt:" << dt << "\n";
					system("pause");
					exit(0);
				}
				//G: Eulerian Velocity gradient tensor
				Matd G = F_delta_[index_i] * F_1_half_[index_i];
				//Update Cauchy Stress by Hughes-Winget algorithm


			}
		}


	}
}
