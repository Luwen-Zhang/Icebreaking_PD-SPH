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
 * @file 	fluid_dynamics_inner.hpp
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details 	We consider here weakly compressible fluids. The algorithms may be
 * 			different for free surface flow and the one without free surface.
 * @author	Chi ZHang and Xiangyu Hu
 */
#pragma once

#include "fluid_dynamics_inner.h"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(fluid_, fluid_) {

			Real h = inner_relation.getSPHBody().sph_adaptation_->getKernel()->CutOffRadius();
			coeff_acoustic_damper_ = 0.3 * fluid_.ReferenceSoundSpeed() * fluid_.ReferenceDensity() * h;

		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressure(rho_[index_i]);
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		Vecd BaseIntegration1stHalf<RiemannSolverType>::computeNonConservativeAcceleration(size_t index_i)
		{
			Vecd acceleration = acc_prior_[index_i] * rho_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];

				acceleration += (p_[index_i] - p_[index_j]) * dW_ijV_j * e_ij;
			}
			return acceleration / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration =  Vecd::Zero();
			Vecd acoustic_damper =  Vecd::Zero();
			Real rho_dissipation(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];

				acceleration -= (p_[index_i] + p_[index_j]) * dW_ijV_j * e_ij;
				acoustic_damper += (u_div_[index_i] + u_div_[index_j]) * dW_ijV_j * e_ij;
				rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
			}
			acc_[index_i] += acceleration / rho_[index_i] + coeff_acoustic_damper_ * acoustic_damper / rho_[index_i];
			drho_dt_[index_i] = rho_dissipation * rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration2ndHalf<RiemannSolverType>::BaseIntegration2ndHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(fluid_, fluid_),
			  Vol_(particles_->Vol_), mass_(particles_->mass_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Real density_change_rate(0);
			Vecd p_dissipation = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

				Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
				density_change_rate += u_jump * dW_ijV_j;
				p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
			}
			drho_dt_[index_i] += density_change_rate * rho_[index_i];
			u_div_[index_i] = -density_change_rate;
			acc_[index_i] = p_dissipation / rho_[index_i];
		};
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseCompressibleIntegration1stHalf<RiemannSolverType>::BaseCompressibleIntegration1stHalf(BaseInnerRelation& inner_relation)
			: BaseCompressibleIntegration(inner_relation), riemann_solver_(fluid_, fluid_) {

			Real h = inner_relation.getSPHBody().sph_adaptation_->getKernel()->CutOffRadius();
			coeff_acoustic_damper_ = 0.3 * fluid_.ReferenceSoundSpeed() * fluid_.ReferenceDensity() * h;

		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseCompressibleIntegration1stHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			E_[index_i] += dE_dt_[index_i] * dt * 0.5;
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressurebyTamann(rho_[index_i], E_[index_i]);
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseCompressibleIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseCompressibleIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Vecd acoustic_damper = Vecd::Zero();
			Real rho_dissipation(0);
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];

				acceleration -= (p_[index_i] + p_[index_j]) * dW_ijV_j * e_ij;
				acoustic_damper += (u_div_[index_i] + u_div_[index_j]) * dW_ijV_j * e_ij;
				rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
			}
			acc_[index_i] += acceleration / rho_[index_i] + coeff_acoustic_damper_ * acoustic_damper / rho_[index_i];
			drho_dt_[index_i] = rho_dissipation * rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseCompressibleIntegration2ndHalf<RiemannSolverType>::BaseCompressibleIntegration2ndHalf(BaseInnerRelation& inner_relation)
			: BaseCompressibleIntegration(inner_relation), riemann_solver_(fluid_, fluid_),
			Vol_(particles_->Vol_), mass_(particles_->mass_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseCompressibleIntegration2ndHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseCompressibleIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			E_[index_i] += dE_dt_[index_i] * dt * 0.5;
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseCompressibleIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Real density_change_rate(0);
			Real energy_change_rate(0);
			Vecd p_dissipation = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

				Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
				density_change_rate += u_jump * dW_ijV_j;
				p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;

				// Energy equation simplified by Haotian_Shi from SJTU
				energy_change_rate += 0.5 * (p_[index_i] + p_[index_j]) * u_jump * dW_ijV_j;
				// end 
			}
			drho_dt_[index_i] += density_change_rate * rho_[index_i];
			dE_dt_[index_i] += energy_change_rate / rho_[index_i];
			u_div_[index_i] = -density_change_rate;
			acc_[index_i] = p_dissipation / rho_[index_i];
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//