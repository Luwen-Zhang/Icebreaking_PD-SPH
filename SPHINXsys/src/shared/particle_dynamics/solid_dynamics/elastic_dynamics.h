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
 * @file 	elastic_dynamics.h
 * @brief 	Here, we define the algorithm classes for elastic solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef ELASTIC_DYNAMICS_H
#define ELASTIC_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for elastic solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<ElasticSolidParticles> ElasticSolidDataSimple;
		typedef DataDelegateInner<ElasticSolidParticles> ElasticSolidDataInner;
		//Added by Haotian Shi from SJTU
		typedef DataDelegateSimple<NosbPDParticles> NosbPDSolidDataSimple;
		typedef DataDelegateInner<NosbPDParticles> NosbPDSolidDataInner;

		/**
		 * @class ElasticDynamicsInitialCondition
		 * @brief  set initial condition for a solid body with different material
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElasticDynamicsInitialCondition : public LocalDynamics, public ElasticSolidDataSimple
		{
		public:
			explicit ElasticDynamicsInitialCondition(SPHBody &sph_body);
			virtual ~ElasticDynamicsInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_, &vel_;
		};

		/**
		 * @class UpdateElasticNormalDirection
		 * @brief update particle normal directions for elastic solid
		 */
		class UpdateElasticNormalDirection : public LocalDynamics, public ElasticSolidDataSimple
		{
		protected:
			StdLargeVec<Vecd> &n_, &n0_;
			StdLargeVec<Matd> &F_;

		public:
			explicit UpdateElasticNormalDirection(SPHBody &sph_body);
			virtual ~UpdateElasticNormalDirection(){};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 * computing time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
									 public ElasticSolidDataSimple
		{
		protected:
			Real CFL_;
			StdLargeVec<Vecd> &vel_, &acc_, &acc_prior_;
			Real smoothing_length_, c0_;

		public:
			explicit AcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
			virtual ~AcousticTimeStepSize(){};

			Real reduce(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class DeformationGradientBySummation
		 * @brief computing deformation gradient tensor by summation
		 */
		class DeformationGradientBySummation : public LocalDynamics, public ElasticSolidDataInner
		{
		public:
			explicit DeformationGradientBySummation(BaseInnerRelation &inner_relation);
			virtual ~DeformationGradientBySummation(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd> &pos_;
			StdLargeVec<Matd> &B_, &F_;
		};

		/**
		 * @class BaseElasticIntegration
		 * @brief base class for elastic relaxation
		 */
		class BaseElasticIntegration : public LocalDynamics, public ElasticSolidDataInner
		{
		public:
			explicit BaseElasticIntegration(BaseInnerRelation &inner_relation);
			virtual ~BaseElasticIntegration(){};

		protected:
			StdLargeVec<Real> &rho_, &mass_;
			StdLargeVec<Vecd> &pos_, &vel_, &acc_;
			StdLargeVec<Matd> &B_, &F_, &dF_dt_;
		};

		/**
		 * @class BaseIntegration1stHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class BaseIntegration1stHalf : public BaseElasticIntegration
		{
		public:
			explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~BaseIntegration1stHalf(){};
			void update(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid &elastic_solid_;
			Real rho0_, inv_rho0_;
			StdLargeVec<Vecd> &acc_prior_;
			Real smoothing_length_;
		};

		/**
		 * @class Integration1stHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class Integration1stHalf : public BaseIntegration1stHalf
		{
		public:
			explicit Integration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~Integration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Matd> stress_PK1_B_;
			Real numerical_dissipation_factor_;
			Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd);
		};

		/**
		 * @class KirchhoffParticleIntegration1stHalf
		 */
		class KirchhoffParticleIntegration1stHalf : public Integration1stHalf
		{
		public:
			explicit KirchhoffParticleIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~KirchhoffParticleIntegration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class KirchhoffIntegration1stHalf
		 * @brief Decompose the stress into particle stress includes isotropic stress
		 * and the stress due to non-homogeneous material properties.
		 * The preliminary shear stress is introduced by particle pair to avoid
		 * spurious stress and deformation.
		 * Note that, for the shear stress term,
		 * due to the mismatch of the divergence contribution between
		 * the pair-wise second-order derivative Laplacian formulation
		 * and particle-wise first-order gradient formulation,
		 * a correction factor slight large than one is introduced.
		 * Note that, if you see time step size goes unusually small,
		 * it may be due to the determinate of deformation matrix become negative.
		 * In this case, you may need decrease CFL number when computing time-step size.
		 */
		class KirchhoffIntegration1stHalf : public BaseIntegration1stHalf
		{
		public:
			explicit KirchhoffIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~KirchhoffIntegration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> J_to_minus_2_over_dimension_;
			StdLargeVec<Matd> stress_on_particle_, inverse_F_T_;
			const Real correction_factor_ = 1.07;
		};

		/**
		 * @class Integration2ndHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the second step
		 */
		class Integration2ndHalf : public BaseElasticIntegration
		{
		public:
			explicit Integration2ndHalf(BaseInnerRelation &inner_relation)
				: BaseElasticIntegration(inner_relation){};
			virtual ~Integration2ndHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @Created by Haotian Shi from SJTU
		 * @class PDTimeStepInitialization
		 * @brief initialize a time step for a PDbody.
		 */
		class PDTimeStepInitialization
			: public BaseTimeStepInitialization,
			public NosbPDSolidDataSimple
		{
		protected:
			StdLargeVec<int>& particleLive_;
			StdLargeVec<Vecd>& pos_, & acc_prior_;

		public:
			PDTimeStepInitialization(SPHBody& sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
			virtual ~PDTimeStepInitialization() {};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * Created by Haotian Shi from SJTU
		 * @class NosbPDShapeMatrix
		 * @brief obtain the shape matrix in non-ordinary state based peridynamics		 
		 */
		class NosbPDShapeMatrix : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			explicit NosbPDShapeMatrix(BaseInnerRelation& inner_relation);
			virtual ~NosbPDShapeMatrix() {};
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Matd>& shape_K_, & shape_K_1_;
		};

		/**
		 * Created by Haotian Shi from SJTU
		 * @class NosbPDFirstStep
		 * @brief time marching into n+1/2 step
		 */
		class NosbPDFirstStep : public LocalDynamics, public NosbPDSolidDataSimple
		{
		public:
			explicit NosbPDFirstStep(SPHBody& sph_body);
			virtual ~NosbPDFirstStep() {};
			
			void update(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid& elastic_solid_;
			Real rho0_;
			StdLargeVec<Real>& rho_;
			StdLargeVec<Vecd>& pos_, & vel_, & acc_, & acc_old_, & acc_prior_;
			StdLargeVec<Matd>& F_;
		};

		/**
		 * Created by Haotian Shi from SJTU
		 * @class NosbPDSecondStep
		 * @brief calculate F_, F_half_, F_1_, F_1_half_, F_delta, and then calculate the Cauchy stress
		 */
		class NosbPDSecondStep : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			explicit NosbPDSecondStep(BaseInnerRelation& inner_relation);
			virtual ~NosbPDSecondStep() {};
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid& elastic_solid_;
			StdLargeVec<int> & particleLive_;
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& pos_, & vel_;
			StdLargeVec<Matd>& F_, & shape_K_1_;
			StdLargeVec<Matd>& N_, & N_deltaU_, & N_half_;
			StdLargeVec<Matd>& F_half_, & F_delta_, & F_1_, & F_1_half_;
			StdLargeVec<Matd>& PK1_, & T0_, & stress_;
		};
		/**
		* @Created by Haotian Shi from SJTU
		* @class NosbPDThirdStep
		* @brief calculate acc_ based on ordinary scheme, and then time march
		*/
		class NosbPDThirdStep : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			explicit NosbPDThirdStep(BaseInnerRelation& inner_relation);
			virtual ~NosbPDThirdStep() {};
			void interaction(size_t index_i, Real dt = 0.0);			

		protected:	
			ElasticSolid& elastic_solid_;
			Real rho0_;
			StdLargeVec<int>& particleLive_;
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& acc_;
			StdLargeVec<Matd>& T0_;
		};
		/**
		 * Created by Haotian Shi from SJTU
		 * @class NosbPDFourthStep
		 * @brief time marching into n+1 step
		 */
		class NosbPDFourthStep : public LocalDynamics, public NosbPDSolidDataSimple
		{
		public:
			explicit NosbPDFourthStep(SPHBody& sph_body);
			virtual ~NosbPDFourthStep() {};

			void update(size_t index_i, Real dt = 0.0);

		protected:			
			StdLargeVec<Vecd>& pos_, & vel_, & acc_;
		};
		/**
		* @Created by Haotian Shi from SJTU
		* @class LittleWoodHourGlassControl
		* @brief Hourglass displacement mode control in bond force format
		* Littlewood, D. J. (2011, January). ASME (Vol. 54945, pp. 567-576).
		*/
		class LittleWoodHourGlassControl : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			//hourglass constant coefficient, range: [1e-3, 1e2]
			Real Chg, Kbulk, horizon;
			Kernel* kernel_ptr;
			explicit LittleWoodHourGlassControl(BaseInnerRelation& inner_relation, Kernel* kernel);
			virtual ~LittleWoodHourGlassControl() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid& elastic_solid_;
			Real rho0_;
			StdLargeVec<int>& particleLive_;
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& pos_, & acc_;
			StdLargeVec<Matd>& F_;
		};
		/**
		* @Created by Haotian Shi from SJTU
		* @class MonaghanArtificialViscosityforPD
		* @brief shear viscosity proposed by Monaghan, 1989
		*/
		class MonaghanArtificialViscosityforPD : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			//hourglass constant coefficient, range: [1e-3, 1e2]
			Real alpha, beta;
			Kernel* kernel_ptr;
			explicit MonaghanArtificialViscosityforPD(BaseInnerRelation& inner_relation, Kernel* kernel);
			virtual ~MonaghanArtificialViscosityforPD() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid& elastic_solid_;
			Real rho0_, c0_;
			StdLargeVec<int>& particleLive_;
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& pos_, & vel_, & acc_;
		};
		/**
		* @Created by Haotian Shi from SJTU
		* @class PairNumericalDampingforPD
		* @brief shear viscosity proposed by Monaghan, 1989
		*/
		class PairNumericalDampingforPD : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			//hourglass constant coefficient, range: [1e-3, 1e2]
			Kernel* kernel_ptr;
			explicit PairNumericalDampingforPD(BaseInnerRelation& inner_relation, Kernel* kernel);
			virtual ~PairNumericalDampingforPD() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid& elastic_solid_;
			Real rho0_, c0_;
			Real numerical_dissipation_factor_;
			Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd);
			Real smoothing_length_;
			StdLargeVec<int>& particleLive_;
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& pos_, & vel_, & acc_;
			StdLargeVec<Matd>& F_;
		};
		/**
		* @Created by Haotian Shi from SJTU
		* @class NosbPDCheckBondLive
		* @brief basic class used to breaking bonds
		*/
		class NosbPDCheckBondLive : public LocalDynamics, public NosbPDSolidDataInner
		{
		public:
			explicit NosbPDCheckBondLive(BaseInnerRelation& inner_relation);
			virtual ~NosbPDCheckBondLive() {};
			
			virtual bool checkBondLive(Matd & EquivalentVar, Real & stretch_rate) = 0;

		protected:	
			Real critical_value_;
			StdLargeVec<int>& particleLive_;
			StdLargeVec<Vecd>& pos_, & vel_, & acc_;
		};
		/**
		* @Created by Haotian Shi from SJTU
		* @class BondBreakByPrinStress
		* @brief fracture criteria based on MAX principal stress
		*/
		class BondBreakByPrinStress : public NosbPDCheckBondLive
		{
		public:
			explicit BondBreakByPrinStress(BaseInnerRelation& inner_relation, Real& cr_value);
			virtual ~BondBreakByPrinStress() {};			

			virtual bool checkBondLive(Matd & stress_eq, Real & stretch_rate) override;

			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<int>& bond_;
			StdLargeVec<Real>& damage_;
			StdLargeVec<Matd>& stress_;
		};
		/**
		 * @Created by Haotian Shi from SJTU
		 * @class ADRFirstStep
		 * @brief computing SUM(U^T * K * U) 
		 */
		class ADRFirstStep : public LocalDynamicsReduce<Real, ReduceSum<Real>>,
			public NosbPDSolidDataSimple
		{
		public:		

			explicit ADRFirstStep(SPHBody& sph_body);
			virtual ~ADRFirstStep() {};

			Real reduce(size_t index_i, Real dt = 0.0);

		protected:			
			StdLargeVec<Vecd>& pos0_, & pos_, & vel_;
			StdLargeVec<Vecd>& acc_, & acc_old_, & acc_prior_;
		};
		/**
		 * @Created by Haotian Shi from SJTU
		 * @class ADRSecondStep
		 * @brief computing SUM(U^T * U)
		 */
		class ADRSecondStep : public LocalDynamicsReduce<Real, ReduceSum<Real>>,
			public NosbPDSolidDataSimple
		{
		public:

			explicit ADRSecondStep(SPHBody& sph_body);
			virtual ~ADRSecondStep() {};

			Real reduce(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd>& pos0_, & pos_;
		};
		/**
		 * Created by Haotian Shi from SJTU
		 * @class NosbPDFourthStepWithADR
		 * @brief ADR time marching into n+1 step
		 */
		class NosbPDFourthStepWithADR : public LocalDynamics, public NosbPDSolidDataSimple
		{
		public:
			explicit NosbPDFourthStepWithADR(SPHBody& sph_body);
			virtual ~NosbPDFourthStepWithADR() {};

			void update(size_t index_i, Real dt = 0.0);
			void getADRcn(Real& ADRcn);

		protected:
			Real ADR_cn_;
			StdLargeVec<Vecd>& pos_, & vel_, & acc_;
		};
	}
}
#endif // ELASTIC_DYNAMICS_H