/*-----------------------------------------------------------------------------*
 *                    SPHinXsys: 3D example J2 Necking for PD                  *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry

Real strip_Y = 0.018;	/**< 0.5 length of the strip. */
Real strip_X = 0.01;	/**< 0.5 width of the strip. */
Real strip_Z = 0.0005;	/**< 0.5 thickness of the strip. */
int resolution(20);
Real resolution_ref = strip_Z / 2;	  // particle spacing
Real BW = resolution_ref * 4; // boundary width

Real gravity_g = 9.81;
Real Load_vel = 10;
//----------------------------------------------------------------------
//	Material parameters of the elastoplastic strip.
//----------------------------------------------------------------------
Real rho0_s = 7850.0;	 /**< Reference density of strip. */
Real poisson = 0.29; /**< Poisson ratio. */
Real Youngs_modulus = 206.9e9;
//Hardening parameters
Real sigma_Y0 = 0.45e9;
Real sigma_Y_infinite = 0.715e9;
Real index_kesi = 16.93;
Real isohardening_modulus_H = 0.12924e9;
//Damage parameters
Real max_tension_stress = 2.2e6;
Real max_shear_stress = 1.7e6;

Vecd halfsize_strip(strip_X, strip_Y, strip_Z);
Vecd translation_strip(0.0, 0.0, 0.0);
Vecd halfsize_holder(strip_X, 0.5 * BW, strip_Z);
Vecd translation_holder(0.0, strip_Y - 0.5 * BW, 0.0);

//	define the strip shape
class StripShape : public ComplexShape
{
public:
	explicit StripShape(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transformd(translation_strip), halfsize_strip);
	}
};
//	define the strip constrain shape
class StripConstrainShape : public ComplexShape
{
public:
	explicit StripConstrainShape(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vecd bottom_holder = -translation_holder;
		add<TransformShape<GeometricShapeBox>>(Transformd(bottom_holder), halfsize_holder);
	}
};
//	define the strip loading shape
class StripLoadingShape : public ComplexShape
{
public:
	explicit StripLoadingShape(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vecd up_holder = translation_holder;
		add<TransformShape<GeometricShapeBox>>(Transformd(up_holder), halfsize_holder);
	}
};
class UpStripLoadingCondition
	: public solid_dynamics::LoadBodyPartConstraint
{
protected:
	Real V_;
public:
	explicit UpStripLoadingCondition(BodyPartByParticle& body_part)
		: solid_dynamics::LoadBodyPartConstraint(body_part),V_(0.0) {};

	inline Real getVel(const Real& time)
	{
		//V_ = 8.75e4 * time;
		//V_ = 3.28e8 * time * time;
		V_ = 2.0 * Pi * sin(Pi * 1e4 * time);
		return V_;
	}

	void update(size_t index_i, Real dt)
	{		
		vel_[index_i][1] = V_;
	};
};
class BottomStripLoadingCondition
	: public solid_dynamics::LoadBodyPartConstraint
{
protected:
	Real V_;
public:
	explicit BottomStripLoadingCondition(BodyPartByParticle& body_part)
		: solid_dynamics::LoadBodyPartConstraint(body_part), V_(0.0) {};

	inline Real getVel(const Real& time)
	{
		//V_ = -8.75e4 * time;
		//V_ = -3.28e8 * time * time;
		V_ = -2.0 * Pi * sin(Pi * 1e4 * time);
		return V_;
	}

	void update(size_t index_i, Real dt)
	{
		vel_[index_i][1] = V_;
	};
};
// the main program with commandline options
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vecd(-strip_X - BW, -strip_Y - BW, -strip_Z - BW), Vecd(strip_X + BW, strip_Y + BW, strip_Z + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);
	system.setRunParticleRelaxation(false);
	system.setReloadParticles(true);
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);

	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	PDBody strip(system, makeShared<StripShape>("PDBody"));
	//SolidBody strip(system, makeShared<StripShape>("SPHBody"));
	strip.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	strip.defineParticlesAndMaterial<NosbPDPlasticParticles, J2PlasticityforPD>(rho0_s, Youngs_modulus, poisson, 
		sigma_Y0, isohardening_modulus_H, 0.0);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? strip.generateParticles<ParticleGeneratorReload>(io_environment, strip.getName())
		: strip.generateParticles<ParticleGeneratorLattice>();
	strip.addBodyStateForRecording<Vecd>("NormalDirection");
	size_t particle_num_s = strip.getBaseParticles().total_real_particles_;

	size_t particle_num = particle_num_s;

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation strip_inner_relation(strip);

	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	if (system.RunParticleRelaxation())
	{
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_column_particles(strip);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_strip_to_vtp(io_environment, strip);
		/** Write the particle reload files. */

		ReloadParticleIO write_particle_reload_files(io_environment, strip);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(strip_inner_relation);
		/**
		 * @brief 	Particle relaxation starts here.
		 */
		random_column_particles.parallel_exec(0.25);
		relaxation_step_inner.SurfaceBounding().parallel_exec();
		write_states.writeToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << setprecision(9) << "Relaxation steps for the column body N = " << ite_p << "\n";
				write_strip_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of cylinder body finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0.0);
		return 0;
	}

	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//	Algorithms of Elastic dynamics.
	//----------------------------------------------------------------------

	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> strip_computing_time_step_size(strip, 0.024);
	/** calculate shape Matrix */
	InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> strip_shapeMatrix(strip_inner_relation);
	//stress relaxation for the beam by Hughes-Winget algorithm
	SimpleDynamics<solid_dynamics::NosbPDFirstStep> NosbPD_firstStep(strip);
	//InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> NosbPD_secondStep(strip_inner_relation);
	InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic> NosbPD_secondStepPlastic(strip_inner_relation);
	InteractionDynamics<solid_dynamics::NosbPDThirdStep> NosbPD_thirdStep(strip_inner_relation);
	SimpleDynamics<solid_dynamics::NosbPDFourthStep> NosbPD_fourthStep(strip);
	// ADR_cn calculation
	ReduceDynamics<solid_dynamics::ADRFirstStep> computing_cn1(strip);
	ReduceDynamics<solid_dynamics::ADRSecondStep> computing_cn2(strip);
	SimpleDynamics<solid_dynamics::NosbPDFourthStepWithADR> NosbPD_fourthStepADR(strip);
	//hourglass displacement mode control by LittleWood method
	InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl>
		hourglass_control(strip_inner_relation, strip.sph_adaptation_->getKernel());
	//Numerical Damping
	InteractionDynamics<solid_dynamics::PairNumericalDampingforPD>
		numerical_damping(strip_inner_relation, strip.sph_adaptation_->getKernel());
	//breaking bonds based on the MAX principal stress criteria
	//InteractionDynamics<solid_dynamics::BondBreakBySigma1andSigma3> check_bondLive(strip_inner_relation, max_tension_stress, max_shear_stress);
	InteractionDynamics<solid_dynamics::BondBreakByPrinStress> check_bondLive(strip_inner_relation, max_tension_stress);

	/** Constrain the holder. */
	BodyRegionByParticle holder(strip, makeShared<StripConstrainShape>("Holder"));
	SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(holder);
	/** Constrain the up loader. */
	BodyRegionByParticle loader(strip, makeShared<StripLoadingShape>("Loader"));
	SimpleDynamics<UpStripLoadingCondition> constraint_loader_up(loader);
	/** Constrain the bottom loader. */	
	SimpleDynamics<BottomStripLoadingCondition> constraint_loader_bottom(holder);

	std::string Logpath = io_environment.output_folder_ + "/SimLog.txt";
	if (fs::exists(Logpath))
	{
		fs::remove(Logpath);
	}
	std::ofstream log_file(Logpath.c_str(), ios::trunc);
	std::cout << "# PARAM SETING #" << "\n" << "\n"
		<< "	particle_num = " << particle_num << "\n"
		<< "	particle_spacing_ref = " << resolution_ref << "\n" << "\n"		
		<< "<---- Struc Domain ---->" << "\n" << "\n"
		<< "	particle_num_s = " << particle_num_s << "\n" << "\n"
		<< "	rho0_s = " << rho0_s << "\n"
		<< "	Youngs_modulus = " << Youngs_modulus << "\n"
		<< "	poisson = " << poisson << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
		<< "# COMPUTATION START # " << "\n" << "\n";
	log_file << "#PARAM SETING#" << "\n" << "\n"
		<< "	particle_num = " << particle_num << "\n" << "\n"
		<< "	particle_spacing_ref = " << resolution_ref << "\n" << "\n"		
		<< "<---- Struc Domain ---->" << "\n" << "\n"
		<< "	particle_num_s = " << particle_num_s << "\n" << "\n"
		<< "	rho0_s = " << rho0_s << "\n"
		<< "	Youngs_modulus = " << Youngs_modulus << "\n"
		<< "	poisson = " << poisson << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
		<< "# COMPUTATION START #" << "\n" << "\n";

	//Time file
	std::string Timepath = io_environment.output_folder_ + "/LogTime.txt";
	if (fs::exists(Timepath))
	{
		fs::remove(Timepath);
	}
	std::ofstream time_file(Timepath.c_str(), ios::trunc);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();	
	strip_shapeMatrix.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.RestartStep();
	size_t number_of_iterations_s = 0;
	int screen_output_interval = 1;
	int ite = 0;
	Real end_time = 4e-4;
	Real output_interval = end_time / 100.0;
	Real dt = 0.0;					// default acoustic time step sizes
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
	//Real dt_s_0 = plate_computing_time_step_size.parallel_exec();
	Real dt_s_0 = 4e-8;
	Real cn1 = 0.0;
	Real cn2 = 0.0;
	Real ADR_cn = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile(0);
	//write_water_mechanical_energy.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			ite++;
			dt = dt_s_0;
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			if (ite % 100 == 0) {
				std::cout << "	N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
				log_file << "	N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			//initialize_a_solid_step.parallel_exec(dt);
			//Loading
			constraint_loader_up.getVel(GlobalStaticVariables::physical_time_);
			constraint_loader_bottom.getVel(GlobalStaticVariables::physical_time_);
			constraint_loader_up.parallel_exec(dt);
			constraint_loader_bottom.parallel_exec(dt);
			//Spatial Integration
			NosbPD_firstStep.parallel_exec(dt);
			NosbPD_secondStepPlastic.parallel_exec(dt);
			//Spatial Numerical dissipation
			hourglass_control.parallel_exec(dt);
			numerical_damping.parallel_exec(dt);
			
			NosbPD_thirdStep.parallel_exec(dt);
			//Time Numerical dissipation
			cn1 = SMAX(TinyReal, computing_cn1.parallel_exec(dt));
			cn2 = computing_cn2.parallel_exec(dt);
			if (cn2 > TinyReal) {
				ADR_cn = 2.0 * sqrt(cn1 / cn2);
			}
			else {
				ADR_cn = 0.0;
			}
			ADR_cn = 0.01 * ADR_cn;
			NosbPD_fourthStepADR.getADRcn(ADR_cn);
			//Time Integration
			NosbPD_fourthStepADR.parallel_exec(dt);
			//NosbPD_fourthStep.parallel_exec(dt);

			//constraint_holder.parallel_exec(dt);
			
		}		
		

		//write_water_mechanical_energy.writeToFile(number_of_iterations);

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		time_file << std::fixed << std::setprecision(9) << GlobalStaticVariables::physical_time_ << "\n";
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tick_count::interval_t tt2;
	tt = t4 - t1 - interval;
	tt2 = t4 - t1;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	log_file << "\n" << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
	std::cout << "\n" << "Total wall time for computation & output: " << tt2.seconds() << " seconds." << endl;
	log_file << "\n" << "Total wall time for computation & output: " << tt2.seconds() << " seconds." << endl;

	return 0;
}