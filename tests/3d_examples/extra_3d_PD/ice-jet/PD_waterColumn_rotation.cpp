/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry

Real LR = 0.2;				  // liquid column radius in XZ
Real LH = 10 * LR;				  // liquid column height in Y
Real inner_circle_radius = LR;
int resolution(20);

Real plate_Y = 0.4;				/**< thickness of the ice plate. */
Real plate_R = 5 * plate_Y;	/**< width of the ice plate. */

Real resolution_ref = LR / 4;	  // particle spacing
Real BW = resolution_ref * 3; // boundary width

// for material properties of the fluid
Real rho0_f = 1000.0;
Real gravity_g = 9.81;
Real U_f = 100;
Real c_f = 10.0 * U_f;
Real mu_f = 0.0;
Real k_f = 0.0;
//----------------------------------------------------------------------
//	Material parameters of the elastic plate.
//----------------------------------------------------------------------
Real rho0_s = 900.0;	 /**< Reference density of plate. */
Real poisson = 0.3; /**< Poisson ratio. */
Real Youngs_modulus = 4.7e9;
//Real critical_stress = 8.436e8;
Real critical_stress = 5.0e6;
//	define the water block shape
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{		
		Vecd translation_column(0.0, 1.5 * LH, 0);
		add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), inner_circle_radius,
			0.5 * LH, resolution, translation_column);
	}
};

//	define the plate shape
class PlateShape : public ComplexShape
{
public:
	explicit PlateShape(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vecd offet_plate(0.0, LH - 0.5 * plate_Y - 0.5 * BW, 0.0);		
		add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), plate_R,
			0.5 * plate_Y, resolution, offet_plate);
	}
};
/**
 * Rotate water initial configuration
 */
class WaterInitialRotation
	: public fluid_dynamics::FluidInitialCondition, BaseSPHBodyRotation
{
public:
	explicit WaterInitialRotation(SPHBody& sph_body, Real theta_X, Real theta_Y, Real theta_Z)
		: fluid_dynamics::FluidInitialCondition(sph_body),
		BaseSPHBodyRotation(theta_X, theta_Y, theta_Z) {};

	void update(size_t index_i, Real dt)
	{
		pos_[index_i] = Q_ * pos_[index_i];
	};
};
/**
 * Rotate water initial configuration
 */
class PlateInitialRotation
	: public solid_dynamics::ElasticDynamicsInitialCondition, BaseSPHBodyRotation
{
public:
	explicit PlateInitialRotation(SPHBody& sph_body, Real theta_X, Real theta_Y, Real theta_Z)
		: solid_dynamics::ElasticDynamicsInitialCondition(sph_body),
		BaseSPHBodyRotation(theta_X, theta_Y, theta_Z) {};

	void update(size_t index_i, Real dt)
	{
		pos_[index_i] = Q_ * pos_[index_i];
	};
};
/**
 * application dependent initial condition
 */
class WaterInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit WaterInitialCondition(SPHBody& sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body) {};

	void update(size_t index_i, Real dt)
	{
		vel_[index_i][1] = -1.0 * U_f;
	};
};

// the main program with commandline options
int main(int ac, char *av[])
{	
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vecd(-plate_R - BW, -BW, -plate_R - BW), Vecd(plate_R + BW, 2.0 * LH + BW, plate_R + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);
	system.setRunParticleRelaxation(false);
	system.setReloadParticles(true);
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	size_t particle_num_f = water_block.getBaseParticles().total_real_particles_;
	
	size_t particle_num_w = 0;

	PDBody plate(system, makeShared<PlateShape>("PDBody"));
	plate.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	plate.defineParticlesAndMaterial<NosbPDParticles, HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? plate.generateParticles<ParticleGeneratorReload>(io_environment, plate.getName())
		: plate.generateParticles<ParticleGeneratorLattice>();
	plate.addBodyStateForRecording<Vecd>("NormalDirection");
	size_t particle_num_s = plate.getBaseParticles().total_real_particles_;

	size_t particle_num = particle_num_f + particle_num_w + particle_num_s;
	
	//ObserverBody fluid_observer(system, "FluidObserver");
	//fluid_observer.generateParticles<WaterObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	//ComplexRelation water_block_complex(water_block, {&wall_boundary});
	//ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
	ComplexRelation water_block_complex(water_block, RealBodyVector{ &plate });
	InnerRelation plate_inner_relation(plate);
	ContactRelation plate_water_contact_relation(plate, { &water_block });
	//ContactRelation plate_observer_contact_relation(plate_observer, { &plate });
	
	BodyStatesRecordingToVtp write_water_block_states(io_environment, system.real_bodies_);
	if (system.RunParticleRelaxation())
	{
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_column_particles(plate);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_plate_to_vtp(io_environment, plate);
		/** Write the particle reload files. */

		ReloadParticleIO write_particle_reload_files(io_environment, plate);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(plate_inner_relation);
		/**
		 * @brief 	Particle relaxation starts here.
		 */
		random_column_particles.parallel_exec(0.25);
		relaxation_step_inner.SurfaceBounding().parallel_exec();
		write_water_block_states.writeToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 2000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << setprecision(9) << "Relaxation steps for the column body N = " << ite_p << "\n";
				write_plate_to_vtp.writeToFile(ite_p);
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
	SimpleDynamics<WaterInitialRotation> jet_pos_rotation(water_block, 0.0, Pi / 9, 0.0);
	SimpleDynamics<WaterInitialCondition> jet_vel_initialization(water_block);
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, -gravity_g, 0.0));
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWallforPD> pressure_relaxation(water_block_complex);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWallforPD> density_relaxation(water_block_complex);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block,0.06);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------	
	SimpleDynamics<NormalDirectionFromBodyShape> plate_normal_direction(plate);
	InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidforPD> fluid_pressure_force_on_plate(plate_water_contact_relation);
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(plate);
	//----------------------------------------------------------------------
	//	Algorithms of Elastic dynamics.
	//----------------------------------------------------------------------
	SimpleDynamics<PlateInitialRotation> plate_pos_rotation(plate, 0.0, Pi / 4, 0.0);
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> plate_computing_time_step_size(plate, 0.024);
	/** calculate shape Matrix */
	InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> plate_shapeMatrix(plate_inner_relation);
	//stress relaxation for the beam by Hughes-Winget algorithm
	SimpleDynamics<solid_dynamics::NosbPDFirstStep> NosbPD_firstStep(plate);
	InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> NosbPD_secondStep(plate_inner_relation);
	InteractionDynamics<solid_dynamics::NosbPDThirdStep> NosbPD_thirdStep(plate_inner_relation);
	SimpleDynamics<solid_dynamics::NosbPDFourthStep> NosbPD_fourthStep(plate);
	//SimpleDynamics<solid_dynamics::NosbPDFourthStepWithADR> NosbPD_fourthStepADR(plate);
	//hourglass displacement mode control by LittleWood method
	InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl>
		hourglass_control(plate_inner_relation, plate.sph_adaptation_->getKernel());
	//Numerical Damping
	InteractionDynamics<solid_dynamics::PairNumericalDampingforPD>
		numerical_damping(plate_inner_relation, plate.sph_adaptation_->getKernel());
	//breaking bonds based on the MAX principal stress criteria
	InteractionDynamics<solid_dynamics::BondBreakByPrinStress> check_bondLive(plate_inner_relation, critical_stress);

	
	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> plate_update_normal(plate);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	
	/*RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, water_block, gravity_ptr);*/
	/*RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);*/
		//Log file
	std::string Logpath = io_environment.output_folder_ + "/SimLog.txt";
	if (fs::exists(Logpath))
	{
		fs::remove(Logpath);
	}
	std::ofstream log_file(Logpath.c_str(), ios::trunc);
	cout << "# PARAM SETING #" << "\n" << "\n"
		<< "	particle_num = " << particle_num << "\n"
		<< "	particle_spacing_ref = " << resolution_ref << "\n" << "\n"
		<< "<---- Fluid Domain ---->" << "\n" << "\n"
		<< "	particle_num_w = " << particle_num_w << "\n"
		<< "	particle_num_f = " << particle_num_f << "\n" << "\n"
		<< "	rho0_f = " << rho0_f << "\n"
		<< "	Characteristic velocity = " << U_f << "\n"
		<< "	Reference sound speed = " << c_f << "\n"
		<< "	Dynamics viscosity = " << mu_f << "\n"
		<< "	Thermal conduction rate = " << k_f << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
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
		<< "<---- Fluid Domain ---->" << "\n" << "\n"
		<< "	particle_num_w = " << particle_num_w << "\n"
		<< "	particle_num_f = " << particle_num_f << "\n" << "\n"
		<< "	rho0_f = " << rho0_f << "\n"
		<< "	Characteristic velocity = " << U_f << "\n"
		<< "	Reference sound speed = " << c_f << "\n"
		<< "	Dynamics viscosity = " << mu_f << "\n"
		<< "	Thermal conduction rate = " << k_f << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
		<< "<---- Struc Domain ---->" << "\n" << "\n"
		<< "	particle_num_s = " << particle_num_s << "\n" << "\n"
		<< "	rho0_s = " << rho0_s << "\n"
		<< "	Youngs_modulus = " << Youngs_modulus << "\n"
		<< "	poisson = " << poisson << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
		<< "# COMPUTATION START #" << "\n" << "\n";
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	jet_pos_rotation.parallel_exec();
	plate_pos_rotation.parallel_exec();
	jet_vel_initialization.parallel_exec();
	plate_normal_direction.parallel_exec();
	plate_shapeMatrix.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.RestartStep();
	size_t number_of_iterations_s = 0;
	int screen_output_interval = 10;
	Real end_time = 0.01;
	Real output_interval = end_time / 100.0;
	Real dt = 0.0;					// default acoustic time step sizes
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
	Real dt_s_0 = plate_computing_time_step_size.parallel_exec();
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_water_block_states.writeToFile(0);
	//write_water_mechanical_energy.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			/** Update normal direction at elastic body surface. */
			plate_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** Fluid relaxation and force computation. */
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_plate.parallel_exec();
				density_relaxation.parallel_exec(dt);
				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					//dt_s = plate_computing_time_step_size.parallel_exec();
					dt_s = dt_s_0;
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

					NosbPD_firstStep.parallel_exec(dt_s);
					check_bondLive.parallel_exec(dt_s);
					NosbPD_secondStep.parallel_exec(dt_s);
					hourglass_control.parallel_exec(dt_s);
					numerical_damping.parallel_exec(dt_s);
					NosbPD_thirdStep.parallel_exec(dt_s);
					
					NosbPD_fourthStep.parallel_exec(dt_s);					

					dt_s_sum += dt_s;
					/*if (number_of_iterations_s % screen_output_interval == 0)
					{
						std::cout << std::fixed << std::setprecision(9) 
							<< "		N_s=" << number_of_iterations_s 
							<< "	dt_s = " << dt_s << "\n";
					}
					number_of_iterations_s++;*/
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "	N=" << number_of_iterations << " Time: "
					<< GlobalStaticVariables::physical_time_ 
					<< "	Dt = " << Dt 
					<< "	dt = " << dt 
					<< "	dt_s = " << dt_s 
					<< "	dt / dt_s = " << dt / dt_s << "\n";

				log_file << std::fixed << std::setprecision(9) << "	N=" << number_of_iterations << " Time: "
					<< GlobalStaticVariables::physical_time_ 
					<< "	Dt = " << Dt 
					<< "	dt = " << dt 
					<< "	dt_s = " << dt_s
					<< "	dt / dt_s = " << dt / dt_s << "\n";
			}
			number_of_iterations++;

			water_block.updateCellLinkedListWithParticleSort(100);
			plate.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			plate_water_contact_relation.updateConfiguration();
			//fluid_observer_contact.updateConfiguration();
			//write_recorded_water_pressure.writeToFile(number_of_iterations);
		}

		//write_water_mechanical_energy.writeToFile(number_of_iterations);

		tick_count t2 = tick_count::now();
		write_water_block_states.writeToFile();
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
	cout << "\n" << "Total wall time for computation & output: " << tt2.seconds() << " seconds." << endl;
	log_file << "\n" << "Total wall time for computation & output: " << tt2.seconds() << " seconds." << endl;

	/*if (system.generate_regression_data_)
	{
		write_water_mechanical_energy.generateDataBase(1.0e-3);
		write_recorded_water_pressure.generateDataBase(1.0e-3);
	}
	else
	{
		write_water_mechanical_energy.newResultTest();
		write_recorded_water_pressure.newResultTest();
	}*/

	return 0;
}
