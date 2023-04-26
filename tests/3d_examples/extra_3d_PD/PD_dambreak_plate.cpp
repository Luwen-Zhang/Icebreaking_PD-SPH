/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry

Real DL = 0.584;			  // tank length
Real DH = 0.365;			  // tank height
Real DW = 0.05;				  // tank width
Real LL = 0.146;				  // liquid length
Real LH = 2 * LL;				  // liquid height
Real LW = 0.05;				  // liquid width
Real Rubber_length = 0.012;			/**< thickness of the rubber plate. */
Real Rubber_height = 20 / 3 * Rubber_length;			/**< Height of the rubber plate. */
Real Rubber_width = 0.048;			/**< width of the rubber plate. */
Real offset = (DW - Rubber_width) / 2;
Real resolution_ref = Rubber_length / 6;	  // particle spacing
Real BW = resolution_ref * 3; // boundary width

// for material properties of the fluid
Real rho0_f = 1000.0;
Real gravity_g = 9.81;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;
Real mu_f = 0.0;
Real k_f = 0.0;
//----------------------------------------------------------------------
//	Material parameters of the elastic plate.
//----------------------------------------------------------------------
Real rho0_s = 2500.0;	 /**< Reference density of plate. */
Real poisson = 0.0; /**< Poisson ratio. */
Real Youngs_modulus = 1e6;
//	define the water block shape
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
		Transformd translation_water(halfsize_water);
		add<TransformShape<GeometricShapeBox>>(Transformd(translation_water), halfsize_water);
	}
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
		Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
		Transformd translation_wall(halfsize_inner);
		add<TransformShape<GeometricShapeBox>>(Transformd(translation_wall), halfsize_outer);
		subtract<TransformShape<GeometricShapeBox>>(Transformd(translation_wall), halfsize_inner);
	}
};
//	define the plate shape
class PlateShape : public ComplexShape
{
public:
	explicit PlateShape(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vecd halfsize_plate(0.5 * Rubber_length, 0.5 * Rubber_height, 0.5 * Rubber_width);
		Vecd offet_plate(0.5 * Rubber_length + 2.0 * LL, 0.5 * Rubber_height, 0.5 * Rubber_width + offset);
		Transformd translation_plate(offet_plate);
		add<TransformShape<GeometricShapeBox>>(Transformd(translation_plate), halfsize_plate);
	}
};
//----------------------------------------------------------------------
// Create the plate constrain shape
//----------------------------------------------------------------------
Vecd halfsize_holder(0.5 * Rubber_length, resolution_ref * 3, 0.5 * Rubber_width);
Vecd translation_holder(0.5 * Rubber_length + 2.0 * LL, 0.5 * resolution_ref * 3, 0.5 * Rubber_width + offset);

//	define an observer particle generator
class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit WaterObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		// add observation points
		positions_.push_back(Vecd(DL, 0.01, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.1, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.2, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.24, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.252, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.266, 0.5 * DW));
	}
};

// the main program with commandline options
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	size_t particle_num_f = water_block.getBaseParticles().total_real_particles_;

	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");
	size_t particle_num_w = wall_boundary.getBaseParticles().total_real_particles_;

	SolidBody plate(system, makeShared<PlateShape>("PDBody"));
	plate.defineParticlesAndMaterial<NosbPDParticles, HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
	plate.generateParticles<ParticleGeneratorLattice>();
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
	ComplexRelation water_block_complex(water_block, RealBodyVector{ &wall_boundary, &plate });
	InnerRelation plate_inner_relation(plate);
	ContactRelation plate_water_contact_relation(plate, { &water_block });
	//ContactRelation plate_observer_contact_relation(plate_observer, { &plate });
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, -gravity_g, 0.0));
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWallforPD> pressure_relaxation(water_block_complex);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWallforPD> density_relaxation(water_block_complex);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> plate_normal_direction(plate);
	InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidforPD> fluid_pressure_force_on_plate(plate_water_contact_relation);
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(plate);
	//----------------------------------------------------------------------
	//	Algorithms of Elastic dynamics.
	//----------------------------------------------------------------------
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> plate_computing_time_step_size(plate, 0.24);
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

	/** Constrain the holder. */
	BodyRegionByParticle holder(plate,
		makeShared<TransformShape<GeometricShapeBox>>(Transformd(translation_holder), halfsize_holder, "Holder"));
	SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(holder);
	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> plate_update_normal(plate);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_water_block_states(io_environment, system.real_bodies_);
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
	wall_boundary_normal_direction.parallel_exec();
	plate_normal_direction.parallel_exec();
	plate_shapeMatrix.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 10;
	Real end_time = 1.0;
	Real output_interval = end_time / 100.0;
	Real dt = 0.0;					// default acoustic time step sizes
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
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
					dt_s = plate_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

					NosbPD_firstStep.parallel_exec(dt_s);
					NosbPD_secondStep.parallel_exec(dt_s);
					hourglass_control.parallel_exec(dt_s);
					numerical_damping.parallel_exec(dt_s);
					NosbPD_thirdStep.parallel_exec(dt_s);
					
					NosbPD_fourthStep.parallel_exec(dt_s);

					constraint_holder.parallel_exec(dt_s);

					dt_s_sum += dt_s;
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
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
				log_file << std::fixed << std::setprecision(9) << "	N=" << number_of_iterations << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
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
