/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry

Real LR = 0.002;				  // liquid column radius in XZ
Real LH = 0.04428;				  // liquid column height in Y

int resolution(20);

Real plate_Y = 0.01;				/**< thickness of the ice plate. */
Real plate_R = 0.02;	/**< width of the ice plate. */

Real resolution_ref = 5e-4;	  // particle spacing
Real BW = resolution_ref * 2; // boundary width
Real inner_circle_radius = LR;
// for material properties of the fluid
Real rho0_f = 1000.0;
//Real gravity_g = 9.81;
Real gravity_g = 0.0;
Real U_f = 62.52;
//Real c_f = 10.0 * U_f;
Real c_f = SMIN(10.0 * U_f, 1580.0);
Real mu_f = 0.0;
Real k_f = 0.0;
//----------------------------------------------------------------------
//	Material parameters of the elastic plate.
//----------------------------------------------------------------------
Real rho0_s = 900.0;	 /**< Reference density of plate. */
Real poisson = 0.33; /**< Poisson ratio. */
Real Youngs_modulus = 6.2e9;
Real sigma_t0 = 2.2e6;
//Real sigma_c0 = 9.4e6 * 3;
Real epsilon_rate = U_f / plate_Y;
Real sigma_c0 = 10.976e6 * pow(epsilon_rate, 0.093783);
//Hardening parameters for J2
Real sigma_Y0 = sigma_t0;
//Hardening parameters for DP
//Real alpha = (sigma_c0 - sigma_t0) / (sigma_c0 + sigma_t0) / sqrt(3);
//Real flow_stress = (2.0 * sigma_c0 * sigma_t0) / (sigma_c0 + sigma_t0) / sqrt(3);
Real alpha = 0.4423 / sqrt(3);
Real flow_stress = 1.1e6 / sqrt(3);
//Shared hardening parameters
Real isohardening_modulus_H = 6.89e6;
Real kinhardening_modulus_H = 0.0;
//Damage parameters
//Real max_tension_stress = SMAX(sigma_t0, flow_stress / alpha / 3);
Real max_tension_stress = 0.7e6;
Real max_pressure = sigma_c0 ;
Real max_shear_stress = 7.0e6;
Real max_plastic_strain = 0.0035;
Real max_stretch = 3.876e-4;

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
		Vecd offet_plate(0.0, LH - 0.5 * plate_Y - 2.0 * resolution_ref, 0.0);
		add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), plate_R,
			0.5 * plate_Y, resolution, offet_plate);
	}
};
//	define the plate constrain shape
class PlateConstrainShape : public ComplexShape
{
public:
	explicit PlateConstrainShape(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vecd offet_plate(0.0, 10 * LR - 0.5 * plate_Y - 0.5 * BW, 0.0);
		add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), plate_R,
			0.6 * plate_Y, resolution, offet_plate);
		subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), plate_R - BW,
			0.65 * plate_Y, resolution, offet_plate);
	}
};
class BoundaryCondition
	: public solid_dynamics::LoadBodyPartConstraint
{
protected:
	const Real coeff;
public:
	explicit BoundaryCondition(BodyPartByParticle& body_part)
		: solid_dynamics::LoadBodyPartConstraint(body_part), coeff(0.001) {};
	

	void update(size_t index_i, Real dt)
	{
		//vel_[index_i][1] = 0.0;
		//acc_[index_i][1] = 0.0;
		acc_[index_i] = coeff * acc_[index_i];
	};
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
	//plate.defineAdaptationRatios(1.5075, 0.8);
	//SolidBody plate(system, makeShared<PlateShape>("SPHBody"));
	//plate.defineAdaptationRatios(1.15, 0.8);
	plate.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	//plate.defineParticlesAndMaterial<NosbPDParticles, HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
	plate.defineParticlesAndMaterial<NosbPDPlasticParticles, DruckerPragerPlasticityforPD>(rho0_s, Youngs_modulus, poisson,
		flow_stress, alpha, isohardening_modulus_H, kinhardening_modulus_H);
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
		while (ite_p < 1000)
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
		write_water_block_states.writeToFile(1200);
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
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
	Real h = water_block.sph_adaptation_->getKernel()->CutOffRadius();
	pressure_relaxation.setcoeffacousticdamper(rho0_f, c_f, h);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block,0.06);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------	
	SimpleDynamics<NormalDirectionFromBodyShape> plate_normal_direction(plate);
	InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidforPD> fluid_pressure_force_on_plate(plate_water_contact_relation);
	fluid_pressure_force_on_plate.setcoeffacousticdamper(rho0_f, c_f, h);
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
	//InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> NosbPD_secondStep(plate_inner_relation);
	InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic> NosbPD_secondStepPlastic(plate_inner_relation);
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
	//InteractionDynamics<solid_dynamics::BondBreakBySigma1andSigma3> check_bondLive(plate_inner_relation, max_tension_stress, max_shear_stress);
	InteractionDynamics<solid_dynamics::BondBreakByPrinStress> check_bondLive(plate_inner_relation, max_tension_stress);
	//InteractionDynamics<solid_dynamics::BondBreakBySigma1andNorm1> check_bondLive(plate_inner_relation, max_tension_stress, max_pressure);		
	//InteractionDynamics<solid_dynamics::BondBreakByPlasticStrain> check_bondLive(plate_inner_relation, max_plastic_strain);

	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> plate_update_normal(plate);
	/** Constrain the boundary */
	BodyRegionByParticle boundary(plate, makeShared<PlateConstrainShape>("Boundary"));
	size_t boundary_num = boundary.SizeOfLoopRange();
	SimpleDynamics<BoundaryCondition> constraint_boundary(boundary);
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
		<< "<---- ICE PROPERTY ---->" << "\n" << "\n"
		<< "	alpha = " << alpha * sqrt(3) << "\n"
		<< "	flow_stress = " << flow_stress * sqrt(3) << "\n"
		<< "	P0 = " << flow_stress / alpha / 3 << "\n"
		<< "	max_tension_stress = " << max_tension_stress << "\n"
		<< "	max_shear_stress = " << max_shear_stress << "\n" << "\n"
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
		<< "<---- ICE PROPERTY ---->" << "\n" << "\n"
		<< "	alpha = " << alpha * sqrt(3) << "\n"
		<< "	flow_stress = " << flow_stress * sqrt(3) << "\n"
		<< "	P0 = " << flow_stress / alpha / 3 << "\n"
		<< "	max_tension_stress = " << max_tension_stress << "\n"
		<< "	max_shear_stress = " << max_shear_stress << "\n" << "\n"
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
	int screen_output_interval = 1;
	int screen_output_interval_s = 3;
	Real end_time = 0.4e-4;
	Real output_interval = end_time / 200.0;
	Real dt = 0.0;					// default acoustic time step sizes
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
	//Real dt_s_0 = plate_computing_time_step_size.parallel_exec();
	Real dt_s_0 = 4e-8;
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
					//check_bondLive.parallel_exec(dt_s);
					//NosbPD_secondStep.parallel_exec(dt_s);
					NosbPD_secondStepPlastic.parallel_exec(dt_s);

					hourglass_control.parallel_exec(dt_s);
					numerical_damping.parallel_exec(dt_s);

					NosbPD_thirdStep.parallel_exec(dt_s);
					//constraint_boundary.parallel_exec(dt_s);
					NosbPD_fourthStep.parallel_exec(dt_s);					

					dt_s_sum += dt_s;
					if (number_of_iterations_s % screen_output_interval_s == 0)
					{
						/*std::cout << std::fixed << std::setprecision(9) 
							<< "		N_s=" << number_of_iterations_s 
							<< "	dt_s = " << dt_s << "\n";*/
						tick_count t2 = tick_count::now();
						write_water_block_states.writeToFile();
						time_file << std::fixed << std::setprecision(9) << GlobalStaticVariables::physical_time_ << "\n";
						tick_count t3 = tick_count::now();
						interval += t3 - t2;
					}
					number_of_iterations_s++;
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

		/*tick_count t2 = tick_count::now();
		write_water_block_states.writeToFile();
		time_file << std::fixed << std::setprecision(9) << GlobalStaticVariables::physical_time_ << "\n";
		tick_count t3 = tick_count::now();
		interval += t3 - t2;*/
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
