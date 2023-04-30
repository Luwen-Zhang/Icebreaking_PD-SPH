/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 500.0;						/**< Tank length. */
Real DH = 200.1;						/**< Tank height. */
Real Dam_L = 100.0;						/**< Water block width. */
Real Dam_H = 140.0;						/**< Water block height. */
Real Gate_width = 5.0;					/**< Width of the gate. */
Real Base_bottom_position = 79.0;		/**< Position of gate base. (In Y direction) */
Real resolution_ref = Gate_width / 8.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0;			/**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Real dp_s = 0.5 * resolution_ref;
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / dp_s) * dp_s);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/** Define the corner points of the water block geometry. */
Vec2d DamP_lb(DL - Dam_L, 0.0);	  /**< Left bottom. */
Vec2d DamP_lt(DL - Dam_L, Dam_H); /**< Left top. */
Vec2d DamP_rt(DL, Dam_H);		  /**< Right top. */
Vec2d DamP_rb(DL, 0.0);			  /**< Right bottom. */
/** Define the corner points of the gate geometry. */
Vec2d GateP_lb(DL - Dam_L - Gate_width, 0.0);		 /**< Left bottom. */
Vec2d GateP_lt(DL - Dam_L - Gate_width, Dam_H + BW); /**< Left top. */
Vec2d GateP_rt(DL - Dam_L, Dam_H + BW);				 /**< Right top. */
Vec2d GateP_rb(DL - Dam_L, 0.0);					 /**< Right bottom. */
/** Define the corner points of the gate constrain. */
Vec2d ConstrainP_lb(DL - Dam_L - Gate_width, Base_bottom_position); /**< Left bottom. */
Vec2d ConstrainP_lt(DL - Dam_L - Gate_width, Dam_H + BW);			/**< Left top. */
Vec2d ConstrainP_rt(DL - Dam_L, Dam_H + BW);						/**< Right top. */
Vec2d ConstrainP_rb(DL - Dam_L, Base_bottom_position);				/**< Right bottom. */
// observer location
StdVec<Vecd> observation_location = {GateP_lb};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						   /**< Reference density of fluid. */
Real gravity_g = 9.8e-3;				   /**< Value of gravity. */
Real U_f = 1.0;							   /**< Characteristic velocity. */
Real c_f = 20.0 * sqrt(140.0 * gravity_g); /**< Reference sound speed. */
Real mu_f = 0.0;
Real k_f = 0.0;
//----------------------------------------------------------------------
//	Material parameters of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 2.0;	 /**< Reference density of gate. */
Real poisson = 0.3; /**< Poisson ratio. */
Real Ae = 7.8e3;	 /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_block_shape;
		water_block_shape.push_back(DamP_lb);
		water_block_shape.push_back(DamP_lt);
		water_block_shape.push_back(DamP_rt);
		water_block_shape.push_back(DamP_rb);
		water_block_shape.push_back(DamP_lb);
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Wall cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	create a gate shape
//----------------------------------------------------------------------
MultiPolygon createGateShape()
{
	std::vector<Vecd> gate_shape;
	gate_shape.push_back(GateP_lb);
	gate_shape.push_back(GateP_lt);
	gate_shape.push_back(GateP_rt);
	gate_shape.push_back(GateP_rb);
	gate_shape.push_back(GateP_lb);

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(gate_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
// Create the gate constrain shape
//----------------------------------------------------------------------
MultiPolygon createGateConstrainShape()
{
	// geometry
	std::vector<Vecd> gate_constraint_shape;
	gate_constraint_shape.push_back(ConstrainP_lb);
	gate_constraint_shape.push_back(ConstrainP_lt);
	gate_constraint_shape.push_back(ConstrainP_rt);
	gate_constraint_shape.push_back(ConstrainP_rb);
	gate_constraint_shape.push_back(ConstrainP_lb);

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(gate_constraint_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBlock"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	size_t particle_num_f = water_block.getBaseParticles().total_real_particles_;

	SolidBody wall_boundary(system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	size_t particle_num_w = wall_boundary.getBaseParticles().total_real_particles_;

	PDBody gate(system, makeShared<MultiPolygonShape>(createGateShape(), "PDBody"));
	gate.defineParticlesAndMaterial<NosbPDParticles, HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
	gate.generateParticles<ParticleGeneratorLattice>();
	size_t particle_num_s = gate.getBaseParticles().total_real_particles_;

	size_t particle_num = particle_num_f + particle_num_w + particle_num_s;
	ObserverBody gate_observer(system, "Observer");
	gate_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_block_complex_relation(water_block, RealBodyVector{&wall_boundary, &gate});
	InnerRelation gate_inner_relation(gate);
	ContactRelation gate_water_contact_relation(gate, {&water_block});
	ContactRelation gate_observer_contact_relation(gate_observer, {&gate});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//	Algorithms of fluid dynamics.
	//----------------------------------------------------------------------
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWallforPD> pressure_relaxation(water_block_complex_relation);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWallforPD> density_relaxation(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex_relation);
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	SimpleDynamics<OffsetInitialPosition> gate_offset_position(gate, offset);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> gate_normal_direction(gate);	
	InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidforPD> fluid_pressure_force_on_gate(gate_water_contact_relation);
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
	//----------------------------------------------------------------------
	//	Algorithms of Elastic dynamics.
	//----------------------------------------------------------------------
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> gate_computing_time_step_size(gate, 0.24);
	/** calculate shape Matrix */
	InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> gate_shapeMatrix(gate_inner_relation);
	//stress relaxation for the beam by Hughes-Winget algorithm
	SimpleDynamics<solid_dynamics::NosbPDFirstStep> NosbPD_firstStep(gate);
	InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> NosbPD_secondStep(gate_inner_relation);
	InteractionDynamics<solid_dynamics::NosbPDThirdStep> NosbPD_thirdStep(gate_inner_relation);
	//SimpleDynamics<solid_dynamics::NosbPDFourthStep> NosbPD_fourthStep(gate);
	 SimpleDynamics<solid_dynamics::NosbPDFourthStepWithADR> NosbPD_fourthStepADR(gate);
	//hourglass displacement mode control by LittleWood method
	InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass_control(gate_inner_relation, gate.sph_adaptation_->getKernel());
	//Numerical Damping
	InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> numerical_damping(gate_inner_relation, gate.sph_adaptation_->getKernel());
	
	// ADR_cn calculation
	ReduceDynamics<solid_dynamics::ADRFirstStep> computing_cn1(gate);
	ReduceDynamics<solid_dynamics::ADRSecondStep> computing_cn2(gate);	

	BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
	SimpleDynamics<solid_dynamics::FixBodyPartConstraint> gate_constraint(gate_constraint_part);
	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> gate_update_normal(gate);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToPlt write_real_body_states_to_plt(io_environment, system.real_bodies_);
	BodyStatesRecordingToVtp write_real_body_states_to_vtp(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_beam_tip_displacement("Position", io_environment, gate_observer_contact_relation);
	//TODO: observing position is not as good observing displacement. 
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
	gate_offset_position.parallel_exec();
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_normal_direction.parallel_exec();
	gate_normal_direction.parallel_exec();
	gate_shapeMatrix.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int number_of_iterations = 0;
	int screen_output_interval = 10;
	Real end_time = 100.0;
	Real output_interval = end_time / 200.0;
	Real dt = 0.0;					/**< Default acoustic time step sizes. */
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */

	Real cn1 = 0.0;
	Real cn2 = 0.0;
	Real ADR_cn = 0.0;

	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states_to_vtp.writeToFile();
	write_beam_tip_displacement.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			/** Acceleration due to viscous force and gravity. */
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			/** Update normal direction at elastic body surface. */
			gate_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** Fluid relaxation and force computation. */
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_gate.parallel_exec();
				density_relaxation.parallel_exec(dt);
				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = gate_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

					NosbPD_firstStep.parallel_exec(dt_s);					
					NosbPD_secondStep.parallel_exec(dt_s);
					hourglass_control.parallel_exec(dt_s);
					numerical_damping.parallel_exec(dt_s);
					NosbPD_thirdStep.parallel_exec(dt_s);
					/*cn1 = SMAX(TinyReal, computing_cn1.parallel_exec(dt));
					cn2 = computing_cn2.parallel_exec(dt);
					if (cn2 > TinyReal) {
						ADR_cn = 2.0 * sqrt(cn1 / cn2);
					}
					else {
						ADR_cn = 0.0;
					}
					ADR_cn = 0.001 * ADR_cn;
					NosbPD_fourthStepADR.getADRcn(ADR_cn);*/
					NosbPD_fourthStepADR.parallel_exec(dt_s);

					gate_constraint.parallel_exec(dt_s);

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
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedListWithParticleSort(100);
			gate.updateCellLinkedList();
			water_block_complex_relation.updateConfiguration();
			gate_water_contact_relation.updateConfiguration();
			/** Output the observed data. */
			write_beam_tip_displacement.writeToFile(number_of_iterations);
		}
		tick_count t2 = tick_count::now();
		write_real_body_states_to_vtp.writeToFile();
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

	//write_beam_tip_displacement.newResultTest();

	return 0;
}
