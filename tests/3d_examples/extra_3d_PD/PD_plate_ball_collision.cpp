/**
 * @file passive_plate.cpp
 * @brief This is the first example of plate 
 * @author Chi Zhang and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 6.0;//X:length 
Real PH = 6.0;//Y:width
Real PW = 1.0;//Z:thickness
Real ZH = 6.0;
Real resolution_ref = PW / 10.0; /**< Initial particle spacing. */
Real BW = resolution_ref * 4;	 /**< Boundary width. */
Vecd halfsize_plate(0.5 * PL, 0.5 * PH, 0.5 * PW);
Vecd translation_plate(0.5 * PL, 0.5 * PH, 0.5 * PW + ZH);
Real ball_radius = PL / 8;
Vecd ball_center(0.5 * PL, 0.5 * PH, ZH - ball_radius - 0.1 * PW);

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(- BW, -BW, -BW),
								 Vecd(PL + BW, PH + BW, PW + ZH + BW));

/** For material properties of the solid. */
Real rho0_s = 1100.0;
Real poisson = 0.45;
Real Youngs_modulus = 1.7e7;
Real critical_stress = 1.0e6;


Real a = Youngs_modulus / (2.0 * (1.0 + poisson));
Real a_f = 0.0 * a;
Real a0[4] = {a, a_f, 0.0, 0.0};
Real b0[4] = {1.0, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real bulk_modulus = Youngs_modulus / 3.0 / (1.0 - 2.0 * poisson);

Real gravity_g = 9.81;

/** Define the plate shape. */
class Plate : public ComplexShape
{
public:
	explicit Plate(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transformd(translation_plate), halfsize_plate);		
	}
};
/**
 * application dependent initial condition 
 */
class PlateInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit PlateInitialCondition(SPHBody &sph_body)
		: solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

	void update(size_t index_i, Real dt)
	{
		Real coff = 1.0;
		if (pos_[index_i][0] > 0.0)
		{
			//vel_[index_i][1] = coff * 5.0 * sqrt(3.0);
			vel_[index_i][2] = coff * (-5.0);
		}
	};
};
/**
 *  The main program
 */
int main(int ac, char* av[])
{
	/** Setup the system. Please the make sure the global domain bounds are correctly defined. */
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.setRunParticleRelaxation(false);
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.setReloadParticles(true);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	/** create a Plate body, corresponding material, particles and reaction model. */
	PDBody plate_body(sph_system, makeShared<Plate>("PDBody"));
	plate_body.defineParticlesAndMaterial<NosbPDParticles, HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
	plate_body.generateParticles<ParticleGeneratorLattice>();

	size_t particle_num_s = plate_body.getBaseParticles().total_real_particles_;

	SolidBody ball(sph_system, makeShared<GeometricShapeBall>(ball_center, ball_radius, "BallBody"));
	ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
	{
		ball.generateParticles<ParticleGeneratorReload>(io_environment, ball.getName());
	}
	else
	{
		ball.defineBodyLevelSetShape()->writeLevelSet(io_environment);
		ball.generateParticles<ParticleGeneratorLattice>();
	}
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		InnerRelation ball_inner(ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation for ball.
		//----------------------------------------------------------------------
		SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
		relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(io_environment, ball);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		ball_random_particles.parallel_exec(0.25);
		write_relaxed_particles.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begins.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			ball_relaxation_step_inner.parallel_exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relaxed_particles.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation plate_body_inner(plate_body);
	SurfaceContactRelation plate_ball_contact(plate_body, { &ball });

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	SimpleDynamics<PlateInitialCondition> plate_vel_initialization(plate_body);
	/** calculate shape Matrix */
	InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix>
		plate_shapeMatrix(plate_body_inner);
	/** Time step size calculation. */
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize>
		computing_time_step_size(plate_body);
	SimpleDynamics<TimeStepInitialization> initialize_a_solid_step(plate_body, makeShared<Gravity>(Vecd(0.0, 0.0, -gravity_g)));
	//stress relaxation for the beam by Hughes-Winget algorithm
	SimpleDynamics<solid_dynamics::NosbPDFirstStep> NosbPD_firstStep(plate_body);
	InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> NosbPD_secondStep(plate_body_inner);
	InteractionDynamics<solid_dynamics::NosbPDThirdStep> NosbPD_thirdStep(plate_body_inner);
	SimpleDynamics<solid_dynamics::NosbPDFourthStep> NosbPD_fourthStep(plate_body);
	//hourglass displacement mode control by LittleWood method
	InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass_control(plate_body_inner, plate_body.sph_adaptation_->getKernel());
	//breaking bonds based on the MAX principal stress criteria
	InteractionDynamics<solid_dynamics::BondBreakByPrinStress> check_bondLive(plate_body_inner,critical_stress);
	/** Algorithms for solid-solid contact. */
	InteractionDynamics<solid_dynamics::ContactDensitySummation> plate_update_contact_density(plate_ball_contact);
	InteractionDynamics<solid_dynamics::ContactForce> plate_compute_solid_contact_forces(plate_ball_contact);
	
	/** Output */
	BodyStatesRecordingToVtp write_states(io_environment, sph_system.real_bodies_);
	
	//Log file
	std::string Logpath = io_environment.output_folder_ + "/SimLog.txt";
	if (fs::exists(Logpath))
	{
		fs::remove(Logpath);
	}
	std::ofstream log_file(Logpath.c_str(), ios::trunc);
	std::cout << "# PARAM SETING #" << "\n" << "\n"
		<< "	particle_spacing_ref = " << resolution_ref << "\n"
		<< "	particle_num_s = " << particle_num_s << "\n" << "\n"
		<< "	rho0_s = " << rho0_s << "\n"
		<< "	Youngs_modulus = " << Youngs_modulus << "\n"
		<< "	poisson = " << poisson << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
		<< "# COMPUTATION START # " << "\n" << "\n";
	log_file << "#PARAM SETING#" << "\n" << "\n"
		<< "	particle_spacing_ref = " << resolution_ref << "\n"
		<< "	particle_num_s = " << particle_num_s << "\n" << "\n"
		<< "	rho0_s = " << rho0_s << "\n"
		<< "	Youngs_modulus = " << Youngs_modulus << "\n"
		<< "	poisson = " << poisson << "\n" << "\n"
		<< "	gravity_g = " << gravity_g << "\n" << "\n"
		<< "#COMPUTATION START HERE# " << "\n" << "\n";
	/**
	 * From here the time stepping begins.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	
	plate_shapeMatrix.parallel_exec();
	/** apply initial condition */
	plate_vel_initialization.parallel_exec();
	
	write_states.writeToFile(0);
	
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.3;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			ite++;
			dt = 0.1 * computing_time_step_size.parallel_exec();
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
			initialize_a_solid_step.parallel_exec(dt);

			plate_update_contact_density.parallel_exec(dt);
			plate_compute_solid_contact_forces.parallel_exec(dt);

			NosbPD_firstStep.parallel_exec(dt);

			check_bondLive.parallel_exec(dt);

			NosbPD_secondStep.parallel_exec(dt);

			hourglass_control.parallel_exec(dt);

			NosbPD_thirdStep.parallel_exec(dt);
			NosbPD_fourthStep.parallel_exec(dt);

			plate_body.updateCellLinkedList();
			plate_ball_contact.updateConfiguration();
			
			
		}
		
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
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
	
	return 0;
}
