/* ---------------------------------------------------------------------------*
 *            SPHinXsys: 2D oscillation beam example-one body version           *
 * ----------------------------------------------------------------------------*
 * This is the one of the basic test cases, also the first case for            *
 * understanding SPH method for solid simulation.                              *
 * In this case, the constraint of the beam is implemented with                *
 * internal constrained subregion.                                             *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;	// beam length
Real PH = 0.02; // for thick plate; =0.01 for thin plate
Real SL = 0.06; // depth of the insert
// reference particle spacing
Real resolution_ref = PH / 10.0;
Real BW = resolution_ref * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
								 Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;		 // reference density
Real Youngs_modulus = 2.0e6; // reference Youngs modulus
Real poisson = 0.3975;		 // Poisson ratio
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
Real gravity_g = 0.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a beam base shape
std::vector<Vecd> beam_base_shape{
	Vecd(-SL - BW, -PH / 2 - BW), Vecd(-SL - BW, PH / 2 + BW), Vecd(0.0, PH / 2 + BW),
	Vecd(0.0, -PH / 2 - BW), Vecd(-SL - BW, -PH / 2 - BW)};
// a beam shape
std::vector<Vecd> beam_shape{
	Vecd(-SL, -PH / 2), Vecd(-SL, PH / 2), Vecd(PL, PH / 2), Vecd(PL, -PH / 2), Vecd(-SL, -PH / 2)};
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
public:
	explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class BeamInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit BeamInitialCondition(SPHBody &sph_body)
		: solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

	void update(size_t index_i, Real dt)
	{
		/** initial velocity profile */
		Real x = pos_[index_i][0] / PL;
		if (x > 0.0)
		{
			vel_[index_i][1] = vf * particles_->elastic_solid_.ReferenceSoundSpeed() *
							   (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
		}
	};
};
//----------------------------------------------------------------------
//	define the beam base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(beam_shape, ShapeBooleanOps::sub);
	return multi_polygon;
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	PDBody beam_body(system, makeShared<Beam>("PDBody"));
	beam_body.defineParticlesAndMaterial<NosbPDParticles, HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
	beam_body.generateParticles<ParticleGeneratorLattice>();

	

	size_t particle_num_s = beam_body.getBaseParticles().total_real_particles_;

	ObserverBody beam_observer(system, "BeamObserver");
	beam_observer.defineAdaptationRatios(1.5075, 2.0);
	beam_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation beam_body_inner(beam_body);
	ContactRelation beam_observer_contact(beam_observer, {&beam_body});
	//-----------------------------------------------------------------------------
	// this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);

	// calculate shape Matrix
	InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> beam_shapeMatrix(beam_body_inner);

	// time step size calculation
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(beam_body);
	
	SimpleDynamics<TimeStepInitialization> initialize_a_solid_step(beam_body, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
	//stress relaxation for the beam by Hughes-Winget algorithm
	SimpleDynamics<solid_dynamics::NosbPDFirstStep> NosbPD_firstStep(beam_body);
	InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> NosbPD_secondStep(beam_body_inner);
	InteractionDynamics<solid_dynamics::NosbPDThirdStep> NosbPD_thirdStep(beam_body_inner);

	// ADR_cn calculation
	ReduceDynamics<solid_dynamics::ADRFirstStep> computing_cn1(beam_body);
	ReduceDynamics<solid_dynamics::ADRSecondStep> computing_cn2(beam_body);
	
	SimpleDynamics<solid_dynamics::NosbPDFourthStepWithADR> NosbPD_fourthStepADR(beam_body);
	//hourglass displacement mode control by LittleWood method
	InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass_control(beam_body_inner, beam_body.sph_adaptation_->getKernel());

	// clamping a solid body part. This is softer than a direct constraint
	BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
	SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_beam_base(beam_base);
	//-----------------------------------------------------------------------------
	// outputs
	//-----------------------------------------------------------------------------
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_beam_states(io_environment, system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Vecd>>
		write_beam_tip_displacement("Position", io_environment, beam_observer_contact);
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
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	beam_shapeMatrix.exec();
	
	//----------------------------------------------------------------------
	//	Setup computing time-step controls.
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real end_time = T0;
	// time step size for output file
	Real output_interval = 0.01 * T0;
	Real Dt = 0.1 * output_interval; /**< Time period for data observing */
	Real dt = 0.0;					 // default acoustic time step sizes

	Real cn1 = 0.0;
	Real cn2 = 0.0;
	Real ADR_cn = 0.0;

	// statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//-----------------------------------------------------------------------------
	// from here the time stepping begins
	//-----------------------------------------------------------------------------
	write_beam_states.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);

	// computation loop starts
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		// integrate time (loop) until the next output time
		while (integration_time < output_interval)
		{

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				ite++;
				dt = 0.1 * computing_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				if (ite % 1000 == 0) {
					std::cout << "	N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
					log_file << "	N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				initialize_a_solid_step.parallel_exec(dt);

				NosbPD_firstStep.parallel_exec(dt);
				NosbPD_secondStep.parallel_exec(dt);

				hourglass_control.parallel_exec(dt);

				NosbPD_thirdStep.parallel_exec(dt);

				cn1 = SMAX(TinyReal, computing_cn1.parallel_exec(dt));
				cn2 = computing_cn2.parallel_exec(dt);

				if (cn2 > TinyReal) {
					ADR_cn = 2.0 * sqrt(cn1 / cn2);
				}
				else {
					ADR_cn = 0.0;
				}
				ADR_cn = 0.005 * ADR_cn;
				NosbPD_fourthStepADR.getADRcn(ADR_cn);
				NosbPD_fourthStepADR.parallel_exec(dt);

				constraint_beam_base.parallel_exec();

				
			}
		}

		write_beam_tip_displacement.writeToFile(ite);

		tick_count t2 = tick_count::now();
		write_beam_states.writeToFile();
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
