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
	beam_body.defineParticlesAndMaterial<NosbPDParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	beam_body.generateParticles<ParticleGeneratorLattice>();

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

	// corrected strong configuration
	InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> beam_shapeMatrix(beam_body_inner);

	// time step size calculation
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(beam_body);
	
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
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	beam_shapeMatrix.exec();
	

	return 0;
}