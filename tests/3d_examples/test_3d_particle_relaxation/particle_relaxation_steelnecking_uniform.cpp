/**
 * @file 	particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate body fitted particles (3D).
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

// general parameters for geometry

Real bar_Y = 2.0 * 0.02667;	/**< length of the metal bar. */
Real bar_R = 0.006413;	/**< radius of the metal bar. */
int resolution(20);
Real resolution_ref = bar_R / 12;	  // particle spacing
Real BW = resolution_ref * 10; // boundary width
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/columnNecking3.stl";
//std::string full_path_to_file = "D:/SPHinXsys-master-build/tests/3d_examples/test_3d_particle_relaxation/bin/input/steelnecking.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-bar_R - BW, -BW, -bar_R - BW);
Vec3d domain_upper_bound(bar_R + BW, bar_Y + BW, bar_R + BW);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 12.5;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	define a body from the imported model.
//----------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
public:
	explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd translation(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	RealBody imported_model(system, makeShared<SolidBodyFromMesh>("SolidBodyFromMesh"));
	std::cout << "body created" << "\n";
	//imported_model.defineAdaptation<ParticleRefinementNearSurface>(1.15, 1.0, 1);
	//imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(io_environment);
	imported_model.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	imported_model.defineParticlesAndMaterial();
	//imported_model.generateParticles<ParticleGeneratorMultiResolution>();
	imported_model.generateParticles<ParticleGeneratorLattice>();
	//imported_model.addBodyStateForRecording<Real>("SmoothingLengthRatio");

	std::cout << "particle_generated" << "\n";
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_imported_model_to_vtp(io_environment, {imported_model});
	MeshRecordingToPlt cell_linked_list_recording(io_environment, imported_model.getCellLinkedList());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	//AdaptiveInnerRelation imported_model_inner(imported_model);
	InnerRelation imported_model_inner(imported_model);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
	//SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(imported_model);
	
	std::cout << "method defined" << "\n";
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	write_imported_model_to_vtp.writeToFile(0);
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.SurfaceBounding().parallel_exec();
	//update_smoothing_length_ratio.parallel_exec();
	write_imported_model_to_vtp.writeToFile(1);
	//imported_model.updateCellLinkedList();
	//cell_linked_list_recording.writeToFile(0);
	 std::cout << "main loop starts" << "\n";
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	int ite_pf = 0;
	while (ite_p < 1000)
	{
		//update_smoothing_length_ratio.parallel_exec();
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			ite_pf = ite_p + 2;
			write_imported_model_to_vtp.writeToFile(ite_pf);
		}
	}
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;

	return 0;
}
