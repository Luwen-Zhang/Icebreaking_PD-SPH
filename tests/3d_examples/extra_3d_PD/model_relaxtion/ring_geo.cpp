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

Real LR = 0.0036;				  // liquid column radius in XZ
Real LH = 12 * LR;				  // liquid column height in Y
Real inner_circle_radius = LR;
int resolution(20);

Real plate_Y = 0.0055;				/**< thickness of the ice plate. */
Real plate_R = 0.0575;	/**< width of the ice plate. */

Real resolution_ref = LR / 4;	  // particle spacing
Real BW = resolution_ref * 3; // boundary width
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
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vecd(-plate_R - BW, -BW, -plate_R - BW), Vecd(plate_R + BW, 2.0 * LH + BW, plate_R + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	RealBody ring_body(system, makeShared<PlateConstrainShape>("RingBody"));
	std::cout << "body created" << "\n";
	ring_body.defineParticlesAndMaterial();
	ring_body.generateParticles<ParticleGeneratorLattice>();

	std::cout << "particle_generated" << "\n";
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_imported_model_to_vtp(io_environment, { ring_body });
	//MeshRecordingToPlt cell_linked_list_recording(io_environment, ring_body.getCellLinkedList());
	
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	write_imported_model_to_vtp.writeToFile(0);
	

	return 0;
}
