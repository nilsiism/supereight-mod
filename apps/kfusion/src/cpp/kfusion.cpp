/*

 Copyright (c) 2014 University of Edinburgh, Imperial College, University of Manchester.
 Developed in the PAMELA project, EPSRC Programme Grant EP/K008730/1

 This code is licensed under the MIT License.

 */
#include <kernels.h>
#include <perfstats.h>
#include <vtk-io.h>
#include <octree.hpp>
#include "continuous/volume_instance.hpp"
#include "algorithms/meshing.hpp"
#include "geometry/octree_collision.hpp"

extern PerfStats Stats;

#ifdef __APPLE__
#include <mach/clock.h>
#include <mach/mach.h>

	
	#define TICK()    {if (print_kernel_timing) {\
		host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);\
		clock_get_time(cclock, &tick_clockData);\
		mach_port_deallocate(mach_task_self(), cclock);\
		}}

	#define TOCK(str,size)  {if (print_kernel_timing) {\
		host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);\
		clock_get_time(cclock, &tock_clockData);\
		mach_port_deallocate(mach_task_self(), cclock);\
		std::cerr<< str << " ";\
		if((tock_clockData.tv_sec > tick_clockData.tv_sec) && (tock_clockData.tv_nsec >= tick_clockData.tv_nsec))   std::cerr<< tock_clockData.tv_sec - tick_clockData.tv_sec << std::setfill('0') << std::setw(9);\
		std::cerr  << (( tock_clockData.tv_nsec - tick_clockData.tv_nsec) + ((tock_clockData.tv_nsec<tick_clockData.tv_nsec)?1000000000:0)) << " " <<  size << std::endl;}}
#else
	
	#define TICK()    {if (print_kernel_timing) {clock_gettime(CLOCK_MONOTONIC, &tick_clockData);}}

	#define TOCK(str,size)  {if (print_kernel_timing) {clock_gettime(CLOCK_MONOTONIC, &tock_clockData); /*std::cerr<< str << " "*/;\
		/*if((tock_clockData.tv_sec > tick_clockData.tv_sec) && (tock_clockData.tv_nsec >= tick_clockData.tv_nsec))*/ { Stats.sample(str, (( tock_clockData.tv_nsec - tick_clockData.tv_nsec) + ((tock_clockData.tv_nsec<tick_clockData.tv_nsec)?1000000000:0)), PerfStats::TIME);  /* std::cerr<< tock_clockData.tv_sec - tick_clockData.tv_sec << std::setfill('0') << std::setw(9);*/}\
		/*std::cerr  << (( tock_clockData.tv_nsec - tick_clockData.tv_nsec) + ((tock_clockData.tv_nsec<tick_clockData.tv_nsec)?1000000000:0)) << " " <<  size << std::endl;*/}}

#endif

bool print_kernel_timing = false;
#ifdef __APPLE__
	clock_serv_t cclock;
	mach_timespec_t tick_clockData;
	mach_timespec_t tock_clockData;
#else
	struct timespec tick_clockData;
	struct timespec tock_clockData;
#endif

#include "preprocessing.cpp"
#include "tracking.cpp"
#include "rendering.cpp"
#include "mapping.hpp"

// input once
float * gaussian;

// inter-frame
Volume<FieldType> volume;
float3 * vertex;
float3 * normal;

float3 bbox_min;
float3 bbox_max;

// intra-frame
TrackData * trackingResult;
float* reductionoutput;
float ** ScaledDepth;
float * floatDepth;
Matrix4 oldPose;
Matrix4 raycastPose;
float3 ** inputVertex;
float3 ** inputNormal;

bool bayesian = false;

  // For debugging purposes, will be deleted once done.
  std::vector<Matrix4> poses;


void Kfusion::languageSpecificConstructor() {

	if (getenv("KERNEL_TIMINGS"))
		print_kernel_timing = true;

	// internal buffers to initialize
	reductionoutput = (float*) calloc(sizeof(float) * 8 * 32, 1);

	ScaledDepth = (float**) calloc(sizeof(float*) * iterations.size(), 1);
	inputVertex = (float3**) calloc(sizeof(float3*) * iterations.size(), 1);
	inputNormal = (float3**) calloc(sizeof(float3*) * iterations.size(), 1);

	for (unsigned int i = 0; i < iterations.size(); ++i) {
		ScaledDepth[i] = (float*) calloc(
				sizeof(float) * (computationSize.x * computationSize.y)
						/ (int) pow(2, i), 1);
		inputVertex[i] = (float3*) calloc(
				sizeof(float3) * (computationSize.x * computationSize.y)
						/ (int) pow(2, i), 1);
		inputNormal[i] = (float3*) calloc(
				sizeof(float3) * (computationSize.x * computationSize.y)
						/ (int) pow(2, i), 1);
	}

	floatDepth = (float*) calloc(
			sizeof(float) * computationSize.x * computationSize.y, 1);
	vertex = (float3*) calloc(
			sizeof(float3) * computationSize.x * computationSize.y, 1);
	normal = (float3*) calloc(
			sizeof(float3) * computationSize.x * computationSize.y, 1);
	trackingResult = (TrackData*) calloc(
			sizeof(TrackData) * computationSize.x * computationSize.y, 1);

	// ********* BEGIN : Generate the gaussian *************
	size_t gaussianS = radius * 2 + 1;
	gaussian = (float*) calloc(gaussianS * sizeof(float), 1);
	int x;
	for (unsigned int i = 0; i < gaussianS; i++) {
		x = i - 2;
		gaussian[i] = expf(-(x * x) / (2 * delta * delta));
	}
	// ********* END : Generate the gaussian *************

  if(config.groundtruth_file != ""){
    parseGTFile(config.groundtruth_file, poses);
    std::cout << "Parsed " << poses.size() << " poses" << std::endl;
  }

  bayesian = config.bayesian;

	volume.init(volumeResolution.x, volumeDimensions.x);
	reset();
}

Kfusion::~Kfusion() {

	free(reductionoutput);
	for (unsigned int i = 0; i < iterations.size(); ++i) {
		free(ScaledDepth[i]);
		free(inputVertex[i]);
		free(inputNormal[i]);
	}
	free(ScaledDepth);
	free(inputVertex);
	free(inputNormal);

	free(vertex);
	free(normal);
	free(gaussian);

	volume.release();
}
void Kfusion::reset() {
}
void init() {
}
;
// stub
void clean() {
}
;

bool Kfusion::preprocessing(const ushort * inputDepth, const uint2 inputSize, 
    const bool filterInput){

    mm2metersKernel(floatDepth, computationSize, inputDepth, inputSize);
    if(filterInput){
	bilateralFilterKernel(ScaledDepth[0], floatDepth, computationSize, gaussian,
			e_delta, radius);
    }
    else {
      unsigned int y;
#pragma omp parallel for \
      shared(ScaledDepth), private(y)
      for (y = 0; y < computationSize.y; y++) {
        for (unsigned int x = 0; x < computationSize.x; x++) {
          ScaledDepth[0][x + y*computationSize.x] = floatDepth[x + y*computationSize.x]; 
        }
      }
    }
	return true;
}

bool Kfusion::tracking(float4 k, float icp_threshold, uint tracking_rate,
		uint frame) {

	if (frame % tracking_rate != 0)
		return false;

  if(!poses.empty()) {
    oldPose = this->pose;
    this->pose = poses[frame];
    this->pose.data[0].w = (this->pose.data[0].w - poses[0].data[0].w) + this->_initPose.x;
    this->pose.data[1].w = (this->pose.data[1].w - poses[0].data[1].w) + this->_initPose.y;
    this->pose.data[2].w = (this->pose.data[2].w - poses[0].data[2].w) + this->_initPose.z;
    return true;
  }

	// half sample the input depth maps into the pyramid levels
	for (unsigned int i = 1; i < iterations.size(); ++i) {
		halfSampleRobustImageKernel(ScaledDepth[i], ScaledDepth[i - 1],
				make_uint2(computationSize.x / (int) pow(2, i - 1),
						computationSize.y / (int) pow(2, i - 1)), e_delta * 3, 1);
	}

	// prepare the 3D information from the input depth maps
	uint2 localimagesize = computationSize;
	for (unsigned int i = 0; i < iterations.size(); ++i) {
		Matrix4 invK = getInverseCameraMatrix(k / float(1 << i));
		depth2vertexKernel(inputVertex[i], ScaledDepth[i], localimagesize,
				invK);
    if(k.y < 0)
      vertex2normalKernel<true>(inputNormal[i], inputVertex[i], localimagesize);
    else
      vertex2normalKernel<false>(inputNormal[i], inputVertex[i], localimagesize);
		localimagesize = make_uint2(localimagesize.x / 2, localimagesize.y / 2);
	}

	oldPose = pose;
	const Matrix4 projectReference = getCameraMatrix(k) * inverse(raycastPose);

	for (int level = iterations.size() - 1; level >= 0; --level) {
		uint2 localimagesize = make_uint2(
				computationSize.x / (int) pow(2, level),
				computationSize.y / (int) pow(2, level));
		for (int i = 0; i < iterations[level]; ++i) {

			trackKernel(trackingResult, inputVertex[level], inputNormal[level],
					localimagesize, vertex, normal, computationSize, pose,
					projectReference, dist_threshold, normal_threshold);

			reduceKernel(reductionoutput, trackingResult, computationSize,
					localimagesize);

			if (updatePoseKernel(pose, reductionoutput, icp_threshold))
				break;

		}
	}
	return checkPoseKernel(pose, oldPose, reductionoutput, computationSize,
			track_threshold);

}

bool Kfusion::raycasting(float4 k, float mu, uint frame) {

  bool doRaycast = false;

  if(frame > 2) {
    raycastPose = pose;
    raycastKernel(volume, vertex, normal, computationSize, 
        raycastPose * getInverseCameraMatrix(k), nearPlane, farPlane, mu, 
        step, step*BLOCK_SIDE);
    doRaycast = true;
  }
  return doRaycast;
}

bool Kfusion::integration(float4 k, uint integration_rate, float mu,
		uint frame) {

	bool doIntegrate = poses.empty() ? checkPoseKernel(pose, oldPose, reductionoutput,
			computationSize, track_threshold) : true;

  if ((doIntegrate && ((frame % integration_rate) == 0)) || (frame <= 3)) {

    volume.updateVolume(pose, getCameraMatrix(k), floatDepth, computationSize, mu, frame);
    // std::stringstream f;
    // f << "./slices/integration_" << frame << ".vtk";
    // save3DSlice(volume._map_index, make_int3(0, volume._size/2, 0),
    //   make_int3(volume._size, volume._size/2 + 1, volume._size), make_int3(volume._size), f.str().c_str());

    // f.str("");
    // f.clear();
    // f << "./slices/collision_" << frame << ".vtk";
    // save3DSlice(volume._map_index, [](const Octree<FieldType>& map,
    //       const int x, const int y, const int z) {
    //       const int3 bbox = make_int3(x, y, z);
    //       const int3 side = make_int3(1);
    //       auto test = [](const Octree<FieldType>::compute_type & val) {
    //       if(val.x == 0.f) return collision_status::unseen;
    //       if(val.x < 5) return collision_status::empty;
    //       return collision_status::occupied;
    //       };
    //       return (float) collides_with(map, bbox, side, test);
    //     }, // end lambda
    //     make_int3(0, volume._size/2, 0),
    //     make_int3(volume._size, volume._size/2 + 1, volume._size), make_int3(volume._size), f.str().c_str());
    doIntegrate = true;
  } else {
    doIntegrate = false;
  }

	return doIntegrate;

}

void Kfusion::dumpVolume(std::string ) {

}

void Kfusion::printStats(){
    int occupiedVoxels = 0;
    for(unsigned int x = 0; x < volume._size; x++){
        for(unsigned int y = 0; y < volume._size; y++){
            for(unsigned int z = 0; z < volume._size; z++){
                if( volume[make_uint3(x,y,z)].x < 1.f){
                    occupiedVoxels++;
                }  
            }
        }
     }
    std::cout << "The number of non-empty voxel is: " <<  occupiedVoxels << std::endl;
}

template <typename FieldType>
void raycastOrthogonal(Volume<FieldType> & volume, std::vector<float4> & points, const float3 origin, const float3 direction,
        const float farPlane, const float step) {

    // first walk with largesteps until we found a hit
    auto select_depth =  [](const auto& val) { return val.x; };
    float t = 0;
    float stepsize = step;
    float f_t = volume.interp(origin + direction * t, select_depth);
    t += step;
    float f_tt = 1.f;

    for (; t < farPlane; t += stepsize) {
      f_tt = volume.interp(origin + direction * t, select_depth);
      if ( (std::signbit(f_tt) != std::signbit(f_t))) {     // got it, jump out of inner loop
        if(f_t == 1.0 || f_tt == 1.0){
          f_t = f_tt;
          continue;
        }
        t = t + stepsize * f_tt / (f_t - f_tt);
        points.push_back(make_float4(origin + direction*t, 1));
      }
      if (f_tt < std::abs(0.8f))               // coming closer, reduce stepsize
        stepsize = step;
      f_t = f_tt;
    }
}


void Kfusion::getPointCloudFromVolume(){

  std::vector<float4> points;

  float x = 0, y = 0, z = 0;

  int3 resolution = make_int3(volume._size);
  float3 incr = make_float3(volume._dim / resolution.x, volume._dim / resolution.y, volume._dim / resolution.z);

  // XY plane

  std::cout << "Raycasting from XY plane.. " << std::endl;
  for(y = 0; y < volume._dim; y += incr.y ){
    for(x = 0; x < volume._dim; x += incr.x){
      raycastOrthogonal(volume, points, make_float3(x, y, 0), make_float3(0,0,1),
          volume._dim, step);
    }        
  }

  // ZY PLANE
  std::cout << "Raycasting from ZY plane.. " << std::endl;
  for(z = 0; z < volume._dim; z += incr.z ){
    for(y = 0; y < volume._dim; y += incr.y){
      raycastOrthogonal(volume, points, make_float3(0, y, z), make_float3(1,0,0),
          volume._dim, step);
    }
  }

  // ZX plane

  for(z = 0; z < volume._dim; z += incr.z ){
    for(x = 0;  x < volume._dim; x += incr.x){
      raycastOrthogonal(volume, points, make_float3(x, 0, z), make_float3(0,1,0),
          volume._dim, step);
    }
  }

  int num_points = points.size();
  std::cout << "Total number of ray-casted points : " << num_points << std::endl;

  if( !getenv("TRAJ")){
    std::cout << "Can't output the model point-cloud, unknown trajectory" << std::endl;
    return;
  }

  int trajectory = std::atoi(getenv("TRAJ"));
  std::stringstream filename;

  filename << "./pointcloud-vanilla-traj" << trajectory << "-" << volume._size << ".ply";

  // Matrix4 flipped = toMatrix4( TooN::SE3<float>(TooN::makeVector(0,0,0,0,0,0)));

  // flipped.data[0].w =  (-1 * this->_initPose.x); 
  // flipped.data[1].w =  (-1 * this->_initPose.y); 
  // flipped.data[2].w =  (-1 * this->_initPose.z); 

  //std::cout << "Generating point-cloud.. " << std::endl;
  //for(std::vector<float4>::iterator it = points.begin(); it != points.end(); ++it){
  //        float4 vertex = flipped * (*it); 
  //    }
}

void Kfusion::renderVolume(uchar4 * out, uint2 outputSize, int frame,
		int raycast_rendering_rate, float4 k, float largestep) {
	if (frame % raycast_rendering_rate == 0)
		renderVolumeKernel(volume, out, outputSize,
        *(this->viewPose) * getInverseCameraMatrix(k), nearPlane, 
        farPlane * 2.0f, _mu, step, largestep, 
        get_translation(*(this->viewPose)), ambient, 
        !compareMatrix4(*(this->viewPose), raycastPose));
}

void Kfusion::renderTrack(uchar4 * out, uint2 outputSize) {
	renderTrackKernel(out, trackingResult, outputSize);
}

void Kfusion::renderDepth(uchar4 * out, uint2 outputSize) {
	renderDepthKernel(out, floatDepth, outputSize, nearPlane, farPlane);
}

void Kfusion::dump_mesh(const char* filename){

  std::vector<Triangle> mesh;
  auto inside = [](const Volume<FieldType>::compute_type& val) {
    // meshing::status code;
    // if(val.y == 0.f) 
    //   code = meshing::status::UNKNOWN;
    // else 
    //   code = val.x < 0.f ? meshing::status::INSIDE : meshing::status::OUTSIDE;
    // return code;
    // std::cerr << val.x << " ";
    return val.x < 0.f;
  };

  auto select = [](const Volume<FieldType>::compute_type& val) {
    return val.x;
  };

  algorithms::marching_cube(volume._map_index, select, inside, mesh);
  writeVtkMesh(filename, mesh);
}

void synchroniseDevices() {
	// Nothing to do in the C++ implementation
}

