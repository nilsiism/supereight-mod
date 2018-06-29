/*

 Copyright (c) 2014 University of Edinburgh, Imperial College, University of Manchester.
 Developed in the PAMELA project, EPSRC Programme Grant EP/K008730/1

 This code is licensed under the MIT License.


 Copyright 2016 Emanuele Vespa, Imperial College London 

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 */

#include <se/DenseSLAMSystem.h>
#include <se/ray_iterator.hpp>
#include <se/algorithms/meshing.hpp>
#include <se/geometry/octree_collision.hpp>
#include <se/vtk-io.h>
#include "timings.h"
#include <perfstats.h>
#include "preprocessing.cpp"
#include "tracking.cpp"
#include "rendering.cpp"
#include "bfusion/mapping_impl.hpp"
#include "kfusion/mapping_impl.hpp"
#include "bfusion/alloc_impl.hpp"
#include "kfusion/alloc_impl.hpp"


extern PerfStats Stats;

// For debugging purposes, will be deleted once done.
std::vector<Matrix4> poses;

DenseSLAMSystem::DenseSLAMSystem(uint2 inputSize, uint3 volumeResolution, float3 volumeDimensions,
			float3 initPose, std::vector<int> & pyramid, Configuration config):
  computation_size_(make_uint2(inputSize.x, inputSize.y)) {

    this->init_pose_ = initPose;
    this->volume_dimension_ = volumeDimensions;
    this->volume_resolution_ = volumeResolution;
    this->voxel_block_size = config.voxel_block_size;
    this->_mu = config.mu;
    this->config = config;

    pose_.data[0] = {1.f, 0.f, 0.f, initPose.x};
    pose_.data[1] = {0.f, 1.f, 0.f, initPose.y};
    pose_.data[2] = {0.f, 0.f, 1.f, initPose.z};
    pose_.data[3] = {0.f, 0.f, 0.f, 1.f};
    this->iterations_.clear();
    for (std::vector<int>::iterator it = pyramid.begin();
        it != pyramid.end(); it++) {
      this->iterations_.push_back(*it);
    }

    step = min(volumeDimensions) / max(volumeResolution);
    viewPose_ = &pose_;
    this->languageSpecificConstructor();
  }

DenseSLAMSystem::DenseSLAMSystem(uint2 inputSize, uint3 volumeResolution, 
    float3 volumeDimensions, Matrix4 initPose, std::vector<int> & pyramid, 
    Configuration config) :
  computation_size_(make_uint2(inputSize.x, inputSize.y)) {
    this->init_pose_ = getPosition();
    this->volume_dimension_ = volumeDimensions;
    this->volume_resolution_ = volumeResolution;
    this->voxel_block_size = config.voxel_block_size;
    this->_mu = config.mu;
    pose_ = initPose;

    this->iterations_.clear();
    for (std::vector<int>::iterator it = pyramid.begin();
        it != pyramid.end(); it++) {
      this->iterations_.push_back(*it);
    }

    step = min(volumeDimensions) / max(volumeResolution);
    viewPose_ = &pose_;
    this->languageSpecificConstructor();
  }

void DenseSLAMSystem::languageSpecificConstructor() {

	if (getenv("KERNEL_TIMINGS"))
		print_kernel_timing = true;

	// internal buffers to initialize
	reductionoutput = (float*) calloc(sizeof(float) * 8 * 32, 1);

	ScaledDepth = (float**) calloc(sizeof(float*) * iterations_.size(), 1);
	inputVertex = (float3**) calloc(sizeof(float3*) * iterations_.size(), 1);
	inputNormal = (float3**) calloc(sizeof(float3*) * iterations_.size(), 1);

	for (unsigned int i = 0; i < iterations_.size(); ++i) {
		ScaledDepth[i] = (float*) calloc(
				sizeof(float) * (computation_size_.x * computation_size_.y)
						/ (int) pow(2, i), 1);
		inputVertex[i] = (float3*) calloc(
				sizeof(float3) * (computation_size_.x * computation_size_.y)
						/ (int) pow(2, i), 1);
		inputNormal[i] = (float3*) calloc(
				sizeof(float3) * (computation_size_.x * computation_size_.y)
						/ (int) pow(2, i), 1);
	}

	floatDepth = (float*) calloc(
			sizeof(float) * computation_size_.x * computation_size_.y, 1);
	vertex = (float3*) calloc(
			sizeof(float3) * computation_size_.x * computation_size_.y, 1);
	normal = (float3*) calloc(
			sizeof(float3) * computation_size_.x * computation_size_.y, 1);
	trackingResult = (TrackData*) calloc(
			sizeof(TrackData) * computation_size_.x * computation_size_.y, 1);

  allocationList = NULL;
  reserved = 0;
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

  volume_ptr = std::make_shared<Volume<FieldType> >();
        volume_ptr->init(volume_resolution_.x, volume_dimension_.x);
}

DenseSLAMSystem::~DenseSLAMSystem() {

	free(reductionoutput);
	for (unsigned int i = 0; i < iterations_.size(); ++i) {
		free(ScaledDepth[i]);
		free(inputVertex[i]);
		free(inputNormal[i]);
	}
	free(ScaledDepth);
	free(inputVertex);
	free(vertex);
	free(normal);
	free(gaussian);
  
  if(allocationList) free(allocationList);
	volume_ptr->release();
}

bool DenseSLAMSystem::preprocessing(const ushort * inputDepth, const uint2 inputSize, 
    const bool filterInput){

    mm2metersKernel(floatDepth, computation_size_, inputDepth, inputSize);
    if(filterInput){
        bilateralFilterKernel(ScaledDepth[0], floatDepth, computation_size_, gaussian,
			e_delta, radius);
    }
    else {
      std::memcpy(ScaledDepth[0], floatDepth, 
          sizeof(float) * computation_size_.x * computation_size_.y);
    }
	return true;
}

bool DenseSLAMSystem::tracking(float4 k, float icp_threshold, uint tracking_rate,
		uint frame) {

	if (frame % tracking_rate != 0)
		return false;

  if(!poses.empty()) {
    oldPose = this->pose_;
    this->pose_ = poses[frame];
    this->pose_.data[0].w = (this->pose_.data[0].w - poses[0].data[0].w) + this->init_pose_.x;
    this->pose_.data[1].w = (this->pose_.data[1].w - poses[0].data[1].w) + this->init_pose_.y;
    this->pose_.data[2].w = (this->pose_.data[2].w - poses[0].data[2].w) + this->init_pose_.z;
    return true;
  }

	// half sample the input depth maps into the pyramid levels
	for (unsigned int i = 1; i < iterations_.size(); ++i) {
		halfSampleRobustImageKernel(ScaledDepth[i], ScaledDepth[i - 1],
				make_uint2(computation_size_.x / (int) pow(2, i - 1),
						computation_size_.y / (int) pow(2, i - 1)), e_delta * 3, 1);
	}

	// prepare the 3D information from the input depth maps
	uint2 localimagesize = computation_size_;
	for (unsigned int i = 0; i < iterations_.size(); ++i) {
		Matrix4 invK = getInverseCameraMatrix(k / float(1 << i));
		depth2vertexKernel(inputVertex[i], ScaledDepth[i], localimagesize,
				invK);
    if(k.y < 0)
      vertex2normalKernel<FieldType, true>(inputNormal[i], inputVertex[i], localimagesize);
    else
      vertex2normalKernel<FieldType, false>(inputNormal[i], inputVertex[i], localimagesize);
		localimagesize = make_uint2(localimagesize.x / 2, localimagesize.y / 2);
	}

	oldPose = pose_;
	const Matrix4 projectReference = getCameraMatrix(k) * inverse(raycastPose);

	for (int level = iterations_.size() - 1; level >= 0; --level) {
		uint2 localimagesize = make_uint2(
				computation_size_.x / (int) pow(2, level),
				computation_size_.y / (int) pow(2, level));
		for (int i = 0; i < iterations_[level]; ++i) {

			trackKernel(trackingResult, inputVertex[level], inputNormal[level],
					localimagesize, vertex, normal, computation_size_, pose_,
					projectReference, dist_threshold, normal_threshold);

			reduceKernel(reductionoutput, trackingResult, computation_size_,
					localimagesize);

			if (updatePoseKernel(pose_, reductionoutput, icp_threshold))
				break;

		}
	}
	return checkPoseKernel(pose_, oldPose, reductionoutput, computation_size_,
			track_threshold);

}

bool DenseSLAMSystem::raycasting(float4 k, float mu, uint frame) {

  bool doRaycast = false;

  if(frame > 2) {
    raycastPose = pose_;
    raycastKernel(*volume_ptr, vertex, normal, computation_size_,
        raycastPose * getInverseCameraMatrix(k), nearPlane, farPlane, mu, 
        step, step*BLOCK_SIDE);
    doRaycast = true;
  }
  return doRaycast;
}

bool DenseSLAMSystem::integration(float4 k, uint integration_rate, float mu,
		uint frame) {

	bool doIntegrate = poses.empty() ? checkPoseKernel(pose_, oldPose, reductionoutput,
			computation_size_, track_threshold) : true;

  if ((doIntegrate && ((frame % integration_rate) == 0)) || (frame <= 3)) {

    float voxelsize =  volume_ptr->_dim/volume_ptr->_size;
    int num_vox_per_pix = volume_ptr->_dim/((se::VoxelBlock<FieldType>::side)*voxelsize);
    size_t total = num_vox_per_pix * computation_size_.x * computation_size_.y;
    if(!reserved) {
      allocationList = (se::key_t* ) calloc(sizeof(se::key_t) * total, 1);
      reserved = total;
    }
    unsigned int allocated = 0;
    if(std::is_same<FieldType, SDF>::value) {
     allocated  = buildAllocationList(allocationList, reserved, 
        volume_ptr->_map_index, pose_, getCameraMatrix(k), floatDepth, computation_size_, volume_ptr->_size,
      voxelsize, 2*mu);  
    } else if(std::is_same<FieldType, OFusion>::value) {
     allocated = buildOctantList(allocationList, reserved, volume_ptr->_map_index,
         pose_, getCameraMatrix(k), floatDepth, computation_size_, voxelsize,
         compute_stepsize, step_to_depth, 6*mu);  
    }

    volume_ptr->_map_index.allocate(allocationList, allocated);

    if(std::is_same<FieldType, SDF>::value) {
      struct sdf_update funct(floatDepth, 
          Eigen::Vector2i(computation_size_.x, computation_size_.y), mu, 100);
      se::functor::projective_map(volume_ptr->_map_index, 
          to_sophus(pose_).inverse(),
          to_eigen(getCameraMatrix(k)), 
          Eigen::Vector2i(computation_size_.x, computation_size_.y),
          funct);
    } else if(std::is_same<FieldType, OFusion>::value) {

      float timestamp = (1.f/30.f)*frame; 
      struct bfusion_update funct(floatDepth, 
          Eigen::Vector2i(computation_size_.x, computation_size_.y), mu, timestamp);

      se::functor::projective_map(volume_ptr->_map_index, 
          to_sophus(pose_).inverse(),
          to_eigen(getCameraMatrix(k)), 
          Eigen::Vector2i(computation_size_.x, computation_size_.y),
          funct);
    }

    // if(frame % 15 == 0) {
    //   std::stringstream f;
    //   f << "./slices/integration_" << frame << ".vtk";
    //   save3DSlice(volume_ptr->_map_index, make_int3(0, volume_ptr->_size/2, 0),
    //       make_int3(volume_ptr->_size, volume_ptr->_size/2 + 1, volume_ptr->_size), make_int3(volume_ptr->_size), f.str().c_str());
    //   f.str("");
    //   f.clear();
    //   }

    // f << "./slices/collision_" << frame << ".vtk";
    // save3DSlice(volume_ptr->_map_index, [](const Octree<FieldType>& map,
    //       const int x, const int y, const int z) {
    //       const int3 bbox = make_int3(x, y, z);
    //       const int3 side = make_int3(1);
    //       auto test = [](const Octree<FieldType>::value_type & val) {
    //       if(val.x == 0.f) return collision_status::unseen;
    //       if(val.x < 5) return collision_status::empty;
    //       return collision_status::occupied;
    //       };
    //       return (float) collides_with(map, bbox, side, test);
    //     }, // end lambda
    //     make_int3(0, volume_ptr->_size/2, 0),
    //     make_int3(volume_ptr->_size, volume_ptr->_size/2 + 1, volume_ptr->_size), make_int3(volume_ptr->_size), f.str().c_str());
    doIntegrate = true;
  } else {
    doIntegrate = false;
  }

	return doIntegrate;

}

void DenseSLAMSystem::dump_volume(std::string ) {

}

void DenseSLAMSystem::renderVolume(uchar4 * out, uint2 outputSize, int frame,
		int raycast_rendering_rate, float4 k, float largestep) {
	if (frame % raycast_rendering_rate == 0)
		renderVolumeKernel(*volume_ptr, out, outputSize,
	*(this->viewPose_) * getInverseCameraMatrix(k), nearPlane,
        farPlane * 2.0f, _mu, step, largestep, 
        get_translation(*(this->viewPose_)), ambient,
        !compareMatrix4(*(this->viewPose_), raycastPose), vertex, normal);
}

void DenseSLAMSystem::renderTrack(uchar4 * out, uint2 outputSize) {
	renderTrackKernel(out, trackingResult, outputSize);
}

void DenseSLAMSystem::renderDepth(uchar4 * out, uint2 outputSize) {
	renderDepthKernel(out, floatDepth, outputSize, nearPlane, farPlane);
}

void DenseSLAMSystem::dump_mesh(const std::string filename){

  std::vector<Triangle> mesh;
  auto inside = [](const Volume<FieldType>::value_type& val) {
    // meshing::status code;
    // if(val.y == 0.f) 
    //   code = meshing::status::UNKNOWN;
    // else 
    //   code = val.x < 0.f ? meshing::status::INSIDE : meshing::status::OUTSIDE;
    // return code;
    // std::cerr << val.x << " ";
    return val.x < 0.f;
  };

  auto select = [](const Volume<FieldType>::value_type& val) {
    return val.x;
  };

  se::algorithms::marching_cube(volume_ptr->_map_index, select, inside, mesh);
  writeVtkMesh(filename.c_str(), mesh);
}

void synchroniseDevices() {
	// Nothing to do in the C++ implementation
}

