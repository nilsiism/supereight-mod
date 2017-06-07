/*

 Copyright (c) 2014 University of Edinburgh, Imperial College, University of Manchester.
 Developed in the PAMELA project, EPSRC Programme Grant EP/K008730/1

 This code is licensed under the MIT License.

 */

#ifndef _KERNELS_
#define _KERNELS_

#include <cstdlib>
#include <commons.h>
#include <perfstats.h>
#include <timings.h>
#include <config.h>

/// OBJ ///

class Kfusion {
private:
	uint2 computationSize;
	float step;
	Matrix4 pose;
	Matrix4 *viewPose;
	float3 volumeDimensions;
	uint3 volumeResolution;
	std::vector<int> iterations;
	bool _tracked;
	bool _integrated;
	float3 _initPose;
  int voxel_block_size;
  float _mu;
  bool shouldRender = false;
  Configuration config;

	void raycast(uint frame, const float4& k, float mu);

public:
	Kfusion(uint2 inputSize, uint3 volumeResolution, float3 volumeDimensions,
			float3 initPose, std::vector<int> & pyramid, Configuration config):
			computationSize(make_uint2(inputSize.x, inputSize.y)) {

		this->_initPose = initPose;
		this->volumeDimensions = volumeDimensions;
		this->volumeResolution = volumeResolution;
    this->voxel_block_size = config.voxel_block_size;
    this->_mu = config.mu;
    this->config = config;

    pose.data[0] = {1.f, 0.f, 0.f, initPose.x};
    pose.data[1] = {0.f, 1.f, 0.f, initPose.y};
    pose.data[2] = {0.f, 0.f, 1.f, initPose.z};
    pose.data[3] = {0.f, 0.f, 0.f, 1.f};
		this->iterations.clear();
		for (std::vector<int>::iterator it = pyramid.begin();
				it != pyramid.end(); it++) {
			this->iterations.push_back(*it);
		}

		step = min(volumeDimensions) / max(volumeResolution);
		viewPose = &pose;
		this->languageSpecificConstructor();
	}
//Allow a kfusion object to be created with a pose which include orientation as well as position
	Kfusion(uint2 inputSize, uint3 volumeResolution, float3 volumeDimensions,
			Matrix4 initPose, std::vector<int> & pyramid, Configuration config) :
			computationSize(make_uint2(inputSize.x, inputSize.y)) {
		this->_initPose = getPosition();
		this->volumeDimensions = volumeDimensions;
		this->volumeResolution = volumeResolution;
        this->voxel_block_size = config.voxel_block_size;
        this->_mu = config.mu;
		pose = initPose;

		this->iterations.clear();
		for (std::vector<int>::iterator it = pyramid.begin();
				it != pyramid.end(); it++) {
			this->iterations.push_back(*it);
		}

		step = min(volumeDimensions) / max(volumeResolution);
		viewPose = &pose;
		this->languageSpecificConstructor();
	}

	void languageSpecificConstructor();
	~Kfusion();

	void reset();
	bool getTracked() {
		return (_tracked);
	}
	bool getIntegrated() {
		return (_integrated);
	}
	float3 getPosition() {
		//std::cerr << "InitPose =" << _initPose.x << "," << _initPose.y  <<"," << _initPose.z << "    ";
		//std::cerr << "pose =" << pose.data[0].w << "," << pose.data[1].w  <<"," << pose.data[2].w << "    ";
		float xt = pose.data[0].w - _initPose.x;
		float yt = pose.data[1].w - _initPose.y;
		float zt = pose.data[2].w - _initPose.z;
		return (make_float3(xt, yt, zt));
	}

    float3 getInitPos(){
        return _initPose;
    }

	inline void computeFrame(const ushort * inputDepth, const uint2 inputSize,
			float4 k, uint integration_rate, uint tracking_rate,
			float icp_threshold, float mu, const uint frame) {
		preprocessing(inputDepth, inputSize, true);
		_tracked = tracking(k, icp_threshold, tracking_rate, frame);
		_integrated = integration(k, integration_rate, mu, frame);
		raycasting(k, mu, frame);
	}

	bool preprocessing(const ushort * inputDepth, const uint2 inputSize, 
                     const bool filterInput);

	bool preprocessing(const ushort * inputDepth, const uchar3 * inputRGB,
                     const uint2 inputSize, const bool filterInput);

	bool tracking(float4 k, float icp_threshold, uint tracking_rate,
			uint frame);
	bool raycasting(float4 k, float mu, uint frame);
	bool integration(float4 k, uint integration_rate, float mu, uint frame);

  void dumpVolume(std::string filename);
  void printStats();

  void getPointCloudFromVolume();

  void renderVolume(uchar4 * out, const uint2 outputSize, int frame, int rate,
      float4 k, float mu);
  void renderTrack(uchar4 * out, const uint2 outputSize);
  void renderDepth(uchar4* out, uint2 outputSize);
	Matrix4 getPose() {
		return pose;
	}
	void setViewPose(Matrix4 *value = NULL) {
		if (value == NULL){
			viewPose = &pose;
      shouldRender = false;
    }
		else {
			viewPose = value;
      shouldRender = true;
    }
	}
	Matrix4 *getViewPose() {
		return (viewPose);
	}
	float3 getModelDimensions() {
		return (volumeDimensions);
	}
	uint3 getModelResolution() {
		return (volumeResolution);
	}
	uint2 getComputationResolution() {
		return (computationSize);
	}

};

void synchroniseDevices(); // Synchronise CPU and GPU

#endif
