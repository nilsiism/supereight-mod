/*

 Copyright (c) 2014 University of Edinburgh, Imperial College, University of Manchester.
 Developed in the PAMELA project, EPSRC Programme Grant EP/K008730/1

 This code is licensed under the MIT License.

 */
#include <kernels.h>
#include <perfstats.h>
#include <vtk-io.h>
#include <volume.hpp>
#include <octree.hpp>
#include <algorithms/raycasting.cpp>

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
	
	#define TICK()    {if (true) {clock_gettime(CLOCK_MONOTONIC, &tick_clockData);}}

	#define TOCK(str,size)  {if (true) {clock_gettime(CLOCK_MONOTONIC, &tock_clockData); /*std::cerr<< str << " "*/;\
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

// input once
float * gaussian;

// inter-frame
Volume volume;
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

void integrateKernel(Volume& vol, const float* depth, uint2 depthSize,
		const Matrix4 invTrack, const Matrix4 K, const float mu,
		const float maxweight) {
	TICK();
	const float3 delta = rotate(invTrack,
			make_float3(0, 0, vol._dim / vol._size));
	const float3 cameraDelta = rotate(K, delta);
	unsigned int y;
#pragma omp parallel for \
        shared(vol), private(y)
	for (y = 0; y < vol._size; y++)
		for (unsigned int x = 0; x < vol._size; x++) {

			uint3 pix = make_uint3(x, y, 0); //pix.x = x;pix.y = y;
			float3 pos = invTrack * vol.pos(pix);
			float3 cameraX = K * pos;

			for (pix.z = 0; pix.z < vol._size;
					++pix.z, pos += delta, cameraX += cameraDelta) {
				if (pos.z < 0.0001f) // some near plane constraint
					continue;
				const float2 pixel = make_float2(cameraX.x / cameraX.z + 0.5f,
						cameraX.y / cameraX.z + 0.5f);
				if (pixel.x < 0 || pixel.x > depthSize.x - 1 || pixel.y < 0
						|| pixel.y > depthSize.y - 1)
					continue;
				const uint2 px = make_uint2(pixel.x, pixel.y);
				if (depth[px.x + px.y * depthSize.x] == 0)
					continue;
				const float diff =
						(depth[px.x + px.y * depthSize.x] - cameraX.z)
								* std::sqrt(
										1 + sq(pos.x / pos.z)
												+ sq(pos.y / pos.z));
				if (diff > -mu) {
					const float sdf = fminf(1.f, diff / mu);
					float2 data = vol[pix];
					data.x = clamp((data.y * data.x + sdf) / (data.y + 1), -1.f,
							1.f);
					data.y = fminf(data.y + 1, maxweight);
					vol.set(pix, data);
				}
			}
		}
	TOCK("integrateKernel", vol._size * vol._size);
}

float4 raycast(const Volume& volume, const uint2 pos, const Matrix4 view,
		const float nearPlane, const float farPlane, const float step,
		const float largestep) {

    float timings[4];
    int raycast_steps[4];

	const float3 origin = get_translation(view);
	const float3 direction = rotate(view, make_float3(pos.x, pos.y, 1.f));

	// intersect ray with a box
	// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
	// compute intersection of ray with all six bbox planes
	const float3 invR = make_float3(1.0f) / direction;
	const float3 tbot = -1 * invR * origin;
	const float3 ttop = invR * (volume._dim - origin);

	// re-order intersections to find smallest and largest on each axis
	const float3 tmin = fminf(ttop, tbot);
	const float3 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	const float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y),
			fmaxf(tmin.x, tmin.z));
	const float smallest_tmax = fminf(fminf(tmax.x, tmax.y),
			fminf(tmax.x, tmax.z));

	// check against near and far plane
	const float tnear = fmaxf(largest_tmin, nearPlane);
	const float tfar = fminf(smallest_tmax, farPlane);

	if (tnear < tfar) {
		// first walk with largesteps until we found a hit
		float t = tnear;
		float stepsize = largestep;
		float f_t = volume.interp(origin + direction * t);
		float f_tt = 0;
		if (f_t > 0) { // ups, if we were already in it, then don't render anything here
			for (; t < tfar; t += stepsize) {
				f_tt = volume.interp(origin + direction * t);
                if (f_tt < 0)                  // got it, jump out of inner loop
                    break;
				if (f_tt < 0.8f)               // coming closer, reduce stepsize
					stepsize = step;
				f_t = f_tt;
			}            
			if (f_tt < 0) {           // got it, calculate accurate intersection
				t = t + stepsize * f_tt / (f_t - f_tt);
				return make_float4(origin + direction * t, t);
			}
		}
	}
	return make_float4(0);

}

float4 raycastBayesian(const Volume& volume, const uint2 pos, const Matrix4 view,
		const float nearPlane, const float farPlane, const float step,
		const float largestep) {

    float timings[4];
    int raycast_steps[4];

	const float3 origin = get_translation(view);
	const float3 direction = rotate(view, make_float3(pos.x, pos.y, 1.f));

	// intersect ray with a box
	// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
	// compute intersection of ray with all six bbox planes
	const float3 invR = make_float3(1.0f) / direction;
	const float3 tbot = -1 * invR * origin;
	const float3 ttop = invR * (volume._dim - origin);

	// re-order intersections to find smallest and largest on each axis
	const float3 tmin = fminf(ttop, tbot);
	const float3 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	const float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y),
			fmaxf(tmin.x, tmin.z));
	const float smallest_tmax = fminf(fminf(tmax.x, tmax.y),
			fminf(tmax.x, tmax.z));

	// check against near and far plane
	const float tnear = fmaxf(largest_tmin, nearPlane);
	const float tfar = fminf(smallest_tmax, farPlane);

	if (tnear < tfar) {
		// first walk with largesteps until we found a hit
		float t = tnear;
		float stepsize = step;
		float f_t = volume.interp(origin + direction * t);
		float f_tt = 0;
		if (f_t >= 0) { // ups, if we were already in it, then don't render anything here
			for (; t < tfar; t += stepsize) {
				f_tt = volume.interp(origin + direction * t);
        if (f_tt > 0.5f )                  // got it, jump out of inner loop
          break;
				f_t = f_tt;
			}            
			if (f_tt > 0.5f) {           // got it, calculate accurate intersection
				t = t + stepsize * f_tt / (f_t - f_tt);
				return make_float4(origin + direction * t, t);
			}
		}
	}
	return make_float4(0);
}

void raycastKernel(float3* vertex, float3* normal, uint2 inputSize, 
    const Matrix4 view, const float nearPlane, const float farPlane, 
    const float mu, const float step, const float largestep, int frame) {
	TICK();
	unsigned int y;
#pragma omp parallel for \
	    shared(normal, vertex), private(y)
	for (y = 0; y < inputSize.y; y++)
		for (unsigned int x = 0; x < inputSize.x; x++) {

			uint2 pos = make_uint2(x, y);

			const float4 hit = algorithms::raycast(volume, pos, view, nearPlane, 
          farPlane, mu, step, largestep, frame);
			if (hit.w > 0.0) {
				vertex[pos.x + pos.y * inputSize.x] = make_float3(hit);
				float3 surfNorm = volume.grad(make_float3(hit));
				if (length(surfNorm) == 0) {
					//normal[pos] = normalize(surfNorm); // APN added
					normal[pos.x + pos.y * inputSize.x].x = INVALID;
				} else {
					normal[pos.x + pos.y * inputSize.x] = normalize(surfNorm);
				}
			} else {
				//std::cerr<< "RAYCAST MISS "<<  pos.x << " " << pos.y <<"  " << hit.w <<"\n";
				vertex[pos.x + pos.y * inputSize.x] = make_float3(0);
				normal[pos.x + pos.y * inputSize.x] = make_float3(INVALID, 0, 0);
			}
		}
	TOCK("raycastKernel", inputSize.x * inputSize.y);
}

bool updatePoseKernel(Matrix4 & pose, const float * output,
		float icp_threshold) {
	bool res = false;
	TICK();
	// Update the pose regarding the tracking result
	TooN::Matrix<8, 32, const float, TooN::Reference::RowMajor> values(output);
	TooN::Vector<6> x = solve(values[0].slice<1, 27>());
	TooN::SE3<> delta(x);
	pose = toMatrix4(delta) * pose;

	// Return validity test result of the tracking
	if (norm(x) < icp_threshold)
		res = true;

	TOCK("updatePoseKernel", 1);
	return res;
}

bool checkPoseKernel(Matrix4 & pose, Matrix4 oldPose, const float * output,
		uint2 imageSize, float track_threshold) {

	// Check the tracking result, and go back to the previous camera position if necessary

	TooN::Matrix<8, 32, const float, TooN::Reference::RowMajor> values(output);

	if ((std::sqrt(values(0, 0) / values(0, 28)) > 2e-2)
			|| (values(0, 28) / (imageSize.x * imageSize.y) < track_threshold)) {
		pose = oldPose;
		return false;
	} else {
		return true;
	}

}

void renderNormalKernel(uchar3* out, const float3* normal, uint2 normalSize) {
	TICK();
	unsigned int y;
#pragma omp parallel for \
        shared(out), private(y)
	for (y = 0; y < normalSize.y; y++)
		for (unsigned int x = 0; x < normalSize.x; x++) {
			uint pos = (x + y * normalSize.x);
			float3 n = normal[pos];
			if (n.x == -2) {
				out[pos] = make_uchar3(0, 0, 0);
			} else {
				n = normalize(n);
				out[pos] = make_uchar3(n.x * 128 + 128, n.y * 128 + 128,
						n.z * 128 + 128);
			}
		}
	TOCK("renderNormalKernel", normalSize.x * normalSize.y);
}

void renderDepthKernel(uchar4* out, float * depth, uint2 depthSize,
		const float nearPlane, const float farPlane) {
	TICK();

	float rangeScale = 1 / (farPlane - nearPlane);

	unsigned int y;
#pragma omp parallel for \
        shared(out), private(y)
	for (y = 0; y < depthSize.y; y++) {
		int rowOffeset = y * depthSize.x;
		for (unsigned int x = 0; x < depthSize.x; x++) {

			unsigned int pos = rowOffeset + x;

			if (depth[pos] < nearPlane)
				out[pos] = make_uchar4(255, 255, 255, 0); // The forth value is a padding in order to align memory
			else {
				if (depth[pos] > farPlane)
					out[pos] = make_uchar4(0, 0, 0, 0); // The forth value is a padding in order to align memory
				else {
					const float d = (depth[pos] - nearPlane) * rangeScale;
					out[pos] = gs2rgb(d);
				}
			}
		}
	}
	TOCK("renderDepthKernel", depthSize.x * depthSize.y);
}

void renderTrackKernel(uchar4* out, const TrackData* data, uint2 outSize) {
	TICK();

	unsigned int y;
#pragma omp parallel for \
        shared(out), private(y)
	for (y = 0; y < outSize.y; y++)
		for (unsigned int x = 0; x < outSize.x; x++) {
			uint pos = x + y * outSize.x;
			switch (data[pos].result) {
			case 1:
				out[pos] = make_uchar4(128, 128, 128, 0);  // ok	 GREY
				break;
			case -1:
				out[pos] = make_uchar4(0, 0, 0, 0);      // no input BLACK
				break;
			case -2:
				out[pos] = make_uchar4(255, 0, 0, 0);        // not in image RED
				break;
			case -3:
				out[pos] = make_uchar4(0, 255, 0, 0);    // no correspondence GREEN
				break;
			case -4:
				out[pos] = make_uchar4(0, 0, 255, 0);        // to far away BLUE
				break;
			case -5:
				out[pos] = make_uchar4(255, 255, 0, 0);     // wrong normal YELLOW
				break;
			default:
				out[pos] = make_uchar4(255, 128, 128, 0);
				break;
			}
		}
	TOCK("renderTrackKernel", outSize.x * outSize.y);
}

void renderVolumeKernel(uchar4* out, const uint2 depthSize, const Matrix4 view, 
    const float nearPlane, const float farPlane, const float mu,
		const float step, const float largestep, const float3 light,
		const float3 ambient, int frame) {
	TICK();
	unsigned int y;
#pragma omp parallel for \
        shared(out), private(y)
	for (y = 0; y < depthSize.y; y++) {
		for (unsigned int x = 0; x < depthSize.x; x++) {
			const uint2 pos = make_uint2(x, y);

			const float4 hit = algorithms::raycast(volume, pos, view, nearPlane, 
          farPlane, mu, step, largestep, frame);
			if (hit.w > 0) {
				const float3 test = make_float3(hit);
				const float3 surfNorm = volume.grad(test);
				if (length(surfNorm) > 0) {
					const float3 diff = normalize(light - test);
					const float dir = fmaxf(dot(normalize(surfNorm), diff),
							0.f);
					const float3 col = clamp(make_float3(dir) + ambient, 0.f,
							1.f) * 255;
					out[x + depthSize.x*y] = make_uchar4(col.x, col.y, col.z, 0); // The forth value is a padding to align memory
				} else {
					out[x + depthSize.x*y] = make_uchar4(0, 0, 0, 0); // The forth value is a padding to align memory
				}
			} else {
				out[x + depthSize.x*y] = make_uchar4(0, 0, 0, 0); // The forth value is a padding to align memory
			}
		}
	}
	TOCK("renderVolumeKernel", depthSize.x * depthSize.y);
}

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

bool Kfusion::preprocessing(const ushort * inputDepth, const uchar3 * inputRGB, 
    const uint2 inputSize, const bool filterInput){

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
      vertex2normalKernel<true>(inputNormal[i], inputVertex[i], localimagesize, k);
    else
      vertex2normalKernel<false>(inputNormal[i], inputVertex[i], localimagesize, k);
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
    raycastKernel(vertex, normal, computationSize, 
        raycastPose * getInverseCameraMatrix(k), nearPlane, farPlane, mu, 
        step, step*BLOCK_SIDE, frame);
    doRaycast = true;
  }
  return doRaycast;
}

bool Kfusion::integration(float4 k, uint integration_rate, float mu,
		uint frame) {

	bool doIntegrate = poses.empty() ? checkPoseKernel(pose, oldPose, reductionoutput,
			computationSize, track_threshold) : true;

  if ((doIntegrate && ((frame % integration_rate) == 0)) || (frame <= 3)) {
    // integrateKernel(volume, floatDepth, computationSize, inverse(pose),
    //     getCameraMatrix(k), mu, maxweight);
    volume.updateVolume(pose, k, floatDepth, computationSize, mu, frame);
    // std::stringstream f;
    // f << "./slices/integration_" << (bayesian ? "bayesian_" : "tsdf_" )<< frame << ".vtk";
    // save3DSlice(volume._data.raw_data(), make_uint3(0, volume._size/2, 0),
    //   make_uint3(volume._size, volume._size/2 + 1, volume._size), make_uint3(volume._size), f.str().c_str());
    doIntegrate = true;
  } else {
    doIntegrate = false;
  }

	return doIntegrate;

}

void Kfusion::dumpVolume(std::string filename) {

	std::ofstream fDumpFile;

	if (filename == "") {
		return;
	}

	std::cout << "Dumping the volumetric representation on file: " << filename
			<< std::endl;
	fDumpFile.open((filename + ".config").c_str(), std::ios::out);
	if (fDumpFile.fail()) {
		std::cout << "Error opening file: " << filename << std::endl;
		exit(1);
	}

    fDumpFile << "/*** CONFIG ***/" << std::endl;
    fDumpFile << "resolution:" << volume._size << std::endl;
    fDumpFile << "dimension:" << volume._dim << std::endl;
    fDumpFile << "offset:" << _initPose.x << "," << _initPose.y << "," << _initPose.z << std::endl;

    fDumpFile.close();

	fDumpFile.open((filename + ".data").c_str(), std::ios::binary);
	if (fDumpFile.fail()) {
		std::cout << "Error opening file: " << filename << std::endl;
		exit(1);
	}
    
  
    for(int z = 0; z < volume._size; z++){
        for(int y = 0; y < volume._size; y++){
            for(int x = 0; x < volume._size; x++){
                uint3 pos = make_uint3(x,y,z);
                float2 data = volume[pos];
		        fDumpFile.write(reinterpret_cast<char *>(&data.x), sizeof(float));
                fDumpFile.write(reinterpret_cast<char *>(&data.y), sizeof(float));
            }
        }
    }

    fDumpFile.close();
}

void Kfusion::printStats(){
    int occupiedVoxels = 0;
    for(int x = 0; x < volume._size; x++){
        for(int y = 0; y < volume._size; y++){
            for(int z = 0; z < volume._size; z++){
                if( volume[make_uint3(x,y,z)].x < 1.f){
                    occupiedVoxels++;
                }  
            }
        }
     }
    std::cout << "The number of non-empty voxel is: " <<  occupiedVoxels << std::endl;
}

void raycastOrthogonal( Volume & volume, std::vector<float4> & points, const float3 origin, const float3 direction,
        const float farPlane, const float step, const float largestep) {

    float timings[4];
    int raycast_steps[4];

    // first walk with largesteps until we found a hit
    float t = 0;
    float stepsize = step;
    float f_t = volume.interp(origin + direction * t);
    t += step;
    float f_tt = 1.f;

    for (; t < farPlane; t += stepsize) {
      f_tt = volume.interp(origin + direction * t);
      if ( (std::signbit(f_tt) != std::signbit(f_t))) {     // got it, jump out of inner loop
        float2 data_t = volume[origin + direction * (t-stepsize)];
        float2 data_tt = volume[origin + direction * t];
        if(f_t == 1.0 || f_tt == 1.0 || data_t.y < 6|| data_tt.y < 6 ){
          f_t = f_tt;
          continue;
        }

        // Else, we have found a good intersection, calculate it.
        //            std::cout << "Interpolating on ray: " << origin.x << ", " << origin.y << ", " << origin.z << 
        //               " at distance " << t <<  "  between f_t: " << f_t << " and f_tt: " << f_tt << " Weights: " << data_t.y << ", " << data_tt.y << std::endl;
        t = t + stepsize * f_tt / (f_t - f_tt);
        points.push_back(make_float4(origin + direction*t, 1));
      }
      if (f_tt < std::abs(0.8f))               // coming closer, reduce stepsize
        stepsize = step;
      f_t = f_tt;
    }
}


void Kfusion::getPointCloudFromVolume(float mu){

    std::vector<float4> points;

    float x = 0, y = 0, z = 0;

    int3 resolution = make_int3(volume._size);
    float3 incr = make_float3(volume._dim / resolution.x, volume._dim / resolution.y, volume._dim / resolution.z);

    //#pragma omp parallel for \
//        shared(points), private(y, posy)

    // XY plane

    std::cout << "Raycasting from XY plane.. " << std::endl;
    for(y = 0; y < volume._dim; y += incr.y ){
        for(x = 0; x < volume._dim; x += incr.x){
            raycastOrthogonal(volume, points, make_float3(x, y, 0), make_float3(0,0,1),
                    volume._dim, step,  step);
        }        
    }

    // ZY PLANE
    std::cout << "Raycasting from ZY plane.. " << std::endl;
    for(z = 0; z < volume._dim; z += incr.z ){
        for(y = 0; y < volume._dim; y += incr.y){
            raycastOrthogonal(volume, points, make_float3(0, y, z), make_float3(1,0,0),
                    volume._dim, step,  step);
        }
    }


    // ZX plane

    for(z = 0; z < volume._dim; z += incr.z ){
        for(x = 0;  x < volume._dim; x += incr.x){
            raycastOrthogonal(volume, points, make_float3(x, 0, z), make_float3(0,1,0),
                    volume._dim, step,  step);
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

    Matrix4 flipped = toMatrix4( TooN::SE3<float>(TooN::makeVector(0,0,0,0,0,0)));

    flipped.data[0].w =  (-1 * this->_initPose.x); 
    flipped.data[1].w =  (-1 * this->_initPose.y); 
    flipped.data[2].w =  (-1 * this->_initPose.z); 

    //std::cout << "Generating point-cloud.. " << std::endl;
    //for(std::vector<float4>::iterator it = points.begin(); it != points.end(); ++it){
    //        float4 vertex = flipped * (*it); 
    //    }
}

void Kfusion::renderVolume(uchar4 * out, uint2 outputSize, int frame,
		int raycast_rendering_rate, float4 k, float largestep) {
	if (frame % raycast_rendering_rate == 0)
		renderVolumeKernel(out, outputSize, 
        *(this->viewPose) * getInverseCameraMatrix(k), nearPlane, 
        farPlane * 2.0f, _mu, step, largestep, light, ambient, frame);
}

void Kfusion::renderTrack(uchar4 * out, uint2 outputSize) {
	renderTrackKernel(out, trackingResult, outputSize);
}

void Kfusion::renderDepth(uchar4 * out, uint2 outputSize) {
	renderDepthKernel(out, floatDepth, outputSize, nearPlane, farPlane);
}

void synchroniseDevices() {
	// Nothing to do in the C++ implementation
}
