/*

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

#ifndef CONFIG_H
#define CONFIG_H

#include <vector_types.h>
#include <vector>
#include <string>

struct Configuration {

	// 
  // KFusion configuration parameters
  // Command line arguments are parsed in default_parameters.h 
  //

	int compute_size_ratio;
	int integration_rate;
	int rendering_rate;
	int tracking_rate;
	uint3 volume_resolution;
	float3 volume_size;
    int voxel_block_size;
	float3 initial_pos_factor;
	std::vector<int> pyramid;
	std::string dump_volume_file;
	std::string input_file;
	std::string log_file;
  std::string groundtruth_file;

	float4 camera;
	bool camera_overrided;

	float mu;
	int fps;
	bool blocking_read;
	float icp_threshold;
	bool no_gui;
	bool render_volume_fullsize;
  bool bilateralFilter;
  bool colouredVoxels;
  bool multiResolution;
  bool bayesian;
};

#endif
