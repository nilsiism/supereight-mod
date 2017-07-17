#include <math_utils.h>
#include "continuous/volume.hpp"


/* Raycasting implementations */ 
#include "bfusion/rendering_impl.hpp"
#include "kfusion/rendering_impl.hpp"

template<typename T>
void raycastKernel(const Volume<T>& volume, float3* vertex, float3* normal, uint2 inputSize, 
    const Matrix4 view, const float nearPlane, const float farPlane, 
    const float mu, const float step, const float largestep) {
	TICK();
	unsigned int y;
#pragma omp parallel for shared(normal, vertex), private(y)
	for (y = 0; y < inputSize.y; y++)
#pragma simd
		for (unsigned int x = 0; x < inputSize.x; x++) {

			uint2 pos = make_uint2(x, y);
      ray_iterator<typename Volume<T>::field_type> ray(volume._map_index, get_translation(view), 
          normalize(rotate(view, make_float3(x, y, 1.f))), nearPlane, farPlane);
      const float t_min = ray.next(); /* Get distance to the first intersected block */
      const float4 hit = t_min > 0.f ? 
        raycast(volume, pos, view, t_min*volume._dim/volume._size, 
          farPlane, mu, step, largestep) : make_float4(0.f);
			if (hit.w > 0.0) {
				vertex[pos.x + pos.y * inputSize.x] = make_float3(hit);
				float3 surfNorm = volume.grad(make_float3(hit), 
            [](const auto& val){ return val.x; });
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

template <typename T>
void renderVolumeKernel(const Volume<T>& volume, uchar4* out, const uint2 depthSize, const Matrix4 view, 
    const float nearPlane, const float farPlane, const float mu,
		const float step, const float largestep, const float3 light,
		const float3 ambient) {
	TICK();
	unsigned int y;
#pragma omp parallel for shared(out), private(y)
	for (y = 0; y < depthSize.y; y++) {
		for (unsigned int x = 0; x < depthSize.x; x++) {
			const uint2 pos = make_uint2(x, y);
      ray_iterator<typename Volume<T>::field_type> ray(volume._map_index, get_translation(view), 
          normalize(rotate(view, make_float3(x, y, 1.f))), nearPlane, farPlane);
      const float t_min = ray.next(); /* Get distance to the first intersected block */
      const float4 hit = t_min > 0.f ? 
        raycast(volume, pos, view, t_min*volume._dim/volume._size, 
          farPlane, mu, step, largestep) : make_float4(0.f);
			if (hit.w > 0) {
				const float3 test = make_float3(hit);
				const float3 surfNorm = volume.grad(test, [](const auto& val){ return val.x; });
				if (length(surfNorm) > 0) {
					const float3 diff = (std::is_same<T, SDF>::value ?
              normalize(light - test) : normalize(test - light));
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

