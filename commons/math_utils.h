#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#ifndef CUDA
#include <vector_types.h>
#else
#include <cuda_runtime.h>
#endif
#include <cutil_math.h>


typedef struct sMatrix4 {
	float4 data[4];
} Matrix4;

inline __host__ __device__ float3 get_translation(const Matrix4& view) {
	return make_float3(view.data[0].w, view.data[1].w, view.data[2].w);
}

////////////////////////////////////////////////////////////////////////////////
// ceilf - missing from cutil_math.h
////////////////////////////////////////////////////////////////////////////////

inline __host__     __device__ float2 ceilf(float2 v) {
	return make_float2(ceilf(v.x), ceilf(v.y));
}
inline __host__     __device__ float3 ceilf(float3 v) {
	return make_float3(ceilf(v.x), ceilf(v.y), ceilf(v.z));
}
inline __host__     __device__ float4 ceilf(float4 v) {
	return make_float4(ceilf(v.x), ceilf(v.y), ceilf(v.z), ceilf(v.w));
}

inline __host__ __device__ bool operator==(const float3 a, float b){
  return((a.x == b) && (a.y == b) && (a.z == b));
}

inline __host__ __device__ bool in(const unsigned int value, const unsigned int lower, 
               const unsigned int upper){
  return value >= lower && value <= upper;
}

inline __host__     __device__ uchar3 operator*(const uchar3 a, float v) {
	return make_uchar3(a.x * v, a.y * v, a.z * v);
}

inline float4 operator*(const Matrix4 & M, const float4 & v) {
	return make_float4(dot(M.data[0], v), dot(M.data[1], v), dot(M.data[2], v),
			dot(M.data[3], v));
}

inline bool operator<=(const uint3 a, const int3 b){
  return((a.x <= b.x) && (a.y <= b.y) && (a.z <= b.z));
}

static inline float __host__ __device__ bspline(float t){
  float value = 0.f;
  if(t >= -3.0f && t <= -1.0f) {
    value = std::pow((3 + t), 3)/48.0f;   
  } else if( t > -1 && t <= 1) {
    value = 0.5f + (t*(3 + t)*(3 - t))/24.f;
  } else if(t > 1 && t <= 3){
    value = 1 - std::pow((3 - t), 3)/48.f;
  } else if(t > 3) {
    value = 1.f;
  }
  return value;
}

static inline __host__ __device__ float H(const float val){
  const float Q_1 = bspline(val);
  const float Q_2 = bspline(val - 3);
  return Q_1 - Q_2 * 0.5f;
}

inline Matrix4 outer(const float4& a, const float4& b){
  Matrix4 mat;
  mat.data[0] = make_float4(a.x*b.x, a.x*b.y, a.x*b.z, a.x*b.w);
  mat.data[1] = make_float4(a.y*b.x, a.y*b.y, a.y*b.z, a.y*b.w);
  mat.data[2] = make_float4(a.z*b.x, a.z*b.y, a.z*b.z, a.z*b.w);
  mat.data[3] = make_float4(a.w*b.x, a.w*b.y, a.w*b.z, a.w*b.w);
  return mat;
}

// Converting quaternion and trans to SE3 matrix 
// Following the implementation provided in TUM scripts.
inline Matrix4 toMatrix4(float4 quat, const float3& trans) {
  const float n = dot(quat, quat);
  quat = quat*(sqrtf(2.0/n));
  Matrix4 mat = outer(quat, quat);
  Matrix4 se3_mat;
  se3_mat.data[0] = make_float4(1.0-mat.data[1].y - mat.data[2].z, mat.data[0].y - mat.data[2].w,     mat.data[0].z + mat.data[1].w, trans.x);
  se3_mat.data[1] = make_float4(mat.data[0].y + mat.data[2].w, 1.0-mat.data[0].x - mat.data[2].z,     mat.data[1].z - mat.data[0].w, trans.y);
  se3_mat.data[2] = make_float4(mat.data[0].z - mat.data[1].w,     mat.data[1].z + mat.data[0].w, 1.0-mat.data[0].x - mat.data[1].y, trans.z);
  se3_mat.data[3] = make_float4(0.0, 0.0, 0.0, 1.0);
  return se3_mat;
}

#endif
