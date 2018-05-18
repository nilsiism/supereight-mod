#ifndef SE_COMMON_H
#define SE_COMMON_H

#include <cmath>

#define SOPHUS_DISABLE_ENSURES

/* 
 * When compiling in debug mode Eigen compilation fails 
 * due to -Wunused-parameter. Disable it if compiling with GCC.
 */
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include <Eigen/Dense>
#include <sophus/se3.hpp>
#pragma GCC diagnostic pop
#else
#include <Eigen/Dense>
#include <sophus/se3.hpp>
#endif

#include "se_math_helper.h"

#endif
