#ifndef VOLUME_INSTANCE_HPP
#define VOLUME_INSTANCE_HPP
#include <se/octree.hpp>
#include "volume_template.hpp"
#include "volume_traits.hpp"

/******************************************************************************
 *
 * Volume Class Instance
 * 
******************************************************************************/

#ifndef INDEX_STRUCTURE
#define INDEX_STRUCTURE Octree
#endif

typedef SDF FieldType;

template <typename T>
using Volume = VolumeTemplate<T, INDEX_STRUCTURE>;

#endif
