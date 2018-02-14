#ifndef VOLUME_INSTANCE_HPP
#define VOLUME_INSTANCE_HPP
#include "volume_template.hpp"
#include "volume_traits.hpp"
#include "octree.hpp"

/******************************************************************************
 *
 * Volume Class Instance
 * 
******************************************************************************/


#ifndef STORAGE_LAYOUT
#define STORAGE_LAYOUT DynamicStorage
#endif

#ifndef INDEX_STRUCTURE
#define INDEX_STRUCTURE Octree
#endif

typedef BFusion FieldType;

template <typename T>
using Volume = VolumeTemplate<T, STORAGE_LAYOUT, INDEX_STRUCTURE>;

#endif
