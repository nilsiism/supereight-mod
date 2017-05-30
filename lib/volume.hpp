#ifndef VOLUME_H
#define VOLUME_H

#include "volume_template.hpp"
#include "octree.hpp"
#include "voxel_traits.hpp"

/**
 * Default typedefs for the Volume Class
 * */

#ifndef FIELD_TYPE
#define FIELD_TYPE SDF
#endif

#ifndef STORAGE_LAYOUT
#define STORAGE_LAYOUT DynamicStorage
#endif

#ifndef INDEX_STRUCTURE
#define INDEX_STRUCTURE Octree
#endif

typedef VolumeTemplate<FIELD_TYPE, STORAGE_LAYOUT, INDEX_STRUCTURE> Volume;
#endif
