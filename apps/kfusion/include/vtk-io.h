#ifndef VTK_IO_H
#define VTK_IO_H
#include <lodepng.h>
#include <sstream>
#include <math_utils.h>

//http://stackoverflow.com/questions/236129/split-a-string-in-c
static inline void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

template <typename T>
void savePointCloud(const T* in, const int num_points, 
    const char* filename){
  std::stringstream points;

  for(int i = 0; i < num_points; ++i ){
    points << in[i].x << " " << in[i].y << " " << in[i].z << std::endl; 
  }   

  std::ofstream f;
  f.open(std::string(filename).append(".vtk").c_str());
  f << "# vtk DataFile Version 1.0" << std::endl;
  f << "vtk mesh generated from KFusion" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET POLYDATA" << std::endl;

  f << "POINTS " << num_points << " FLOAT" << std::endl;
  f << points.str();
  f.close();
}

// 
// Function to save a slice of the 3D volume to vtk file. The convention
// for the iteration limits is that it expets them in the format [lower, upper),
// i.e. the upper bound is not included in the set. 
//
template <typename T>
void save3DSlice(const T* in, const uint3 lower, const uint3 upper, 
    const uint3 size, const char* filename){
  std::stringstream x_coordinates, y_coordinates, z_coordinates, scalars;
  std::ofstream f;
  f.open(filename);
 
  const int dimX = upper.x - lower.x;
  const int dimY = upper.y - lower.y;
  const int dimZ = upper.z - lower.z;

  f << "# vtk DataFile Version 1.0" << std::endl;
  f << "vtk mesh generated from KFusion" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET RECTILINEAR_GRID" << std::endl;
  f << "DIMENSIONS " << dimX << " " << dimY << " " << dimZ << std::endl;

  for(int x = lower.x; x < upper.x; ++x)
    x_coordinates << x << " ";  
  for(int y = lower.y; y < upper.y; ++y)
    y_coordinates << y << " ";  
  for(int z = lower.z; z < upper.z; ++z)
    z_coordinates << z << " ";  

  for(int z = lower.z; z < upper.z; ++z)
    for(int y = lower.y; y < upper.y; ++y)
      for(int x = lower.x; x < upper.x; ++x) {
        const float data = in[x + y*size.x + z*size.x*size.y].x * 0.00003051944088f;
        scalars << data  << std::endl;
      }

  f << "X_COORDINATES " << dimX << " int " << std::endl;
  f << x_coordinates.str() << std::endl;

  f << "Y_COORDINATES " << dimY << " int " << std::endl;
  f << y_coordinates.str() << std::endl;

  f << "Z_COORDINATES " << dimZ << " int " << std::endl;
  f << z_coordinates.str() << std::endl;

  f << "POINT_DATA " << dimX*dimY*dimZ << std::endl;
  f << "SCALARS scalars float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  f << scalars.str() << std::endl;
  f.close();
} 

template <typename MapType>
void save3DSlice(const MapType& in, const int3 lower, const int3 upper, 
    const int3, const char* filename){
  std::stringstream x_coordinates, y_coordinates, z_coordinates, scalars;
  std::ofstream f;
  f.open(filename);
 
  const int dimX = upper.x - lower.x;
  const int dimY = upper.y - lower.y;
  const int dimZ = upper.z - lower.z;

  f << "# vtk DataFile Version 1.0" << std::endl;
  f << "vtk mesh generated from KFusion" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET RECTILINEAR_GRID" << std::endl;
  f << "DIMENSIONS " << dimX << " " << dimY << " " << dimZ << std::endl;

  for(int x = lower.x; x < upper.x; ++x)
    x_coordinates << x << " ";  
  for(int y = lower.y; y < upper.y; ++y)
    y_coordinates << y << " ";  
  for(int z = lower.z; z < upper.z; ++z)
    z_coordinates << z << " ";  

  for(int z = lower.z; z < upper.z; ++z)
    for(int y = lower.y; y < upper.y; ++y)
      for(int x = lower.x; x < upper.x; ++x) {
        const float data = in[make_uint3(x, y, z)].x;
        scalars << data  << std::endl;
      }

  f << "X_COORDINATES " << dimX << " int " << std::endl;
  f << x_coordinates.str() << std::endl;

  f << "Y_COORDINATES " << dimY << " int " << std::endl;
  f << y_coordinates.str() << std::endl;

  f << "Z_COORDINATES " << dimZ << " int " << std::endl;
  f << z_coordinates.str() << std::endl;

  f << "POINT_DATA " << dimX*dimY*dimZ << std::endl;
  f << "SCALARS scalars float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  f << scalars.str() << std::endl;
  f.close();
} 

void printNormals(const float3* in, const unsigned int xdim, 
                 const unsigned int ydim, const char* filename) {
  unsigned char* image = new unsigned char [xdim * ydim * 4];
  for(unsigned int y = 0; y < ydim; ++y)
    for(unsigned int x = 0; x < xdim; ++x){
      const float3 n = in[x + y*xdim];
      image[4 * xdim * y + 4 * x + 0] = (n.x/2 + 0.5) * 255; 
      image[4 * xdim * y + 4 * x + 1] = (n.y/2 + 0.5) * 255;
      image[4 * xdim * y + 4 * x + 2] = (n.z/2 + 0.5) * 255;
      image[4 * xdim * y + 4 * x + 3] = 255;
    }
  lodepng_encode32_file(std::string(filename).append(".png").c_str(),
      image, xdim, ydim);
} 

void parseGTFile(const std::string& filename, std::vector<Matrix4>& poses){

  std::ifstream file;
  file.open(filename.c_str());

  if(!file.is_open()) {
    std::cout << "Failed to open GT file " << filename << std::endl;
    return;
  }

  std::string line;
  while (getline(file,line))
  {
    std::vector<std::string> data;
    split(line, ' ', data);

    if (data.size() != 8)
    {
      continue;
    }
    const float3 trans = make_float3(std::stof(data[1]), std::stof(data[2]),
        std::stof(data[3]));
    const float4 quat = make_float4(std::stof(data[4]), std::stof(data[5]),
        std::stof(data[6]), std::stof(data[7]));
    poses.push_back(toMatrix4(quat, trans));
  }
}

#endif
