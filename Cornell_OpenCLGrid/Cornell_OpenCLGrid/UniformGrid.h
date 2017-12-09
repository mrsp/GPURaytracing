//
//  UniformGrid.h
//  
//
//  Created by Stylianos Piperakis on 5/24/15.
//  Copyright (c) 2015 Stylianos Piperakis. All rights reserved.
//

#ifndef __UniformGrid__
#define __UniformGrid__

#include "primitives.h"


/* Voxel Struct */
typedef struct Voxel_
{
        cl_float3 Max;
        cl_float3 Min;
        cl_float3 Max_;
        cl_float3 Min_;  //For comparison reasons
        int MaxTriangleLimit;
        int MaxSphereLimit;

        Triangle** Triangles;
        Sphere** Spheres;
        int TrianglesCount;
        int SpheresCount;
        float Size;

         Voxel_();
    
       /*Function Prototypes */
        void Delete();
        void SetAttributes(const cl_float3 &Min, float Size);
        bool checkInside(const cl_float3 &Point);
        void AddTriangle(Triangle *Triangle);
        void AddSphere(Sphere *Sphere);
        bool IntersectTriangle(Triangle *Triangle);
        bool IntersectSphere(Sphere *Sphere);

    
    bool triBoxOverlap(cl_float3 boxcenter,float boxhalfsize, Triangle* Triangle);
    bool AXISTEST_Z0(cl_float3 V0, cl_float3 V1, float a, float b, float fa, float fb);
    bool AXISTEST_Z12(cl_float3 V1, cl_float3 V2, float a, float b, float fa, float fb);

    bool AXISTEST_Y1(cl_float3 V0, cl_float3 V1, float a, float b, float fa, float fb);
    bool AXISTEST_Y02(cl_float3 V0, cl_float3 V2, float a, float b, float fa, float fb);
    bool AXISTEST_X2(cl_float3 V0, cl_float3 V1, float a, float b, float fa, float fb);
    bool AXISTEST_X01(cl_float3 V0, cl_float3 V2, float a, float b, float fa, float fb);
    bool planeBoxOverlap(cl_float3 Normal, cl_float3 Vertex, cl_float3 maxbox);
    

}Voxel;

/* Grid Struct */
typedef struct UniformGrid_
{
    
        cl_float3 d;
        int k;
        float Volume;
        cl_float3 SphereVertex[8]; // For bounding box of sphere

        cl_float3 Min,Max;
        Voxel *Voxels;
    
        Triangle *triag;
        Sphere *spheres;
    
        float VoxelSize;
        int numPrim, numTriangle, numSpheres;
        int Nx,Ny,Nz;
        int Nxy,Nxyz, Nxmin, Nymin, Nzmin;

    /* Prototypes */
    UniformGrid_();
    void Delete();
    void computeAABB(Triangle *mtriag, int numTriangle_,  Sphere *mspheres, int numSphere_);
    void computeDimensions(float ratio);
    void computeSphereVertices(Sphere sph);
    
  

    
    
    
}UniformGrid;




/*---- World/Voxel Transformation ----*/
float VoxelToWorldX(int x, float VoxelSize, float MinX);

float VoxelToWorldY(int y, float VoxelSize, float MinY);

float VoxelToWorldZ(int z, float VoxelSize, float MinZ);

cl_float3 VoxelToWorld(const cl_int3 &Voxel ,float VoxelSize, cl_float3 Min);

int WorldToVoxelX(float x, float VoxelSize, float MinX);

int WorldToVoxelY(float y, float VoxelSize, float MinY);

int WorldToVoxelZ(float z, float VoxelSize, float MinZ);

cl_int3 WorldToVoxel(const cl_float3 &World, float VoxelSize, cl_float3 Min);





//Triangle* Traverse(const cl_float3i &voxel, const cl_float3 &pos, const Ray &r, float &intersected);
bool Traverse(UniformGrid Grid,const cl_int3 &voxel, const cl_float3 &pos, const Ray &r, float &intersected, Triangle* &triag, Sphere* &sph);
float RayTriangleIntersect(Ray r, Triangle *tri);
float RaySphereIntersect(Ray r, Sphere sph);














#endif 