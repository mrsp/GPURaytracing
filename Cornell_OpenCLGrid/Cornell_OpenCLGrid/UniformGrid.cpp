//
//  UniformGrid.cpp
//  
//
//  Created by Stylianos Piperakis on 5/24/15.
//  Copyright (c) 2015 Stylianos Piperakis. All rights reserved.
//

#include "UniformGrid.h"


/* ----------------------------------------------------------------------------------------------------*/

Voxel::Voxel_()
{
    Triangles = NULL;
    Spheres = NULL;
    
    Size = 0.0;
    TrianglesCount =0;
    MaxTriangleLimit=1;
    SpheresCount = 0;
    MaxSphereLimit=1;
    
    
    
    
}

/* ----------------------------------------------------------------------------------------------------*/


void Voxel::Delete()
{
    
    if(Triangles !=NULL)
    {
        delete[] Triangles;
        Triangles = NULL;
        TrianglesCount = 0;
        MaxTriangleLimit = 1;
        
    }
    
    if(Spheres !=NULL)
    {
        delete[] Spheres;
        Spheres = NULL;
        SpheresCount =0;
        MaxSphereLimit= 1;
       
    }
    Size = 0.0;
    Min = Min_ = Max = Max_ = setVec(0.0);
    

}

/* ----------------------------------------------------------------------------------------------------*/

void Voxel::SetAttributes(const cl_float3 &Min, float Size)
{
    this->Size = Size;
    this->Min = Min;
    this->Max = this->Min;
    this->Max.x+= Size;
    this->Max.y+= Size;
    this->Max.z+= Size;

    
    
    this->Min_ = this->Min;
    this->Min_.x -=0.001;
    this->Min_.y -=0.001;
    this->Min_.z -=0.001;

    
    this->Max_ = this->Max;
    this->Max_.x +=0.001;
    this->Max_.y +=0.001;
    this->Max_.z +=0.001;

}




bool Voxel::checkInside(const cl_float3 &Point)
{
    
    if (Min_.x < Point.x && Point.x < Max_.x)
    {
        if (Min_.y < Point.y && Point.y < Max_.y)
        {
            
            if (Min_.z < Point.z && Point.z < Max_.z)
            {
                
                return true;
            }
        }
        
    }
    return false;
}

void Voxel::AddTriangle(Triangle* Triangle)
{
    if (TrianglesCount % MaxTriangleLimit == 0)
    {
       
        class Triangle** Triangles_ =  Triangles;
        
        MaxTriangleLimit *= 2;
        
        
       // if(Triangles !=NULL) delete[] Triangles;
        Triangles = new class Triangle*[MaxTriangleLimit];
        
        for (int i = 0; i < TrianglesCount; i++)
        {
            Triangles[i] = Triangles_[i];
        }

        if (Triangles_ !=NULL) delete[] Triangles_;
        
    }
    Triangles[TrianglesCount] = Triangle;
    TrianglesCount++;
    
}
/* ----------------------------------------------------------------------------------------------------*/
void Voxel::AddSphere(Sphere* Sphere)
{
    if (SpheresCount % MaxSphereLimit == 0)
    {
        
        class Sphere** Spheres_ =  Spheres;
        
        MaxSphereLimit *= 2;
        
        
        Spheres = new class Sphere*[MaxSphereLimit];
        
        for (int i = 0; i < SpheresCount; i++)
        {
            Spheres[i] = Spheres_[i];
        }
        
        if (Spheres_ !=NULL) delete[] Spheres_;
        
    }
    Spheres[SpheresCount] = Sphere;
    SpheresCount++;
    
}

/* ----------------------------------------------------------------------------------------------------*/

bool Voxel::IntersectTriangle(Triangle *Triangle)
{
    cl_float3 Vec;
    
    Vec=Min;
    Vec.x += (Max.x - Min.x)/2.0;
    Vec.y += (Max.y - Min.y)/2.0;
    Vec.z += (Max.z - Min.z)/2.0;

      if(triBoxOverlap(Vec, Size/2.0, Triangle)) return true;
    return false;
    
    
    
}
bool Voxel::IntersectSphere(Sphere *Sphere)
{
    float d=0;
    
    float e1=Sphere->p.x - Min.x;
    float e2=Sphere->p.x - Max.x;
    if (e1<0){
        if(e1 < - Sphere->rad) return false;
        d+=e1*e1;
    }
    else if (e2 >0){
        if (e2>Sphere->rad) return false;
        d+=e2*e2;
    }
    
     e1=Sphere->p.y - Min.y;
     e2=Sphere->p.y - Max.y;
    if (e1<0){
        if(e1 < - Sphere->rad) return false;
        d+=e1*e1;
    }
    else if (e2 >0){
        if (e2>Sphere->rad)  return false;
        d+=e2*e2;
    }

    
     e1=Sphere->p.z - Min.z;
     e2=Sphere->p.z - Max.z;
    if (e1<0){
        if(e1 < - Sphere->rad) return false;
        d+=e1*e1;
    }
    else if (e2 >0){
        if (e2>Sphere->rad) return false;
        d+=e2*e2;
    }

    
    
    
    if(d<=Sphere->rad*Sphere->rad) return true;
    return false;
}

/* ----------------------------------------------------------------------------------------------------*/
bool Voxel::planeBoxOverlap(cl_float3 Normal, cl_float3 Vertex, cl_float3 maxbox)
{
    

  
    cl_float3 vmin,vmax;
    
    /* X axis */
        if (Normal.x >0.0)
        {
            vmin.x = -maxbox.x -Vertex.x;
            vmax.x =  maxbox.x -Vertex.x;
            
        }
        
        else
        {
            vmin.x = maxbox.x - Vertex.x;
            vmax.x = -maxbox .x - Vertex.x;
            
        }
     /* Y axis */
        if (Normal.y >0.0)
        {
            vmin.y = -maxbox.y -Vertex.y;
            vmax.y =  maxbox.y -Vertex.y;
        
        }
    
        else
        {
            vmin.y = maxbox.y - Vertex.y;
            vmax.y = -maxbox.y - Vertex.y;
        
        }
    
    /* Z axis */
        if (Normal.z > 0.0)
        {
            vmin.z = -maxbox.z -Vertex.z;
            vmax.z =  maxbox.z -Vertex.z;
        
        }
    
        else
        {
            vmin.z = maxbox.z - Vertex.z;
            vmax.z = -maxbox.z - Vertex.z;
        
        }
    
    
    if(dot(Normal,vmin)>0.0) return false;
    if(dot(Normal,vmax)>=0.0) return true;

    return false;
}

/* -------------------------------- X - Tests ------------------------------------*/
bool Voxel::AXISTEST_X01(cl_float3 V0, cl_float3 V2, float a, float b, float fa, float fb)
{
    
    float p0 = a * V0.y - b * V0.z;
    
    float p2 = a * V2.y - b * V2.z;
    
    float min,max;

    if (p0 < p2) {min=p0; max=p2;} else {min=p2;max=p0;}
    
    float rad = fa * Size/2.0 + fb * Size/2.0;
    
    if(min>rad || max<-rad) return false;
    return true;
    
}

bool Voxel::AXISTEST_X2(cl_float3 V0, cl_float3 V1, float a, float b, float fa, float fb)
{
    
    float p0 = a * V0.y - b * V0.z;
    
    float p1 = a * V1.y - b * V1.z;
    
    float min,max;
    
    if (p0 < p1) {min=p0; max=p1;} else {min=p1;max=p0;}
    
    float rad = fa * Size/2.0 + fb * Size/2.0;
    
    if(min>rad || max<-rad) return false;
    return true;
    
}
/* -------------------------------- Y - Tests ------------------------------------*/

bool Voxel::AXISTEST_Y02(cl_float3 V0, cl_float3 V2, float a, float b, float fa, float fb)
{
    
    float p0 = -a * V0.x + b * V0.z;
    
    float p2 = -a * V2.x + b * V2.z;
    
    float min,max;
    
    if (p0 < p2) {min=p0; max=p2;} else {min=p2;max=p0;}
    
    float rad = fa * Size/2.0 + fb * Size/2.0;
    
    if(min>rad || max<-rad) return false;
    return true;
    
}

bool Voxel::AXISTEST_Y1(cl_float3 V0, cl_float3 V1, float a, float b, float fa, float fb)
{
    
    float p0 = -a * V0.x + b * V0.z;
    
    float p1 = -a * V1.x + b * V1.z;
    
    float min,max;
    
    if (p0 < p1) {min=p0; max=p1;} else {min=p1;max=p0;}
    
    float rad = fa * Size/2.0 + fb * Size/2.0;
    
    if(min>rad || max<-rad) return false;
    return true;
    
}

/* -------------------------------- Z - Tests ------------------------------------*/

bool Voxel::AXISTEST_Z12(cl_float3 V1, cl_float3 V2, float a, float b, float fa, float fb)
{
    
    float p1 = a * V1.x - b * V1.y;
    
    float p2 = a * V2.x - b * V2.y;
    
    float min,max;
    
    if (p2 < p1) {min=p2; max=p1;} else {min=p1;max=p2;}
    
    float rad = fa * Size/2.0 + fb * Size/2.0;
    
    if(min>rad || max<-rad) return false;
    return true;
    
}


bool Voxel::AXISTEST_Z0(cl_float3 V0, cl_float3 V1, float a, float b, float fa, float fb)
{
    
    float p0 = a * V0.x - b * V0.y;
    
    float p1 = a * V1.x - b * V1.y;
    
    float min,max;
    
    if (p0 < p1) {min=p0; max=p1;} else {min=p1;max=p0;}
    
    float rad = fa * Size/2.0 + fb * Size/2.0;
    
    if(min>rad || max<-rad) return false;
    return true;
    
}
/* -------------------------------- ------- ------------------------------------*/


bool Voxel::triBoxOverlap(cl_float3 boxcenter,float boxhalfsize, Triangle * Triangle)// triverts[3][3])

{
    
    
    
    /*    use separating axis theorem to test overlap between triangle and box */
    
    /*    need to test for overlap in these directions: */
    
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    
    /*       we do not even need to test these) */
    
    /*    2) normal of the triangle */
    
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    
    /*       this gives 3x3=9 more tests */
    
    cl_float3 v0,v1,v2;
    
    //   float axis[3];
    
    float min,max,fex,fey,fez;		// -NJMP- "d" local variable removed
    
    cl_float3 normal,e0,e1,e2;
    
    
    
    /* This is the fastest branch on Sun */
    
    /* move everything so that the boxcenter is in (0,0,0) */
    
    v0 = Triangle->Vertex[0];
    v0.x -= boxcenter.x;
    v0.y -= boxcenter.y;
    v0.z -= boxcenter.z;

    
    v1 = Triangle->Vertex[1];
    v1.x -= boxcenter.x;
    v1.y -= boxcenter.y;
    v1.z -= boxcenter.z;

    
    v2 = Triangle->Vertex[2];
    v2.x -= boxcenter.x;
    v2.y -= boxcenter.y;
    v2.z -= boxcenter.z;
    
    
    
    /* compute triangle edges */
    
    
    
    
    e0 = v1;        /* tri edge 0 */
    e0.x-=v0.x;
    e0.y-=v0.y;
    e0.z-=v0.z;

    e1 = v2;
    e1.x-=v1.x;      /* tri edge 1 */
    e1.y-=v1.y;
    e1.z-=v1.z;
    
   
    e2 = v0;       /* tri edge 2 */
    e2.x-=v2.x;
    e2.y-=v2.y;
    e2.z-=v2.z;

    
    
    /* Bullet 3:  */
    
    /*  test the 9 tests first (this was faster) */
    
    fex = fabsf(e0.x);
    
    fey = fabsf(e0.y);
    
    fez = fabsf(e0.z);
    
    if(!AXISTEST_X01(v0,v2,e0.z, e0.y, fez, fey)) return false;
    
    if(!AXISTEST_Y02(v0,v2,e0.z, e0.x, fez, fex)) return false;
    
    if(!AXISTEST_Z12(v1,v2,e0.y, e0.x, fey, fex)) return false;
    
    
    
    fex = fabsf(e1.x);
    
    fey = fabsf(e1.y);
    
    fez = fabsf(e1.z);
    
     if(!AXISTEST_X01(v0,v2, e1.z, e1.y, fez, fey)) return false;
    
     if(!AXISTEST_Y02(v0,v2, e1.z, e1.x, fez, fex)) return false;
     if(!AXISTEST_Z0(v0, v1, e1.y, e1.x, fey, fex)) return false;
    
    
    
    fex = fabsf(e2.x);
    
    fey = fabsf(e2.y);
    
    fez = fabsf(e2.z);
    
     if(!AXISTEST_X2(v0,v1,e2.z, e2.y, fez, fey)) return false;
    
     if(!AXISTEST_Y1(v0,v1,e2.z, e2.x, fez, fex)) return false;
    
     if(!AXISTEST_Z12(v1,v2,e2.y, e2.x, fey, fex)) return false;
    
    
    
    /* Bullet 1: */
    
    /*  first test overlap in the {x,y,z}-directions */
    
    /*  find min, max of the triangle each direction, and test for overlap in */
    
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    
    /*  the triangle against the AABB */
    
    
    
    /* test in X-direction */
    min = max  = v0.x;
    
    if(v1.x<min) min=v1.x;
    if(v1.x>max) max=v1.x;
    if(v2.x<min) min=v2.x;
    if(v2.x>max) max=v2.x;
    
    
    if(min>boxhalfsize || max < -boxhalfsize) return false;
    
    
    
    /* test in Y-direction */
    min = max  = v0.y;
    
    if(v1.y<min) min=v1.y;
    if(v1.y>max) max=v1.y;
    if(v2.y<min) min=v2.y;
    if(v2.y>max) max=v2.y;

    
    if(min > boxhalfsize || max < -boxhalfsize) return false;
    
    
    
    /* test in Z-direction */
    min = max  = v0.x;
    
    if(v1.z<min) min=v1.z;
    if(v1.z>max) max=v1.z;
    if(v2.z<min) min=v2.z;
    if(v2.z>max) max=v2.z;

    
    if(min>boxhalfsize || max < -boxhalfsize) return false;
    
    
    
    /* Bullet 2: */
    
    /*  test if the box intersects the plane of the triangle */
    
    /*  compute plane equation of triangle: normal*x+d=0 */
    
    normal=cross(e0,e1);
    
    // -NJMP- (line removed here)
    
    if(!planeBoxOverlap(normal,v0,setVec(boxhalfsize,boxhalfsize,boxhalfsize))) return false;	// -NJMP-
    
    
    
    return true;   /* box and triangle overlaps */
    
}



/* ----------------------------------------------------------------------------------------------------*/
UniformGrid::UniformGrid_()
{
    
    /* For AABB Computation */
    Nx = Ny = Nz = Nxmin = Nymin = Nzmin= Nxy= Nxyz = 0;
    k = 5;
    Min.x=1e20;
    Min.y=1e20;
    Min.z=1e20;
    Max.x = -Min.x;
    Max.y = -Min.y;
    Max.z = -Min.z;

    Voxels = NULL;
    VoxelSize= 0.0;
    
    
}
/*----------------------------------------------------------------------------------------------------*/


/* ----------------------------------------------------------------------------------------------------*/

void UniformGrid::Delete()
{
    
    if(Voxels!= NULL)
    {    for (int i=0; i<Nxyz ; i++)
        Voxels[i].Delete();
        
        Min=Max=setVec(0.0);
        Nx = Ny = Nz = Nxmin = Nymin = Nzmin= Nxy= Nxyz = 0;
        VoxelSize=0.0;
        delete [] Voxels;
        Voxels = NULL;
    }

    
}
/* ----------------------------------------------------------------------------------------------------*/
void UniformGrid::computeAABB(Triangle *triag, int numTriangle,  Sphere *spheres, int numSphere)
{
   // Delete();
    
    /* Compute Primitives */
    this->triag = triag;
    this->spheres = spheres;
    
    
    this->numTriangle =  numTriangle;
    this->numSpheres =  numSphere;

    numPrim = numTriangle+ numSpheres;
    
   /* Bounding box */
    
    
    
    for (unsigned int i=0; i<numSpheres; i++)
    {
        computeSphereVertices(spheres[i]);
        for (unsigned int j=0; j<8; j++)
        {
            if(SphereVertex[j].x< Min.x) Min.x = SphereVertex[j].x;
            if(SphereVertex[j].y< Min.y) Min.y = SphereVertex[j].y;
            if(SphereVertex[j].z< Min.z) Min.z = SphereVertex[j].z;
            
            if(SphereVertex[j].x > Max.x) Max.x = SphereVertex[j].x;
            if(SphereVertex[j].y > Max.y) Max.y = SphereVertex[j].y;
            if(SphereVertex[j].z > Max.z) Max.z = SphereVertex[j].z;
        }
        
    }
    
  
  
    
    for (unsigned int i =0; i< numTriangle; i++)
    {
        
        
        
        if(triag[i].Vertex[0].x< Min.x) Min.x = triag[i].Vertex[0].x;
        if(triag[i].Vertex[0].y< Min.y) Min.y = triag[i].Vertex[0].y;
        if(triag[i].Vertex[0].z< Min.z) Min.z = triag[i].Vertex[0].z;
        
        if(triag[i].Vertex[0].x > Max.x) Max.x = triag[i].Vertex[0].x;
        if(triag[i].Vertex[0].y > Max.y) Max.y = triag[i].Vertex[0].y;
        if(triag[i].Vertex[0].z > Max.z) Max.z = triag[i].Vertex[0].z;

        
        
        if(triag[i].Vertex[1].x< Min.x) Min.x = triag[i].Vertex[1].x;
        if(triag[i].Vertex[1].y< Min.y) Min.y = triag[i].Vertex[1].y;
        if(triag[i].Vertex[1].z< Min.z) Min.z = triag[i].Vertex[1].z;
        
        if(triag[i].Vertex[2].x< Min.x) Min.x = triag[i].Vertex[2].x;
        if(triag[i].Vertex[2].y< Min.y) Min.y = triag[i].Vertex[2].y;
        if(triag[i].Vertex[2].z< Min.z) Min.z = triag[i].Vertex[2].z;
        
        
        if(triag[i].Vertex[1].x > Max.x) Max.x = triag[i].Vertex[1].x;
        if(triag[i].Vertex[1].y > Max.y) Max.y = triag[i].Vertex[1].y;
        if(triag[i].Vertex[1].z > Max.z) Max.z = triag[i].Vertex[1].z;
        
        
        if(triag[i].Vertex[2].x > Max.x) Max.x = triag[i].Vertex[2].x;
        if(triag[i].Vertex[2].y > Max.y) Max.y = triag[i].Vertex[2].y;
        if(triag[i].Vertex[2].z > Max.z) Max.z = triag[i].Vertex[2].z;
    }

    
}



void UniformGrid::computeDimensions(float ratio)
{
    d.x= (Max.x-Min.x);
    d.y= (Max.y-Min.y);
    d.z= (Max.z-Min.z);
    Volume = d.x * d.y * d.z;
    VoxelSize = pow( k * numPrim/Volume, -1/3.)/ratio;
    std::cout<< " Voxel Size " << VoxelSize<< std::endl;
    
    
    Min.x /= VoxelSize;
    Min.y /= VoxelSize;
    Min.z /= VoxelSize;

    Max.x /= VoxelSize;
    Max.y /= VoxelSize;
    Max.z /= VoxelSize;
    

    
    
    Min.x = floor(Min.x); Max.x = ceil(Max.x);
    Min.y = floor(Min.y); Max.y = ceil(Max.y);
    Min.z = floor(Min.z); Max.z = ceil(Max.z);
    
    if(Min.x == Max.x) Max.x += 1.0;
    if(Min.y == Max.y) Max.y += 1.0;
    if(Min.z == Max.z) Max.z += 1.0;
    
    

    
    Nx = (int) (Max.x-Min.x);
    Ny = (int) (Max.y-Min.y);
    Nz = (int) (Max.z-Min.z);
    Nxmin= Nx-1;
    Nymin= Ny-1;
    Nzmin= Nz-1;
    
    
    Min.x *= VoxelSize;
    Min.y *= VoxelSize;
    Min.z *= VoxelSize;
    
    Max.x *= VoxelSize;
    Max.y *= VoxelSize;
    Max.z *= VoxelSize;
    
    Nxy= Nx * Ny;
    Nxyz=Nxy * Nz;
    
    std::cout<< " Grid Volume : " << Volume << std::endl;
    std::cout<< " Dimension X: " << Nx << " Dimension Y: " << Ny << " Dimension Z: "<< Nz <<std::endl;
    std::cout<<" Grid Initialized " <<std::endl;

    Voxels= new Voxel[Nxyz];
    
    for(int z = 0; z < Nz; z++)
    {
        for(int y = 0; y < Ny; y++)
        {
            for(int x = 0; x < Nx; x++)
            {
                Voxels[Nxy * z + Nx * y + x].SetAttributes(VoxelToWorld(setVecI(x,y,z), VoxelSize, Min), VoxelSize);
            }
        }
    }
    
    
    /* Bounding box */
    for (unsigned int i=0; i<numSpheres; i++)
    {
        computeSphereVertices(spheres[i]);
        cl_float3 tmin = SphereVertex[0], tmax=tmin;
        int vminx,vminy,vminz,vmaxx,vmaxy,vmaxz;
        for (unsigned int j=0; j<8; j++)
        {
            if(SphereVertex[j].x< tmin.x) tmin.x = SphereVertex[j].x;
            if(SphereVertex[j].y< tmin.y) tmin.y = SphereVertex[j].y;
            if(SphereVertex[j].z< tmin.z) tmin.z = SphereVertex[j].z;
            
            if(SphereVertex[j].x > tmax.x) tmax.x = SphereVertex[j].x;
            if(SphereVertex[j].y > tmax.y) tmax.y = SphereVertex[j].y;
            if(SphereVertex[j].z > tmax.z) tmax.z = SphereVertex[j].z;
        }
        
        vminx = WorldToVoxelX(tmin.x, VoxelSize, Min.x); vmaxx = WorldToVoxelX(tmax.x, VoxelSize, Min.x);
        vminy = WorldToVoxelY(tmin.y, VoxelSize, Min.y); vmaxy = WorldToVoxelY(tmax.y, VoxelSize, Min.y);
        vminz = WorldToVoxelZ(tmin.z, VoxelSize, Min.z); vmaxz = WorldToVoxelZ(tmax.z, VoxelSize, Min.z);
        
        
        if(vminx>=Nx) vminx = Nxmin; if (vmaxx >= Nx) vmaxx= Nxmin;
        if(vminy>=Ny) vminy = Nymin; if (vmaxy >= Ny) vmaxy= Nymin;
        if(vminz>=Nz) vminz = Nzmin; if (vmaxz >= Nz) vmaxz= Nzmin;
        
        
        for(int z = vminz; z<=vmaxz; z++)
        {
            for(int y = vminy; y<=vmaxy; y++)
            {
                for(int x = vminx; x<=vmaxx; x++)
                {
                    if(Voxels[Nxy*z+Nx*y+x].IntersectSphere(spheres+i))
                    {
                        Voxels[Nxy*z+Nx*y+x].AddSphere(spheres+i);

                    }
                }
            }
            
        }

        
    }
   
   /* Avoid Looping, Loops are bad */
    
    for (unsigned int i =0; i< numTriangle; i++)
    {
        
        cl_float3 tmin = triag[i].Vertex[0], tmax=tmin;
        int vminx,vminy,vminz,vmaxx,vmaxy,vmaxz;
        
            if(triag[i].Vertex[1].x< tmin.x) tmin.x = triag[i].Vertex[1].x;
            if(triag[i].Vertex[1].y< tmin.y) tmin.y = triag[i].Vertex[1].y;
            if(triag[i].Vertex[1].z< tmin.z) tmin.z = triag[i].Vertex[1].z;
        
            if(triag[i].Vertex[2].x< tmin.x) tmin.x = triag[i].Vertex[2].x;
            if(triag[i].Vertex[2].y< tmin.y) tmin.y = triag[i].Vertex[2].y;
            if(triag[i].Vertex[2].z< tmin.z) tmin.z = triag[i].Vertex[2].z;

        
            if(triag[i].Vertex[1].x > tmax.x) tmax.x = triag[i].Vertex[1].x;
            if(triag[i].Vertex[1].y > tmax.y) tmax.y = triag[i].Vertex[1].y;
            if(triag[i].Vertex[1].z > tmax.z) tmax.z = triag[i].Vertex[1].z;
        
        
            if(triag[i].Vertex[2].x > tmax.x) tmax.x = triag[i].Vertex[2].x;
            if(triag[i].Vertex[2].y > tmax.y) tmax.y = triag[i].Vertex[2].y;
            if(triag[i].Vertex[2].z > tmax.z) tmax.z = triag[i].Vertex[2].z;
        
        
            vminx = WorldToVoxelX(tmin.x, VoxelSize, Min.x); vmaxx = WorldToVoxelX(tmax.x, VoxelSize, Min.x);
            vminy = WorldToVoxelY(tmin.y, VoxelSize, Min.y); vmaxy = WorldToVoxelY(tmax.y, VoxelSize, Min.y);
            vminz = WorldToVoxelZ(tmin.z, VoxelSize, Min.z); vmaxz = WorldToVoxelZ(tmax.z, VoxelSize, Min.z);
    
   
            if(vminx>=Nx) vminx = Nxmin; if (vmaxx >= Nx) vmaxx= Nxmin;
            if(vminy>=Ny) vminy = Nymin; if (vmaxy >= Ny) vmaxy= Nymin;
            if(vminz>=Nz) vminz = Nzmin; if (vmaxz >= Nz) vmaxz= Nzmin;
        
      
        for(int z = vminz; z<=vmaxz; z++)
        {
            for(int y = vminy; y<=vmaxy; y++)
            {
                for(int x = vminx; x<=vmaxx; x++)
                {
                    if(Voxels[Nxy*z+Nx*y+x].IntersectTriangle(triag+i))
                    {
                        Voxels[Nxy*z+Nx*y+x].AddTriangle(triag+i);
                    }
                }
            }
            
        }
        
    }
    

    
    
    
    std::cout<< " Voxels Initialized " << std::endl;

    
  


    
    
}

void UniformGrid::computeSphereVertices(Sphere sph)
{
    SphereVertex[0] = setVec(sph.p.x + sph.rad, sph.p.y + sph.rad, sph.p.z + sph.rad);
    SphereVertex[1] = setVec(sph.p.x + sph.rad, sph.p.y + sph.rad, sph.p.z - sph.rad);
    SphereVertex[2] = setVec(sph.p.x + sph.rad, sph.p.y - sph.rad, sph.p.z - sph.rad);
    SphereVertex[3] = setVec(sph.p.x - sph.rad, sph.p.y - sph.rad, sph.p.z - sph.rad);
    SphereVertex[4] = setVec(sph.p.x + sph.rad, sph.p.y - sph.rad, sph.p.z + sph.rad);
    SphereVertex[5] = setVec(sph.p.x - sph.rad, sph.p.y - sph.rad, sph.p.z + sph.rad);
    SphereVertex[6] = setVec(sph.p.x - sph.rad, sph.p.y + sph.rad, sph.p.z - sph.rad);
    SphereVertex[7] = setVec(sph.p.x - sph.rad, sph.p.y + sph.rad, sph.p.z + sph.rad);

}









float VoxelToWorldX(int x, float VoxelSize, float MinX)
{
    return x * VoxelSize + MinX;
}

float VoxelToWorldY(int y, float VoxelSize, float MinY)
{
    return y * VoxelSize + MinY;
}

float VoxelToWorldZ(int z, float VoxelSize, float MinZ)
{
    return z * VoxelSize + MinZ;
}

cl_float3 VoxelToWorld(const cl_int3 &Voxel, float VoxelSize, cl_float3 Min)
{
    cl_float3 V;
    V.x = (float) Voxel.x * VoxelSize + Min.x;
    V.y = (float) Voxel.y * VoxelSize + Min.y;
    V.z = (float) Voxel.z * VoxelSize + Min.z;
    
    
    return V;
}

int WorldToVoxelX(float x, float VoxelSize, float MinX)
{
    return (int)floor((x - MinX) / VoxelSize);
}

int WorldToVoxelY(float y, float VoxelSize, float MinY)
{
    return (int)floor((y - MinY) / VoxelSize);
}

int WorldToVoxelZ(float z, float VoxelSize, float MinZ)
{
    return (int)floor((z - MinZ) / VoxelSize);
}

cl_int3 WorldToVoxel(const cl_float3 &World, float VoxelSize, cl_float3 Min)
{
    cl_int3 voxel;
    
    voxel.x =  (World.x - Min.x)/VoxelSize;
    voxel.y =  (World.y - Min.y)/VoxelSize;
    voxel.z =  (World.z - Min.z)/VoxelSize;
    
    voxel.x = floor(voxel.x);
    voxel.y = floor(voxel.y);
    voxel.z = floor(voxel.z);
    
    return voxel;
}



float RayTriangleIntersect(Ray r, Triangle *tri)
{
    float det, inv_det, u,v;
    cl_float3 edge1, edge2, tvec, pvec, qvec;
    
    edge1 = tri->Vertex[1];
    edge1.x-= tri->Vertex[0].x;
    edge1.y-= tri->Vertex[0].y;
    edge1.z-= tri->Vertex[0].z;
    
    edge2 = tri->Vertex[2];
    edge2.x -= tri->Vertex[0].x;
    edge2.y -= tri->Vertex[0].y;
    edge2.z -= tri->Vertex[0].z;
    
    pvec = cross(r.d,edge2);
    det  = dot(edge1, pvec);
    
    
    if ( det > 1e-5 && det < 1e-5)
    return 0;
    inv_det = 1.0 /det;
    
    
    /* calculate distance from vert0 to ray origin */
    tvec = r.o;
    tvec.x -= tri->Vertex[0].x;
    tvec.y -= tri->Vertex[0].y;
    tvec.z -= tri->Vertex[0].z;
    
    /* calculate U parameter and test bounds */
    
    u = dot(tvec,pvec) * inv_det;
    if ( u < 0.0 || u>1.0)
    return 0;
    
    /* prepare to test V parameter */
    qvec = cross(tvec, edge1);
    
    /* calculate V parameter and test bounds */
    v = dot(r.d,qvec) * inv_det;
    if (v <0.0 || u + v > 1.0)
    return 0;
    /* calculate t, ray intersects triangle */
    return dot(edge2, qvec) * inv_det;
    
}

float RaySphereIntersect(Ray r, Sphere *sph){ // returns distance, 0 if nohit
    cl_float3 op = sph->p; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    op.x-=r.o.x;
    op.y-=r.o.y;
    op.z-=r.o.z;
    
    float t, eps=1e-4, b=dot(op,r.d), det=b*b-dot(op,op)+sph->rad*sph->rad;
    if (det<0.0) return 0.0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0.0);
}

bool RayBoxIntersect(const Ray r,Box box,float &tmax, float &tmin){
    
    
    
    float  tymin, tymax, tzmin, tzmax;
    cl_float3 inv_d;
    inv_d.x = 1.0f/r.d.x;
    inv_d.y = 1.0f/r.d.y;
    inv_d.z = 1.0f/r.d.z;
    
    int sign[3];
    sign[0] = (inv_d.x<0.0);
    sign[1] = (inv_d.y<0.0);
    sign[2] = (inv_d.z<0.0);
    
    
    
    
    tmin = (box.bounds[sign[0]].x - r.o.x) * inv_d.x;
    tmax = (box.bounds[1-sign[0]].x - r.o.x) * inv_d.x;
    tymin = (box.bounds[sign[1]].y - r.o.y) * inv_d.y;
    tymax = (box.bounds[1-sign[1]].y - r.o.y) * inv_d.y;
    if ( (tmin > tymax) || (tymin > tmax) )
    return false;
    if (tymin > tmin)
    tmin = tymin;
    if (tymax < tmax)
    tmax = tymax;
    tzmin = (box.bounds[sign[2]].z - r.o.z) * inv_d.z;
    tzmax = (box.bounds[1-sign[2]].z - r.o.z) * inv_d.z;
    if ( (tmin > tzmax) || (tzmin > tmax) )
    return false;
    if (tzmin > tmin)
    tmin = tzmin;
    if (tzmax < tmax)
    
    tmax = tzmax;
    return true;
    
    
}



bool Traverse(UniformGrid Grid, const cl_int3 &voxel, const cl_float3 &pos, const Ray &r, float &intersected, Triangle* &triag, Sphere* &sph)
{
    cl_float3 Voxel= setVec(voxel.x,voxel.y,voxel.z),step, out, t, delta = setVec(Grid.VoxelSize,Grid.VoxelSize,Grid.VoxelSize);
    delta.x /= r.d.x;
    delta.y /= r.d.y;
    delta.z /= r.d.z;

    
    
    /* Compute necessary data for Traversal */
    
    /* Determine if the X,Y,Z of the voxel will be incremented or decremented based on the sign of ray direction */
    /* Determine the value of t, t is the computed by the first voxel boundary hit by the ray, also the minimum of t indicates how far can we move in the current voxel */
    if(r.d.x < 0.0)
    {
        step.x = -1.0;
        out.x= Voxel.x <= 0.0 ? Voxel.x - 1.0 : -1.0;
        t.x = (VoxelToWorldX(Voxel.x,Grid.VoxelSize, Grid.Min.x) - pos.x) / r.d.x;
    }
    else
    {
        step.x = 1.0;
        out.x = Voxel.x >= Grid.Nx ? Voxel.x + 1 : Grid.Nx;
        t.x = (VoxelToWorldX(Voxel.x + 1.0,Grid.VoxelSize, Grid.Min.x) - pos.x) / r.d.x;
    }
    
    if(r.d.y < 0.0)
    {
        step.y = -1.0;
        out.y = Voxel.y <= 0.0 ? Voxel.y - 1.0 : -1.0;
        t.y = (VoxelToWorldY(Voxel.y,Grid.VoxelSize, Grid.Min.y) - pos.y) / r.d.y;
    }
    else
    {
        step.y = 1.0;
        out.y = Voxel.y >= Grid.Ny ? Voxel.y + 1 : Grid.Ny;
        t.y = (VoxelToWorldY(Voxel.y + 1.0,Grid.VoxelSize, Grid.Min.y) - pos.y) / r.d.y;
    }
    
    if(r.d.z < 0.0)
    {
        step.z = -1.0;
        out.z = Voxel.z <= 0.0 ? Voxel.z - 1.0 : -1.0;
        t.z = (VoxelToWorldZ(Voxel.z,Grid.VoxelSize, Grid.Min.z) - pos.z) / r.d.z;
    }
    else
    {
        step.z = 1.0;
        out.z = Voxel.z >= Grid.Nz ? Voxel.z + 1 : Grid.Nz;
        t.z = (VoxelToWorldZ(Voxel.z + 1.0,Grid.VoxelSize, Grid.Min.z) - pos.z) / r.d.z;
    }
    
    /* Compute DeltaX,DeltaY,DeltaZ, which indicate how far can we move on the ray in units of t, to the momevement to be equal to the width, height, depth of the voxel */
    
    delta.x *= step.x;
    delta.y *= step.y;
    delta.z *= step.z;

    
    
    
    
    /* Traversal Algorithm by Amanatides and Woo */
    /* Computation of the next voxel */
    
    float interminT,interminS,intersectedS,intersectedT;
    
    triag = NULL;
    sph = NULL;
  

    
    
    
    while(1)
    {
        float min = t.x;
        
        if(t.y < min) min = t.y;
        if(t.z < min) min = t.z;
        
        if(t.x == min)
        {
            Voxel.x += step.x;
            if(Voxel.x == out.x) break;
            t.x += delta.x;
        }
        
        if(t.y == min)
        {
            Voxel.y += step.y;
            if(Voxel.y == out.y) break;
            t.y += delta.y;
        }
        
        if(t.z == min)
        {
            Voxel.z += step.z;
            if(Voxel.z == out.z) break;
            t.z += delta.z;
        }

        

        int x = (int)(Voxel.x), y = (int)Voxel.y, z = (int)Voxel.z;
        if(x >= 0 && x < Grid.Nx && y >= 0 && y < Grid.Ny && z >= 0 && z < Grid.Nz)
        {

            int TrianglesCount = Grid.Voxels[Grid.Nxy * z + Grid.Nx * y + x].TrianglesCount;
            int SpheresCount = Grid.Voxels[Grid.Nxy * z + Grid.Nx * y + x].SpheresCount;
       
            intersectedT=1e20;
            intersectedS=1e20;
            if(TrianglesCount>0) // Hit Happened
            {
                for (int i=0; i<TrianglesCount; i++){
                    /* Mailbox */
                    //if (Voxels[Nxy * z + Nx * y + x].Triangles[i]->id !=r.id)
                    //{
                      //  Voxels[Nxy * z + Nx * y + x].Triangles[i]->id = r.id;
                   
                        if( (interminT = RayTriangleIntersect(r,  Grid.Voxels[Grid.Nxy * z + Grid.Nx * y + x].Triangles[i])) && interminT <  intersectedT)
                        {intersectedT=interminT;  triag = Grid.Voxels[Grid.Nxy * z + Grid.Nx * y + x].Triangles[i];}
                        
                    //}
                }
            }
            
            if(SpheresCount>0)
            {
                for (int i=0; i<SpheresCount; i++){
                    /* Mailbox */
                     //if (Voxels[Nxy * z + Nx * y + x].Spheres[i]->id !=r.id)
                    //{
                      //  Voxels[Nxy * z + Nx * y + x].Spheres[i]->id = r.id;
                       if((interminS = RaySphereIntersect(r,Grid.Voxels[Grid.Nxy * z + Grid.Nx * y + x].Spheres[i])) &&interminS<intersectedS)
                        {intersectedS = interminS;  sph = Grid.Voxels[Grid.Nxy * z +Grid.Nx * y + x].Spheres[i];}
                    
                    //}
                }
            }
            
            if (intersectedS > intersectedT){intersected = intersectedT; sph=NULL; return true;}
            if (intersectedT > intersectedS){intersected = intersectedS; triag=NULL; return true;}

           
        }
    }
    

    return false;
}


