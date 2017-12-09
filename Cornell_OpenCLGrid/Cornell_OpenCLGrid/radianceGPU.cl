//
//  radianceGPU.cl
//
//
//  Created by Stylianos Piperakis on 5/23/15.
//  Copyright (c) 2015 Stylianos Piperakis. All rights reserved.
//

/* ------------------------------------------------------------------------------------------- */
typedef struct Ray_ {
 
 float3 o, d, inv_d;
 int3 sign;
 
 
 }Ray;
/* ------------------------------------------------------------------------------------------- */

typedef struct Box_{
    
    float3 bounds[2];
    
}Box;

 typedef struct Sphere_{
 
 float rad;       // radius
 float3 p, e, c;      // position, emission, color
 unsigned int refl;      // reflection type (DIFFuse, SPECular, REFRactive)
 
 }Sphere;
/* ------------------------------------------------------------------------------------------- */

 typedef struct Triangle_
 {
 float3 Vertex[3], Normal, Color, e;
 unsigned int  refl;
 }Triangle;
 
/* ------------------------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------------------------- */

 
/* Voxel Struct */
typedef struct Voxel_
    {
        float3 Max;
        float3 Min;
        float3 Max_;
        float3 Min_;  //For comparison reasons
        int MaxTriangleLimit;
        int MaxSphereLimit;
        
        Triangle** Triangles;
        Sphere** Spheres;
        int TrianglesCount;
        int SpheresCount;
        float Size;
    }Voxel;

typedef struct UniformGrid_
    {
        
        float3 d;
        int k;
        float Volume;
        float3 SphereVertex[8]; // For bounding box of sphere
        
        float3 Min,Max;
        Voxel *Voxels;
        
        Triangle *triag;
        Sphere *spheres;
        
        float VoxelSize;
        int numPrim, numTriangle, numSpheres;
        int Nx,Ny,Nz;
        int Nxy,Nxyz, Nxmin, Nymin, Nzmin;
        
            
        
        
        
        
    }UniformGrid;





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

float3 VoxelToWorld(int3 Voxel, float VoxelSize, float3 Min)
{
    float3 V;
    V.x = (float) Voxel.x * VoxelSize + Min.x;
    V.y = (float) Voxel.y * VoxelSize + Min.y;
    V.z = (float) Voxel.z * VoxelSize + Min.z;
    
    
    return V;
}

int WorldToVoxelX(float x, float VoxelSize, float MinX)
{
    return (int)floor( (float) ((x - MinX) / VoxelSize));
}

int WorldToVoxelY(float y, float VoxelSize, float MinY)
{
    return (int)floor( (float)((y - MinY) / VoxelSize));
}

int WorldToVoxelZ(float z, float VoxelSize, float MinZ)
{
    return (int)floor( (float)((z - MinZ) / VoxelSize));
}

int3 WorldToVoxel(float3 World, float VoxelSize, float3 Min)
{
    int3 voxel;
    
    voxel.x =  (World.x - Min.x)/VoxelSize;
    voxel.y =  (World.y - Min.y)/VoxelSize;
    voxel.z =  (World.z - Min.z)/VoxelSize;
    
    voxel.x = floor((float) (voxel.x));
    voxel.y = floor((float) (voxel.y));
    voxel.z = floor((float) (voxel.z));
    
    return voxel;
}



float RayTriangleIntersect(Ray r, Triangle *tri)
{
    float det, inv_det, u,v;
    float3 edge1, edge2, tvec, pvec, qvec;
    
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
    float3 op = sph->p; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    op.x-=r.o.x;
    op.y-=r.o.y;
    op.z-=r.o.z;
    
    float t, eps=1e-4, b=dot(op,r.d), det=b*b-dot(op,op)+sph->rad*sph->rad;
    if (det<0.0) return 0.0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0.0);
}

bool RayBoxIntersect( Ray r,Box box,float** tmaxf, float** tminf){
    
    
    
    float  tymin, tymax, tzmin, tzmax, tmax, tmin;
    float3 inv_d;
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
    
    *tmaxf = &tmax;
    *tminf = &tmin;
    return true;
    
    
}





/* ------------------------------------------------------------------------------------------- */
                            /*  Random function */
/* ------------------------------------------------------------------------------------------- */
 float GetRandom(unsigned int *seed0, unsigned int *seed1) {
 *seed0 = 36969 * ((*seed0) & 65535) + ((*seed0) >> 16);
 *seed1 = 18000 * ((*seed1) & 65535) + ((*seed1) >> 16);
 
 unsigned int ires = ((*seed0) << 16) + (*seed1);
 
 
 union {
 float f;
 unsigned int ui;
 } res;
 res.ui = (ires & 0x007fffff) | 0x40000000;
 
 return (res.f - 2.f) / 2.f;
 }
 

/* ------------------------------------------------------------------------------------------- */
                        /* Main Radiance function */
/* ------------------------------------------------------------------------------------------- */
float3 radiance(float3 currentPos,Ray r, __constant UniformGrid* Grid, unsigned int *seed0, unsigned int *seed1){
 
    unsigned int depth = 0;
    float3 color={0.0,0.0,0.0};
    float3 cf = { 1.0,1.0,1.0};
    
    Triangle objT;
    Sphere objS;
    
    float3 step, out, t, delta, pos;
    int3 Voxel;
    int idx;
    /* For recursion */

    
 while(1){
     
    /* Check ray-primitive intersection */
    
     pos = currentPos;
     Voxel = WorldToVoxel(currentPos, Grid->VoxelSize, Grid->Min);
     idx=-1;

 /* -------------* Traversal Initiated *------------------*/
    
     delta.x = Grid->VoxelSize;
     delta.y = Grid->VoxelSize;
     delta.z = Grid->VoxelSize;
     delta.x /= r.d.x;
     delta.y /= r.d.y;
     delta.z /= r.d.z;
     
     
     
     if(r.d.x < 0.0)
     {
         step.x = -1.0;
         out.x= Voxel.x <= 0.0 ? Voxel.x - 1.0 : -1.0;
         t.x = (VoxelToWorldX(Voxel.x,Grid->VoxelSize, Grid->Min.x) - pos.x) / r.d.x;
     }
     else
     {
         step.x = 1.0;
         out.x = Voxel.x >= Grid->Nx ? Voxel.x + 1 : Grid->Nx;
         t.x = (VoxelToWorldX(Voxel.x + 1.0,Grid->VoxelSize, Grid->Min.x) - pos.x) / r.d.x;
     }
     
     if(r.d.y < 0.0)
     {
         step.y = -1.0;
         out.y = Voxel.y <= 0.0 ? Voxel.y - 1.0 : -1.0;
         t.y = (VoxelToWorldY(Voxel.y,Grid->VoxelSize, Grid->Min.y) - pos.y) / r.d.y;
     }
     else
     {
         step.y = 1.0;
         out.y = Voxel.y >= Grid->Ny ? Voxel.y + 1 : Grid->Ny;
         t.y = (VoxelToWorldY(Voxel.y + 1.0,Grid->VoxelSize, Grid->Min.y) - pos.y) / r.d.y;
     }
     
     if(r.d.z < 0.0)
     {
         step.z = -1.0;
         out.z = Voxel.z <= 0.0 ? Voxel.z - 1.0 : -1.0;
         t.z = (VoxelToWorldZ(Voxel.z,Grid->VoxelSize, Grid->Min.z) - pos.z) / r.d.z;
     }
     else
     {
         step.z = 1.0;
         out.z = Voxel.z >= Grid->Nz ? Voxel.z + 1 : Grid->Nz;
         t.z = (VoxelToWorldZ(Voxel.z + 1.0,Grid->VoxelSize, Grid->Min.z) - pos.z) / r.d.z;
     }
     
     
     delta.x *= step.x;
     delta.y *= step.y;
     delta.z *= step.z;
     
     
     
     
     
     float interminT,interminS,intersectedS,intersectedT,intersected;
     
     /* ------------* Traversal  Algorithm *------------------*/
     
     
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
         if(x >= 0 && x < Grid->Nx && y >= 0 && y < Grid->Ny && z >= 0 && z < Grid->Nz)
         {
             
             int TrianglesCount = Grid->Voxels[Grid->Nxy * z + Grid->Nx * y + x].TrianglesCount;
             int SpheresCount = Grid->Voxels[Grid->Nxy * z + Grid->Nx * y + x].SpheresCount;
             intersectedT=1e20;
             intersectedS=1e20;
             if(TrianglesCount>0) // Hit Happened
             {
                 for (int i=0; i<TrianglesCount; i++){
                     
                     if( (interminT = RayTriangleIntersect(r,  Grid->Voxels[Grid->Nxy * z + Grid->Nx * y + x].Triangles[i])) && interminT <  intersectedT)
                     {
                         intersectedT=interminT;

                         objT = *(Grid->Voxels[Grid->Nxy * z + Grid->Nx * y + x].Triangles[i]);
                     
                     }

                 }
             }
             
             if(SpheresCount>0)
             {
                 for (int i=0; i<SpheresCount; i++){
                     if((interminS = RaySphereIntersect(r,Grid->Voxels[Grid->Nxy * z + Grid->Nx * y + x].Spheres[i])) &&interminS<intersectedS)
                     {
                         intersectedS = interminS;


                         objS = *(Grid->Voxels[Grid->Nxy * z +Grid->Nx * y + x].Spheres[i]);
                     }
                     
                 }
             }
             if (intersectedS > intersectedT){intersected = intersectedT;  idx=1; }
             if (intersectedT > intersectedS){intersected = intersectedS;  idx=0; }

             
             if(idx != -1) break;
         }
     }
     
    /* ------------------------------------------------------------------------------------------*/
     
    if (idx==-1)  return color;
     
     /* ------------------------------------------------------------------------------------------*/

     
     if(idx==0)
     {

         Sphere obj=objS;
         
         float3 x;
         
         float3 rdir;
         rdir.x = r.d.x * intersected;
         rdir.y = r.d.y * intersected;
         rdir.z = r.d.z * intersected;

         x= r.o + rdir;
         currentPos=x;

         float3 n = (x - obj.p);
         n=normalize(n);
         float3 nl = dot(r.d,n);
         if (nl.x>=0 && nl.y>=0 && nl.z>0) n=-n;
         float3 f=obj.c;
         
         
         float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
         color = color + cf*obj.e;
         
         
         if (++depth>5)
             if (GetRandom(seed0,seed1)<p) f=f*(1/p);
             else {return color;}
         cf=cf*f;
       
        if(obj.refl == 0)
        {
         float r1 = 2.0f * 3.14159265358979323846f * GetRandom(seed0, seed1);
         float r2 = GetRandom(seed0, seed1);
         float r2s = sqrt(r2);
         
         float3 w; {(w).x = (nl).x; (w).y = (nl).y; (w).z = (nl).z; };
         
         float3 u, a;
         if (fabs(w.x) > 0.1f) {
             { (a).x = 0.0f; (a).y = 1.0f; (a).z = 0.0f; };
             u=a;
         } else {
             { (a).x = 1.0f; (a).y = 0.0f; (a).z = 0.0f; };
             a = cross(a,w);
             u=a;
         }
         u=normalize(u);
         float3 v;
         v= cross(w,u);
         
         float3 d;
         float temp1 = cos(r1)*r2s;
         float temp2 = sin(r1)*r2s;
         float temp3 = sqrt(1-r2);
         
         d.x = u.x * temp1 + v.x * temp2 + w.x * temp3;
         d.y = u.y * temp1 + v.y * temp2 + w.y * temp3;
         d.z = u.z * temp1 + v.z * temp2 + w.z * temp3;
         
         d=normalize(d);
         
         r.o = x;
         r.d = d;
         
         
         continue;
         
 
     }
        
     else if(obj.refl == 1)
     {
         r.o = x;
         float3 d;
         d= r.d;
         float ref=2.0 * dot(r.d,n);
         float3 tempV;
         tempV.x = n.x * ref;
         tempV.y = n.y * ref;
         tempV.z = n.z * ref;
         d = r.d - tempV;
         r.d = d;
         
         continue;
     }
         // Ideal dielectric REFRACTION
         Ray reflRay;
         reflRay.o = x;
         float3 d = r.d;
         float ref=2.0 * dot(r.d,n);
         float3 tempV;
         tempV.x = n.x * ref;
         tempV.y = n.y * ref;
         tempV.z = n.z * ref;
         d = r.d - tempV;
         reflRay.d = d;
         // Ray from outside going in?
         bool into = dot(n,nl) > 0;
         float nc=1.0, nt = 1.5, nnt=into?nc/nt:nt/nc, ddn =dot(r.d,nl), cos2t;
         // Total internal reflection
         if ((cos2t=1.0-nnt*nnt*(1.0-ddn*ddn))<0)
         {    // Total internal reflection
             //return obj.e + f.mult(radiance(reflRay,depth,Xi));
             r=reflRay;
             continue;
         }
         
         float3  tdir;
         tdir.x = r.d.x * nnt - n.x* ((into?1.0:-1.0)* ddn*nnt+sqrt(cos2t));
         tdir.y = r.d.y * nnt - n.y* ((into?1.0:-1.0)* ddn*nnt+sqrt(cos2t));
         tdir.z = r.d.z * nnt - n.z* ((into?1.0:-1.0)* ddn*nnt+sqrt(cos2t));
         tdir=normalize(tdir);
         float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1.0 - (into?-ddn: dot(n,tdir));
         float Re=R0+(1.0 - R0) * c * c * c * c *c, Tr=1-Re, P=0.25+0.5 *Re, RP=Re/P, TP=Tr/(1.0-P);
         
         
         if (GetRandom(seed0, seed1)<P)
         {
             cf=cf*RP;
             r = reflRay;
         }
         else
         {
             
             cf=cf*TP;
             r.o = x;
             r.d = tdir;
         }
         continue;
         
     }
     
 
     
     
     
     else
     {
       
         Triangle obj=objT;
         float3 x;
         
         float3 rdir;
         rdir.x = r.d.x * intersected;
         rdir.y = r.d.y * intersected;
         rdir.z = r.d.z * intersected;
         
         x= r.o + rdir;
         currentPos=x;
         float3 n = cross(obj.Vertex[1] - obj.Vertex[0],obj.Vertex[2]-obj.Vertex[0]);
         n=normalize(n);
         float3 nl = dot(r.d,n);
         if (nl.x>=0 && nl.y>=0 && nl.z>0) n=-n;
         float3 f=obj.Color;
         
         
         float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
         color = color + cf*obj.e;
         
         
         if (++depth>5)
             if (GetRandom(seed0,seed1)<p) f=f*(1/p);
             else return color;
         cf=cf*f;
         
         if(obj.refl == 0)
         {
             float r1 = 2.0f * 3.14159265358979323846f * GetRandom(seed0, seed1);
             float r2 = GetRandom(seed0, seed1);
             float r2s = sqrt(r2);
             
             float3 w; {(w).x = (nl).x; (w).y = (nl).y; (w).z = (nl).z; };
             
             float3 u, a;
             if (fabs(w.x) > 0.1f) {
                 { (a).x = 0.0f; (a).y = 1.0f; (a).z = 0.0f; };
                 u=a;
             } else {
                 { (a).x = 1.0f; (a).y = 0.0f; (a).z = 0.0f; };
                 a = cross(a,w);
                 u=a;
             }
             u=normalize(u);
             float3 v;
             v= cross(w,u);
             
             float3 d;
             float temp1 = cos(r1)*r2s;
             float temp2 = sin(r1)*r2s;
             float temp3 = sqrt(1-r2);
             
             d.x = u.x * temp1 + v.x * temp2 + w.x * temp3;
             d.y = u.y * temp1 + v.y * temp2 + w.y * temp3;
             d.z = u.z * temp1 + v.z * temp2 + w.z * temp3;
             
             d=normalize(d);
             
             r.o = x;
             r.d = d;
             
             
             continue;
             
             
         }
         
         else if(obj.refl == 1)
         {
             r.o = x;
             float3 d;
             d= r.d;
             float ref=2.0 * dot(r.d,n);
             float3 tempV;
             tempV.x = n.x * ref;
             tempV.y = n.y * ref;
             tempV.z = n.z * ref;
             d = r.d - tempV;
             r.d = d;
             
             continue;
         }
         // Ideal dielectric REFRACTION
         Ray reflRay;
         reflRay.o = x;
         float3 d = r.d;
         float ref=2.0 * dot(r.d,n);
         float3 tempV;
         tempV.x = n.x * ref;
         tempV.y = n.y * ref;
         tempV.z = n.z * ref;
         d = r.d - tempV;
         reflRay.d = d;
         // Ray from outside going in?
         bool into = dot(n,nl) > 0;
         float nc=1.0, nt = 1.5, nnt=into?nc/nt:nt/nc, ddn =dot(r.d,nl), cos2t;
         // Total internal reflection
         if ((cos2t=1.0-nnt*nnt*(1.0-ddn*ddn))<0)
         {    // Total internal reflection
             //return obj.e + f.mult(radiance(reflRay,depth,Xi));
             r=reflRay;
             continue;
         }
         
         float3  tdir;
         tdir.x = r.d.x * nnt - n.x* ((into?1.0:-1.0)* ddn*nnt+sqrt(cos2t));
         tdir.y = r.d.y * nnt - n.y* ((into?1.0:-1.0)* ddn*nnt+sqrt(cos2t));
         tdir.z = r.d.z * nnt - n.z* ((into?1.0:-1.0)* ddn*nnt+sqrt(cos2t));
         tdir=normalize(tdir);
         float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1.0 - (into?-ddn: dot(n,tdir));
         float Re=R0+(1.0 - R0) * c * c * c * c *c, Tr=1-Re, P=0.25+0.5 *Re, RP=Re/P, TP=Tr/(1.0-P);
         
         
         if (GetRandom(seed0, seed1)<P)
         {
             cf=cf*RP;
             r = reflRay;
         }
         else
         {
             
             cf=cf*TP;
             r.o = x;
             r.d = tdir;
         }
         continue;
         
        
     }
 }
}

/* ------------------------------------------------------------------------------------------- */
                                // Generate the Ray //
/* ------------------------------------------------------------------------------------------- */

 Ray genRay(const int width,const int height,const int x,const int y,int sx,int sy, unsigned int *seed0, unsigned int *seed1)
 {
     Ray r;
 
  
     
     
     
      float r1 =2.0 * GetRandom(seed0, seed1), dx=r1<1? sqrt(r1)-1: 1-sqrt(2-r1);
      float r2 =2.0 * GetRandom(seed0, seed1), dy=r2<1? sqrt(r2)-1: 1-sqrt(2-r2);
     
    
     float3 camOrigin={50,52,295.6};
     float3 camDir = {0,-0.042612,-1};
     camDir=normalize(camDir);
     float3 cx;
     cx.x=width*0.5135/height;
     cx.y=0.0;
     cx.z=0.0;
     
     float3 cy;
     cy=cross(cx,camDir);
     cy=normalize(cy);
     cy.x = cy.x *0.5135;
     cy.y = cy.y *0.5135;
     cy.z = cy.z *0.5135;

     
     
     float3 temp1;
     temp1.x = ((sx+0.5 + dx)/2.0 + x)/width - 0.5;
     temp1.y = ((sx+0.5 + dx)/2.0 + x)/width - 0.5;
     temp1.z = ((sx+0.5 + dx)/2.0 + x)/width - 0.5;
     
    
     float3 temp2;
     temp2.x =  ( (sy+0.5 + dy)/2.0 + y)/height - 0.5;
     temp2.y =  ( (sy+0.5 + dy)/2.0 + y)/height - 0.5;
     temp2.z =  ( (sy+0.5 + dy)/2.0 + y)/height - 0.5;

     float3 d = cx*temp1 + cy*temp2 + camDir;
     
     d.x =d.x *140.0;
     d.y =d.y *140.0;
     d.z =d.z *140.0;

     r.o=camOrigin+d;
     r.d=normalize(d);

     
     
     
 return r;
 }
 

__kernel void radianceGPU(__constant UniformGrid* Grid,__global float3* colors, __constant Sphere* spheres, __constant Triangle* triag,  __global unsigned int *randomInput, const int width, const int height, const unsigned int sphereCount, const unsigned int triangleCount)
{
    
    /* ------------ */
    const int samps =2;
    /* ------------ */

    
        
    const int gid = get_global_id(0);
    const int gid2 = 2 * gid;
    const int x = gid % width;
    const int y = gid / width;
    
    
    if (y >= height)
        return;
    
    unsigned int seed0 = randomInput[gid2];
    unsigned int seed1 = randomInput[gid2+1];
    unsigned int i=(height-y-1)*width+x;
    float3 currentPos;
    /* Init the color */
    colors[i].x = 0.0;
    colors[i].y = 0.0;
    colors[i].z = 0.0;

    
    float* tmin;
    float* tmax;
    Box BoundingBox;
    BoundingBox.bounds[0]=Grid->Min;
    BoundingBox.bounds[1]=Grid->Max;
    for (int sy=0; sy<2; sy++){    // 2x2 subpixel rows
        float3 r = {0.0,0.0,0.0};
        for (int sx=0; sx<2; sx++){        // 2x2 subpixel cols
            for (int s=0; s<samps; s++){
                
                
                
                Ray ray=genRay(width,height,x,y,sx,sy, &seed0, &seed1);
                if(RayBoxIntersect(ray,BoundingBox,&tmax,&tmin)){
                    if((*tmin)<0) (*tmin)=0;
                
                    currentPos =ray.o;
                    currentPos.x += (*tmin) * ray.d.x;
                    currentPos.y += (*tmin) * ray.d.y;
                    currentPos.z += (*tmin) * ray.d.z;
                    
                r = r +  radiance(currentPos,ray,Grid, &seed0, &seed1)* (float3)(1.0/samps,1.0/samps,1.0/samps);
                }
            } // Camera rays are pushed ^^^^^ forward to start in interior
            colors[i].x =colors[i].x + clamp((float)r.x,(float)0.0,(float)1.0)*0.25;
            colors[i].y =colors[i].y + clamp((float)r.y,(float)0.0,(float)1.0)*0.25;
            colors[i].z =colors[i].z + clamp((float)r.z,(float)0.0,(float)1.0)*0.25;
        }
    }
     randomInput[gid2] = seed0;
     randomInput[gid2+1] = seed1;

}
