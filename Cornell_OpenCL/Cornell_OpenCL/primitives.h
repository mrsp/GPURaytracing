//
//  primitives.h
//  
//
//  Created by Stylianos Piperakis on 5/23/15.
//  Copyright (c) 2015 Stylianos Piperakis. All rights reserved.
//

#ifndef __mySmallPT__prim__
#define __mySmallPT__prim__

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif


float clamp(float x){ return x<0 ? 0 : x>1 ? 1 : x; };
int toInt(float x){ return int(pow(clamp(x),1/2.2)*255+.5); };

typedef struct Ray {
    
    cl_float3 o, d, inv_d;
    cl_int3 sign;
    
    Ray(cl_float3 o_, cl_float3 d_)
    {
        o=o_;
        d=d_;
        inv_d.x=1.0f/d.x;
        inv_d.y=1.0f/d.y;
        inv_d.z=1.0f/d.z;
        sign.x = (inv_d.x<0.0);
        sign.y = (inv_d.y<0.0);
        sign.z = (inv_d.z<0.0);
    }

    
}Ray;



typedef struct Sphere {
    
    cl_float rad;       // radius
    cl_float3 p, e, c;      // position, emission, color
    unsigned int refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(float rad_, cl_float3 p_, cl_float3 e_, cl_float3 c_, unsigned int refl_)
    {
        rad=rad_;
        p=p_;
        e=e_;
        c=c_;
        refl=refl_;
    }
    
}Sphere;





typedef struct Triangle
{
    cl_float3 Vertex[3], Normal, Color, e;
    unsigned int  refl;
    Triangle(cl_float3 Vertex_[3], cl_float3 e_, cl_float3 c_, unsigned int refl_)
    {
        
        Vertex[0]=Vertex_[0];
        Vertex[1]=Vertex_[1];
        Vertex[2]=Vertex_[2];
        e=e_;
        Color=c_;
        refl=refl_;
    }
    
}Triangle;

cl_float3 setVec(float x=0.0,float y=0.0,float z=0.0)
{
    cl_float3 V;
    V.x=x;
    V.y=y;
    V.z=z;
    
    return V;
}




#endif
