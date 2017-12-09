//
//  scene.h
//
//
//  Created by Stylianos Piperakis on 5/23/15.
//  Copyright (c) 2015 Stylianos Piperakis. All rights reserved.
//

#ifndef __mySmallPT__scene__
#define __mySmallPT__scene__



#include "primitives.h"





Sphere  spheres[]={
    Sphere(16.5,setVec(27,16.5,47),setVec(0.0) ,setVec(0.999,0.999,0.999), 1),//Mirr
    Sphere(12.5,setVec(73,10.5,120),setVec(0.0),setVec(0.999,0.999,0.999), 2)//Glas
};




/* Cube */
/* Back Face */
cl_float3 Vertex1[3]={setVec(17.0,0.0,110.0),setVec(37.0,0.0 ,110.0),  setVec(37,20.5,110)};
cl_float3 Vertex2[3]={setVec(17,0.0,110.0),setVec(17,20.5 ,110.0),setVec(37,20.5,110.0)};



/* Front Face */
cl_float3 Vertex3[3]={setVec(37,20.5,130.0),setVec(37,0.0 ,130.0),setVec(17,0.0,130.0)};
cl_float3 Vertex4[3]={setVec(37,20.5,130.0),setVec(17,0.0 ,130.0),setVec(17,20.5,130.0)};

/*Left  Face*/
cl_float3 Vertex5[3]={setVec(17,20.5,130.0),setVec(17,0.0 ,130.0),setVec(17,0.0,110.0)};
cl_float3 Vertex6[3]={setVec(17,20.5,130.0),setVec(17,0.0 ,110.0),setVec(17,20.5,110.0)};

/* Right Face */
cl_float3 Vertex7[3]={setVec(37,20.5,130.0),setVec(37,0.0 ,130.0),setVec(37,0.0,110.0)};
cl_float3 Vertex8[3]={setVec(37,20.5,130.0),setVec(37,0.0 ,110.0),setVec(37,20.5,110.0)};

/* Botom Face */
cl_float3 Vertex9[3]={setVec(17,0.0,130.0),setVec(37,0.0 ,130.0),setVec(17,0.0,110.0)};
cl_float3 Vertex10[3]={setVec(17,0.0,110.0),setVec(37,0.0 ,130.0),setVec(37,0.0,110.0)};

/* Top Face */
cl_float3 Vertex29[3]={setVec(17,20.5,130.0),setVec(37,20.5 ,130.0),setVec(17,20.5,110.0)};
cl_float3 Vertex30[3]={setVec(17,20.5,110.0),setVec(37,20.5 ,130.0),setVec(37,20.5,110.0)};


/* Pyramid */
/*Bottom Face */
cl_float3 Vertex11[3]={ setVec(85,0.0,10.0),setVec(95.0,0.0 ,25.0),setVec(85.0, 0.0,40.0)};
cl_float3 Vertex12[3]={setVec(85.0,0.0 ,10.0), setVec(75.0,0.0,25.0),setVec(85.0,0.0 ,40.0)};
/* Apex */
cl_float3 Vertex13[3]={setVec(85.0,46.5,25.0),setVec(85.0,0.0 ,40.0),setVec(75.0,0.0,25.0)};
cl_float3 Vertex14[3]={setVec(85.0,46.5,25.0),setVec(85.0,0.0,10.0), setVec(75.0,0.0,25.0)};

cl_float3 Vertex15[3]={ setVec(85,46.5,25.0), setVec(95,0.0,25.0),setVec(85,0.0,40.0)};
cl_float3 Vertex16[3]={ setVec(85,46.5,25.0),setVec(95,0.0,25.0), setVec(85,0.0,10.0)};


/* Back*/
cl_float3 Vertex17[3]={setVec(1.0,0.0,1.0),setVec(99.0,0.0 ,1.0),setVec(99.0, 81.6,1.0)};
cl_float3 Vertex18[3]={setVec(1.0,81.6 ,1.0),setVec(99.0,81.6,1.0),setVec(1.0,0.0,1.0)};

/* Front */
cl_float3 Vertex19[3]={setVec(1.0,0.0,170.0),setVec(99.0,0.0 ,170.0),setVec(99.0, 81.6,170.0)};
cl_float3 Vertex20[3]={setVec(1.0,81.6 ,170.0),setVec(99.0,81.6,170.0),setVec(1.0,0.0,170.0)};

/* Bottom */
cl_float3 Vertex21[3]={setVec(1.0,0.0,1.0),setVec(99.0,0.0 ,1.0),setVec(1.0, 0.0,170.0)};
cl_float3 Vertex22[3]={setVec(1.0,0.0 ,170.0),setVec(99.0,0.0,170.0),setVec(99.0,0.0,1.0)};

/* Top */
cl_float3 Vertex23[3]={setVec(1.0,81.6,1.0),setVec(99.0,81.6 ,1.0),setVec(1.0, 81.6,170.0)};
cl_float3 Vertex24[3]={setVec(1.0,81.6 ,170.0),setVec(99.0,81.6,170.0),setVec(99.0,81.6,1.0)};

/* Left */
cl_float3 Vertex25[3]={setVec(1.0,81.6,170.0),setVec(1.0,81.6 ,1.0),setVec(1.0,0.0,1.0)};
cl_float3 Vertex26[3]={setVec(1.0,81.6 ,170.0),setVec(1.0,0.0,170.0),setVec(1.0,0.0,1.0)};

/* Right */
cl_float3 Vertex27[3]={setVec(99.0,81.6,170.0),setVec(99.0,81.6 ,1.0),setVec(99.0,0.0,1.0)};
cl_float3 Vertex28[3]={setVec(99.0,81.6 ,170.0),setVec(99.0,0.0,170.0),setVec(99.0,0.0,1.0)};

/* Light */
cl_float3 Vertex31[3]={setVec(66 , 81.6-0.5,56),setVec(66, 81.6-0.5,112),setVec(33,81.6-0.5,112)}; //Top Face
cl_float3 Vertex32[3]={setVec(33,81.6-0.5,112), setVec(33,81.6-0.5,56), setVec(66,81.6-0.5,56)};
cl_float3 Vertex33[3]={setVec(66 , 75.6-0.5,56),setVec(66, 75.6-0.5,112),setVec(33,75.6-0.5,112)}; // Bottom Face
cl_float3 Vertex34[3]={setVec(33,75.6-0.5,112), setVec(33,75.6-0.5,56), setVec(66,75.6-0.5,56)};

cl_float3 Vertex35[3]={setVec(66 , 81.6-0.5,56),setVec(66, 75.6-0.5,56),setVec(33,75.6-0.5,56)}; // Front Face
cl_float3 Vertex36[3]={setVec(33,75.6-0.5,56),setVec(33, 81.6-0.5,56),setVec(66,81.6-0.5,56)};

cl_float3 Vertex37[3]={setVec(66 , 81.6-0.5,112),setVec(66, 75.6-0.5,112),setVec(33,75.6-0.5,112)}; // Back Face
cl_float3 Vertex38[3]={setVec(33,75.6-0.5,112),setVec(33, 81.6-0.5,112),setVec(66,81.6-0.5,112)};



cl_float3 Vertex39[3]={setVec(66 , 81.6-0.5,56),setVec(66, 75.6-0.5,56),setVec(66,75.6-0.5,112)}; // Left Face
cl_float3 Vertex40[3]={setVec(66,75.6-0.5,112),setVec(66, 81.6-0.5,112),setVec(66,81.6-0.5,56)};

cl_float3 Vertex41[3]={setVec(33 , 81.6-0.5,56),setVec(33, 75.6-0.5,56),setVec(33,75.6-0.5,112)}; // Right Face
cl_float3 Vertex42[3]={setVec(33,75.6-0.5,112),setVec(33, 81.6-0.5,112),setVec(33,81.6-0.5,56)};




Triangle triag[] = {
    
    /* CUBE */
    //setVec(0.4,0.0,0.8)
    Triangle(Vertex1,  setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex2 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex3 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex4 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex5 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex6 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex7 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex8 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex9 , setVec(0), setVec(0.4,0.0,0.8),0),
    Triangle(Vertex10 , setVec(0),setVec(0.4,0.0,0.8),0),
    
    /* Pyramid */
    
    Triangle(Vertex11 , setVec(0), setVec(0.0,0.6,0.0),0),
    Triangle(Vertex12 , setVec(0), setVec(0.0,0.6,0.0),0),
    Triangle(Vertex13 , setVec(0), setVec(0.0,0.6,0.0),0),
    Triangle(Vertex14 , setVec(0), setVec(0.0,0.6,0.0),0),
    Triangle(Vertex15 , setVec(0), setVec(0.0,0.6,0.0),0),
    Triangle(Vertex16 , setVec(0), setVec(0.0,0.6,0.0),0),
    
    
    /* Back */
    Triangle(Vertex17, setVec(0),setVec(.75,.75,.75),0),
    Triangle(Vertex18, setVec(0),setVec(.75,.75,.75),0),
    
    /* Front */
    //Triangle(Vertex19, setVec(0),setVec(0),0),
    //Triangle(Vertex20, setVec(0),setVec(0),0),
    
    /* Bottom */
    Triangle(Vertex21, setVec(0),setVec(.75,.75,.75),0),
    Triangle(Vertex22, setVec(0),setVec(.75,.75,.75),0),
    
    /* Top */
    Triangle(Vertex23, setVec(0),setVec(.75,.75,.75),0),
    Triangle(Vertex24, setVec(0),setVec(.75,.75,.75),0),
    /* Left */
    Triangle(Vertex25, setVec(0),setVec(.75,.25,.25),0),
    Triangle(Vertex26, setVec(0),setVec(.75,.25,.25),0),
    
    /*Right*/
    Triangle(Vertex27, setVec(0),setVec(.25,.25,.75),0),
    Triangle(Vertex28, setVec(0),setVec(.25,.25,.75),0),
    
    /* Cube */
    Triangle(Vertex29, setVec(0),setVec(0.4,0.0,0.8),0),
    Triangle(Vertex30, setVec(0),setVec(0.4,0.0,0.8),0),
    
    /* Light */
    Triangle(Vertex31, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex32, setVec(24.0,12.0,24.0),setVec(0),0),
    Triangle(Vertex33, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex34, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex35, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex36, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex37, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex38, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex39, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex40, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex41, setVec(24.0,24.0,24.0),setVec(0),0),
    Triangle(Vertex42, setVec(24.0,24.0,24.0),setVec(0),0)
    
    
};




#endif /* defined(__mySmallPT__scene__) */
