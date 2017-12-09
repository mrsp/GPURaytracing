//
//  main.cpp
//
//
//  Created by Stylianos Piperakis on 5/23/15.
//  Copyright (c) 2015 Stylianos Piperakis. All rights reserved.
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include "scene.h"
#include <chrono>




unsigned int *pixels;
int width = 600;
int height = 600;

cl_float3 *colors;
unsigned int *seeds;
size_t pixelCount = (size_t) width * (size_t) height;

unsigned int sphereCount = sizeof(spheres)/sizeof(Sphere);
unsigned int triangleCount = sizeof(triag)/sizeof(Triangle);

bool USE_CPU = true;
int main()
{
    
  
    
    /* OpenCL varianbles */
    
    cl_uint numPlatforms;

    cl_platform_id platform;
    
    cl_device_id device,devices[32];
    
    cl_context context;
    
    cl_command_queue commandQueue;
    
    cl_uint  deviceCount,addr_data;
    
    cl_int  err;
    
    char name_data[48], ext_data[4096];
    
    cl_program program;
    
    FILE *program_handle;
    
    char *program_buffer;
    
    
    size_t program_size;
    
    
    cl_kernel kernel;
    
    
    /* OpenCL Buffers */
    cl_mem colorBuffer, randomBuffer, sphereBuffer, triangleBuffer;
    
    
    /* --------------------------------------------------------------- */

    
    /* Get the Vendor */
    
   err = clGetPlatformIDs(0, NULL, &numPlatforms);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to get OpenCL platforms\n");
        exit(-1);
    }

    
    if (numPlatforms > 0) {
        cl_platform_id *platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id) * numPlatforms);
        err = clGetPlatformIDs(numPlatforms, platforms, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "Failed to get OpenCL platform IDs\n");
            exit(-1);
        }
        
        for (unsigned int i = 0; i < numPlatforms; ++i) {
            char pbuf[100];
            err = clGetPlatformInfo(platforms[i],
                                       CL_PLATFORM_VENDOR,
                                       sizeof(pbuf),
                                       pbuf,
                                       NULL);
            
            err = clGetPlatformIDs(numPlatforms, platforms, NULL);
            if (err != CL_SUCCESS) {
                fprintf(stderr, "Failed to get OpenCL platform IDs\n");
                exit(-1);
            }
            
            fprintf(stderr, "OpenCL Platform %d: %s\n", i, pbuf);
        }
        
        platform = platforms[0];
        free(platforms);
    }
    
    
    /* --------------------------------------------------------------- */

    /* Get the GPU Device (One device in our case) */
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 32, devices, &deviceCount);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to get OpenCL device IDs\n");
        exit(-1);
    }
    
    
    for (unsigned int i = 0; i < deviceCount; ++i) {
        cl_device_type type = 0;
        err = clGetDeviceInfo(devices[i],
                                 CL_DEVICE_TYPE,
                                 sizeof(cl_device_type),
                                 &type,
                                 NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "Failed to get OpenCL device info: %d\n", err);
            exit(-1);
        }
    }
    
    
    for (unsigned int i=0; i<deviceCount; ++i)
    {
        
    
    /* Check the Address Width and the Extensions the Device Support */
    err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME,
                          sizeof(name_data), name_data, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to get OpenCL device info\n");
        exit(-1);
    }
    err = clGetDeviceInfo(devices[i], CL_DEVICE_ADDRESS_BITS,
                          sizeof(ext_data), &addr_data, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to get OpenCL device info\n");
        exit(-1);
    }
    err = clGetDeviceInfo(devices[i], CL_DEVICE_EXTENSIONS,
                          sizeof(ext_data), ext_data, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to get OpenCL device info\n");
        exit(-1);
    }
    
    std::cout<<"NAME :" << name_data<< std::endl;
    std::cout<<"ADDRESS_WIDTH :"<< addr_data << std::endl;
    std::cout<<"EXTENSIONS :"<< ext_data<<std::endl<<std::endl;
        
    }
    /* --------------------------------------------------------------- */


    if (USE_CPU) device = devices[0];
    else device = devices[1];        // OR 2,3,4 ... depending on which GPU to use
    
    
    /* --------------------------------------------------------------- */

    
    
    
    /* Create the Context to hold the device, queue and buffers */
    context = clCreateContext(NULL, 1, &device, NULL,
                              NULL, &err);
    
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to open OpenCL context\n");
        exit(-1);
    }
    
    /* --------------------------------------------------------------- */


    /* Read the Kernel */
    program_handle = fopen("/Users/Master/Desktop/ICG_PROJECT/Cornell_OpenCL/Cornell_OpenCL/radianceGPU.cl", "r");       // YOUR PATH TO KERNEL
    if (!program_handle)
        perror("Failed to open kernel");
    
    fseek(program_handle, 0, SEEK_END);
    program_size = ftell(program_handle);
    rewind(program_handle);
    program_buffer = (char*)malloc(program_size + 1);
    program_buffer[program_size] = '\0';
    fread(program_buffer, sizeof(char), program_size,
          program_handle);
    fclose(program_handle);
    
    
    /* --------------------------------------------------------------- */

    /* Create the Kernel Program */
    program = clCreateProgramWithSource(context, 1,
                                        (const char**)&program_buffer, &program_size, &err);
    
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to open OpenCL kernel sources: %d\n", err);
        exit(-1);
    }
    free(program_buffer);
    
    
    /* --------------------------------------------------------------- */

    /* Build the Kernel Program */
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to build OpenCL kernel: %d\n", err);
    }
    /* If Success, create the Kernel*/
    kernel = clCreateKernel(program, "radianceGPU", &err);
    if (err != CL_SUCCESS) {
            fprintf(stderr, "Failed to create OpenCL kernel: %d\n", err);
            exit(-1);
    }
    /* --------------------------------------------------------------- */

    /* Create the command Queue to hold it */
    commandQueue = clCreateCommandQueue(context, device, 0, &err);
    if (err != CL_SUCCESS) {
            fprintf(stderr, "Failed to create OpenCL command queue: %d\n", err);
            exit(-1);
    }
    
    
    /* --------------------------------------------------------------- */

    /* Create the Data buffers */
    
    
    colors = (cl_float3 *)malloc(sizeof(cl_float3) * pixelCount);
    
    
    /* Random Seeds Buffer */
    
    seeds = (unsigned int *)malloc(sizeof(unsigned int) * pixelCount * 2);
    for (unsigned int i = 0; i < pixelCount * 2; i++) {
        seeds[i] = rand();
        if (seeds[i] < 2)
            seeds[i] = 2;
    }
    
    
    
    colorBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY
                              ,  sizeof(cl_float3) * pixelCount, NULL, &err);
    
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to create OpenCL color buffer: %d\n", err);
        exit(-1);
    }

    
    randomBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY |
                              CL_MEM_COPY_HOST_PTR, sizeof(unsigned int) * pixelCount * 2, seeds, &err);
    
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to create OpenCL random buffer: %d\n", err);
        exit(-1);
    }
    sphereBuffer = clCreateBuffer(
                                    context,
                                    CL_MEM_READ_ONLY,
                                    sizeof(Sphere) * (size_t)sphereCount,
                                    NULL,
                                    &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to create OpenCL scene buffer: %d\n", err);
        exit(-1);
    }
    
    err = clEnqueueWriteBuffer(
                                  commandQueue,
                                  sphereBuffer,
                                  CL_TRUE,
                                  0,
                                  sizeof(Sphere) * (size_t)sphereCount,
                                  spheres,
                                  0,
                                  NULL,
                                  NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to write the OpenCL scene Sphere buffer: %d\n", err);
        exit(-1);
    }
    
    
    triangleBuffer = clCreateBuffer(
                                    context,
                                    CL_MEM_READ_ONLY,
                                    sizeof(Triangle) * (size_t)triangleCount,
                                    NULL,
                                    &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to create OpenCL scene Triangle buffer: %d\n", err);
        exit(-1);
    }
    
    err = clEnqueueWriteBuffer(
                                  commandQueue,
                                  triangleBuffer,
                                  CL_TRUE,
                                  0,
                                  sizeof(Triangle) * (size_t)triangleCount,
                                  triag,
                                  0,
                                  NULL,
                                  NULL);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to write the OpenCL scene Triangle buffer: %d\n", err);
        exit(-1);
    }
    
    /* --------------------------------------------------------------- */

    /* Set the Kernel Arguments */
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem),(void *) &colorBuffer);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. colorBuffer: %d\n", err);
        exit(-1);
    }
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem),(void *) &sphereBuffer);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. sphereBuffer: %d\n", err);
        exit(-1);
    }

    
    
    err = clSetKernelArg(kernel, 2, sizeof(cl_mem),(void *) &triangleBuffer);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. triangleBuffer: %d\n", err);
        exit(-1);
    }

    err = clSetKernelArg(kernel, 3, sizeof(cl_mem),(void *) &randomBuffer);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. randomBuffer: %d\n", err);
        exit(-1);
    }

    
    err = clSetKernelArg(kernel, 4, sizeof(int),(void *) &width);
    
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. width: %d\n", err);
        exit(-1);
    }

    err = clSetKernelArg(kernel, 5, sizeof(int),(void *) &height);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. height: %d\n", err);
        exit(-1);
    }

    err = clSetKernelArg(kernel, 6, sizeof(unsigned int),(void *) &sphereCount);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. sphereCount: %d\n", err);
        exit(-1);
    }

    err = clSetKernelArg(kernel, 7, sizeof(unsigned int),(void *) &triangleCount);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Failed to set OpenCL kernel arg. triangleCount: %d\n", err);
        exit(-1);
    }
    /* --------------------------------------------------------------- */

    /* Execute Kernel */
    const size_t work_units_per_kernel = pixelCount;

    
    std::cout<< "Kernel is being executed " << std::endl;
    //size_t    local_units_per_kernel = 1;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL,
                           &work_units_per_kernel, NULL, 0, NULL, NULL);
    
    
 
    /* Read the Result from the output buffer */
    clEnqueueReadBuffer(commandQueue, colorBuffer, CL_TRUE, 0,
                        sizeof(cl_float3)* pixelCount, colors, 0, NULL, NULL);
    
    
    std::cout<< " "<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout<<"Computation Time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() <<" ms"<<std::endl;
    std::cout<< "Writing to file " << std::endl;

    FILE *f = fopen("/Users/Master/Desktop/ICG_PROJECT/Cornell_OpenCL/image.ppm", "w");         // Write image to PPM file.   // YOUR PATH TO IMAGE
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (int i=0; i<width*height; i++)
        fprintf(f,"%d %d %d ", toInt(colors[i].x), toInt(colors[i].y),toInt(colors[i].z));
    
    
    
    /* --------------------------------------------------------------- */

    /* Always Deallocate the resources used */
    clReleaseMemObject(colorBuffer);
    clReleaseMemObject(randomBuffer);
    clReleaseMemObject(sphereBuffer);
    clReleaseMemObject(triangleBuffer);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commandQueue);
    clReleaseProgram(program);
    clReleaseContext(context);
    clReleaseDevice(device);
    for (unsigned i= 0 ; i < deviceCount; ++i) clReleaseDevice(devices[i]);
    free(seeds);
    free(colors);
    std::cout<< "Resources Deallocated " << std::endl;

    
    /* --------------------------------------------------------------- */

    return 0;
    
}