/*
   File name: test_device.cu
   Date:      2009/03/31 23:53
   Author:    Lukas Vlcek

   Copyright (C) 2009 Lukas Vlcek

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   in a file called COPYING along with this program; if not, write to
   the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
   02139, USA.
*/

#include <stdio.h>

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
} 

int main( int argc, char** argv)
{
    int count;

    cudaGetDeviceCount(&count);
    checkCUDAError("cudaGetDeviceCount");
    printf("Device count: %d\n", count);

    cudaGetDevice(&count);
    checkCUDAError("cudaGetDevice");
    printf("Active device: %d\n", count);

    cudaDeviceProp* prop;
    prop = (cudaDeviceProp *)malloc(sizeof(cudaDeviceProp));
    cudaGetDeviceProperties(prop, count);
    checkCUDAError("cudaGetDeviceProperties");

    printf("Device no.%d properties:\n", count);
    printf("Name: %s\n", prop->name); 
    printf("GlobMem: %d\n", prop->totalGlobalMem); 
    printf("ShMem/block: %d\n", prop->sharedMemPerBlock); 
    printf("Regs/block: %d\n", prop->regsPerBlock); 
    printf("WarpSize: %d\n", prop->warpSize); 
    printf("memPitch: %d\n", prop->memPitch); 
    printf("maxThreadsPerBlock: %d\n", prop->maxThreadsPerBlock); 
    printf("maxThreadsDim: %d %d %d\n", prop->maxThreadsDim[0], prop->maxThreadsDim[1], prop->maxThreadsDim[2]); 
    printf("maxGridSize: %d %d %d\n", prop->maxGridSize[0], prop->maxGridSize[1], prop->maxGridSize[2]); 
    printf("ConstMem: %d\n", prop->totalConstMem); 
    printf("Compute capability: %d.%d\n", prop->major, prop->minor); 
//    printf("Major: %d\n", prop->minor); 
//    printf("Minor: %d\n", prop->minor); 
    printf("clockRate: %d\n", prop->clockRate); 
    printf("textureAlignment: %d\n", prop->textureAlignment); 
    printf("deviceOverlap: %d\n", prop->deviceOverlap); 
    printf("multiProcessorCount: %d\n", prop->multiProcessorCount);

    return 0;
}

/* end of test_device.cu */
