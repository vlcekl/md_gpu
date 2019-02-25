/*
   File name: gpu_properties.cu
   Date:      2009/03/31 22:39
   Author:    Aaron Thompson and Lukas Vlcek

   Copyright (C) 2009 Aaron Thompson and Lukas Vlcek

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
#include <cuda.h>

#include "gpu_properties.h"


// related to DEVICE but allocated on HOST
gpu_struct* gpu_properties_init( config* conf )
{
    gpu_struct* g_prop;
    g_prop = (gpu_struct *)malloc(sizeof(gpu_struct));

    cudaGetDeviceCount(&g_prop->numdevices);

    cudaGetDevice(&g_prop->thisdevice);

    g_prop->prop = (cudaDeviceProp *)malloc(sizeof(cudaDeviceProp));

    cudaGetDeviceProperties(g_prop->prop, g_prop->thisdevice);

    if (g_prop->prop->major == 1 && g_prop->prop->minor == 1) {
        g_prop->block_size = 128;
    } else {
        printf("Not tested for compute capability other than 1.1\n"); // should be an error message
        g_prop->block_size = 256;
    }

    g_prop->dimBlock.x = g_prop->block_size;	// in nbody demo p=256, q=1
    g_prop->dimBlock.y = 1;	// in nbody demo p=256, q=1
    g_prop->dimBlock.z = 1;	// in nbody demo p=256, q=1
    g_prop->dimGrid.x = conf->natom/g_prop->block_size + ( !(conf->natom % g_prop->block_size)?0:1);
    g_prop->dimGrid.y = 1;	// in nbody demo p=256, q=1
    g_prop->dimGrid.z = 1;	// in nbody demo p=256, q=1
    g_prop->sharedMemSize = g_prop->block_size * sizeof(float4); // 4 floats for position

    {
        printf( "block_size: %i\n", g_prop->block_size );
        printf( "dimGrid:   (%i, %i, %i)\n", g_prop->dimGrid.x, g_prop->dimGrid.y, g_prop->dimGrid.z );
        printf( "dimBlock:  (%i, %i, %i)\n", g_prop->dimBlock.x, g_prop->dimBlock.y, g_prop->dimBlock.z );
        printf( "shareMem:   %i\n", g_prop->sharedMemSize );
    }

    return g_prop;
}

void gpu_properties_finish( gpu_struct* g_prop )
{
    free(g_prop->prop);
    free(g_prop);
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
} 

/* end of gpu_properties.cu */
