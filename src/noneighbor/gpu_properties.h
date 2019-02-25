/*
   File name: gpu_properties.h
   Date:      2009/03/31 22:40
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


#ifndef __GPU_PROPERTIES_H__
#define __GPU_PROPERTIES_H__

#include <cuda.h>

#include "config.h"

typedef struct {
    int numdevices;
    int thisdevice;
    int block_size;
    dim3 dimGrid;
    dim3 dimBlock;
    size_t sharedMemSize;
    dim3 sum_dimGrid;
    dim3 sum_dimBlock;
    size_t sum_sharedMemSize;
    unsigned int nblocks;
    cudaDeviceProp *prop;
} gpu_struct;

gpu_struct* gpu_properties_init( config* conf );
void gpu_properties_finish( gpu_struct* g_prop );
void checkCUDAError(const char *msg);

#endif

/* end of gpu_properties.h */
