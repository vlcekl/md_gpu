/*
   File name: pair_forces.h
   Date:      2009/04/01 01:18
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


#ifndef __PAIR_FORCES_H__
#define __PAIR_FORCES_H__

#include <vector_types.h>
#include "simulation.h"
#include "gpu_properties.h"
#include "config.h"

// LJ parameters & atomic masses
#define NUM_PARTICLE_TYPES      8
__constant__ float lj1[NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES];
__constant__ float lj2[NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES];

// General constant parameters
__constant__ float qfac = 1389.35427611722f;

typedef struct {
    float something;
//    float2* lj;
} gpu_ffield;

// pair forces
gpu_ffield* pair_forces_init(ffield* fld, params* pars, config* conf, gpu_struct* g_prop );
void pair_forces_finish(ffield* fld, gpu_ffield* g_fld);
void pair_forces_HostToDevice (ffield* fld, gpu_ffield* g_fld );
void pair_forces_DeviceToHost (ffield* fld, gpu_ffield* g_fld);

// Kernel
__global__ void pair_forces( float4* pos, float4* force );

#endif

/* end of pair_forces.h */

