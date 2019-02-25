/*
   File name: pair_forces.cu
   Date:      2009/04/01 01:17
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

#ifdef __DEVICE_EMULATION__
#include <stdio.h>
#endif

#include "config.h"
#include "pair_forces.h"
#include "gpu_properties.h"

// Macros to simplify shared memory addressing
#define SX(i) sharedPos[i+blockDim.x*threadIdx.y]
// This macro is only used when multithreadBodies is true (below)
#define SX_SUM(i,j) sharedPos[i+blockDim.x*j]
// WRAP is used to force each block to start working on a different 
// chunk (and wrap around back to the beginning of the array) so that
// not all multiprocessors try to read the same memory locations at once.
#define WRAP(x,m) (((x)<m)?(x):(x-m))  // Mod without divide, works on values from 0 up to 2m

__constant__ int natoms;  // number of atoms
__constant__ int ntype;  // number of atomtypes
__constant__ float4 fbox;  // box dimensions
__constant__ float4 fboxi; // inverse box dimensions
__constant__ float rcutsq;  // cutoff distance squared
__constant__ float reacf;  // cutoff distance squared

texture<float4, 1, cudaReadModeElementType> pos_tex;

gpu_ffield* pair_forces_init(ffield* fld, params* pars, config* conf, gpu_struct* g_prop )
{
    gpu_ffield* g_fld;
    g_fld = (gpu_ffield *)malloc(sizeof(gpu_ffield));

    // config constants
    cudaMemcpyToSymbol( "natoms", &conf->natom, sizeof(int) );
    cudaMemcpyToSymbol( "fbox", &conf->box, sizeof(float4) ); // may not be constant for npt
    cudaMemcpyToSymbol( "fboxi", &conf->boxi, sizeof(float4) ); // may not be constant for npt
    cudaMemcpyToSymbol( "ntype", &fld->ntype, sizeof(int) );
    // ffield constants
    cudaMemcpyToSymbol( "rcutsq", &pars->rcutsq, sizeof(float) );
    cudaMemcpyToSymbol( "reacf", &pars->reacf, sizeof(float) ); checkCUDAError("lj0"); 
    cudaMemcpyToSymbol( lj1, fld->lj1, NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES*sizeof(float), 0 ); checkCUDAError("lj1");
    cudaMemcpyToSymbol( lj2, fld->lj2, NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES*sizeof(float), 0 ); checkCUDAError("lj2");

    return g_fld;
}

void pair_forces_finish(ffield* fld, gpu_ffield* g_fld)
{
}

// Copy particle property arrays to device
void pair_forces_HostToDevice (ffield* fld, gpu_ffield* g_fld )
{
}

void pair_forces_DeviceToHost (ffield* fld, gpu_ffield* g_fld)
{
}

// Computes the forces between one particle and all others
__global__ void pair_forces( float4* force, float4* pos )
{
    extern __shared__ float4 sharedPos[];

    // Load LJ parameters to shared memory in a way to prevent bank conflicts
    for (int cur_offset = 0; cur_offset < coeff_width*coeff_width; cur_offset += blockDim.x)
        {
            if (cur_offset + threadIdx.x < coeff_width*coeff_width)
                s_coeffs[cur_offset + threadIdx.x] = d_coeffs[cur_offset + threadIdx.x];
        }
    __syncthreads();

    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    float4 posi = tex1Dfetch(pos_tex, i); // float4 posi = pos[i];
    float4 f = {0.0f, 0.0f, 0.0f, 0.0f}; //make_float4(rcutsq, (float)i, (float)natoms, 4.0f);

    unsigned int it = ntype*0;

    // Each thread loads a particle position from global memory into  shared // load a neighborlist?
    sharedPos[threadIdx.x + blockDim.x*threadIdx.y] = pos[WRAP(blockIdx.x + tile, gridDim.x) * blockDim.x + threadIdx.x];

    for(unsigned int counter=0; counter < blockDim.x; counter++)
    {
        if( i < natoms )
            return;

        float4 posj = tex1Dfetch(pos_tex, j); // float4 posi = pos[i];

        // Distance for the minimum image
        posj.x = posi.x - posj.x;
        posj.x -= (fbox.x) * rintf(posj.x*fboxi.x);
        posj.y = posi.y - posj.y;
        posj.y -= (fbox.y) * rintf(posj.y*fboxi.y);
        posj.z = posi.z - posj.z;
        posj.z -= (fbox.z) * rintf(posj.z*fboxi.z);

        // Square and inverse distance 
        float rsq = posj.x*posj.x + posj.y*posj.y + posj.z*posj.z;
        float ri = rsqrtf(rsq);

        // Forces + virial
        float force_added = qfac*posj.w*posi.w*(ri - reacf*rsq); // Coulomb
        ri *= ri; ri *= ri; ri *= ri;  // r^-8
        unsigned int jt = it + 0;      // particle-particle type (only 0 for now)
        force_added += (lj1[jt] * ri * rsq - lj2[jt]) * ri; // 12-6 VdW // consider putting lj to shared memory
        if (rsq > rcutsq) { force_added = 0.0f; }// no self-interaction - not needed for neighbor list

        f.w += force_added * rsq;     // virial*6
        f.x += force_added * posj.x;  // forces
        f.y += force_added * posj.y;
        f.z += force_added * posj.z;
    }
    force[i] = f;
}

