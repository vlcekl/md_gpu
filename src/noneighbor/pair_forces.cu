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
__constant__ unsigned int numTiles;  // number of atoms

texture<float4, 1, cudaReadModeElementType> posDTex;
//texture<unsigned int, 1, cudaReadModeElementType> neighborlistDTex;
//extern __shared__ float4 sharedPos[];

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

    unsigned int numi = conf->natom / g_prop->block_size;
    cudaMemcpyToSymbol( "numTiles", &numi, sizeof(unsigned int) );


    cudaBindTexture(0, posDTex, conf->pos, sizeof(float4)*conf->natom);

//    { // Ewald parameters
//        float faux;
//        faux = alfa*p;
//        cudaMemcpyToSymbol( "palfa", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "spiialfa", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "alfasq", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "a1", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "a2", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "a3", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "a4", &pars->rcutsq, sizeof(float) );
//        cudaMemcpyToSymbol( "a5", &pars->rcutsq, sizeof(float) );
//    }

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
//__global__ void pair_forces( float4* force ) // maybe enough
__global__ void pair_forces( float4* force, float4* pos )
{
    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if( i < natoms )
    {
	float4 posi = tex1Dfetch(posDTex, i);
        float4 posj;
	int posTexIndex;

	unsigned int num_neighbors = neighborlistD[i*nnlist_max];

        unsigned int it = ntype*0;

        //sharedPos[threadIdx.x + blockDim.x*threadIdx.y] = pos[WRAP(blockIdx.x + tile, gridDim.x) * blockDim.x + threadIdx.x];

        float4 f = {0.0f, 0.0f, 0.0f, 0.0f}; //make_float4(rcutsq, (float)i, (float)natoms, 4.0f);
  	for(unsigned int counter=0; counter < blockDim.x; counter++)
  	{
            posDTexIndex = neighborlistD[i*NEIGHBORLIST_CAPACITY + j];
            posj = tex1Dfetch(posDTex, posTexIndex);

            unsigned int jt = it + 0;      // particle-particle type (only 0 for now)

            float lja = lj1[jt];
            float ljc = lj2[jt];

            // Distance for the minimum image
            posj.x = posi.x - posj.x;
            posj.y = posi.y - posj.y;
            posj.z = posi.z - posj.z;

            posj.x -= fbox.x * rintf(posj.x*fboxi.x);
            posj.y -= fbox.y * rintf(posj.y*fboxi.y);
            posj.z -= fbox.z * rintf(posj.z*fboxi.z);

            // Square and inverse distance 
            float rsq = posj.x*posj.x + posj.y*posj.y + posj.z*posj.z; // r^2
            float ri = rsqrtf(rsq); // r^-1

            // Forces + virial
            float force_added = qfac*posj.w*posi.w*(ri - reacf*rsq); // Coulomb
            ri *= ri;            // r^-2
            force_added *= ri;
            ri *= ri; ri *= ri;  // r^-8
            force_added += (lja * ri * rsq - ljc) * ri; // 12-6 VdW // consider putting lj to shared memory

            f.w += force_added * rsq;     // virial*6
            f.x += force_added * posj.x;  // forces
            f.y += force_added * posj.y;
            f.z += force_added * posj.z;
  	}
        //    __syncthreads();
        force[i] = f;
    }	// end if i < natoms
}

//// Computes the forces between one particle and all others
//__global__ void pair_forces_ewald( float4* force, float4* pos )
//{
//    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
//
//    if( i < natoms )
//    {
//        extern __shared__ float4 sharedPos[];
//	float4 posi = pos[i];
//        float4 f = {0.0f, 0.0f, 0.0f, 0.0f}; //make_float4(rcutsq, (float)i, (float)natoms, 4.0f);
//
//        unsigned int it = ntype*0;
//
////        unsigned int numTiles = natoms / (blockDim.x * blockDim.y);
//
//        for(unsigned int tile = blockIdx.y; tile < numTiles + blockIdx.y; tile++)
//	{
//            // Each thread loads a particle position from global memory into  shared // load a neighborlist?
//            sharedPos[threadIdx.x + blockDim.x*threadIdx.y] = pos[WRAP(blockIdx.x + tile, gridDim.x) * blockDim.x + threadIdx.x];
//
//            __syncthreads();
//
//  	    for(unsigned int counter=0; counter < blockDim.x; counter++)
//  	    {
//                float4 posj = SX(counter);
//
//                // Distance for the minimum image
//                posj.x = posi.x - posj.x;
//                posj.x -= (fbox.x) * rintf(posj.x*fboxi.x);
//                posj.y = posi.y - posj.y;
//                posj.y -= (fbox.y) * rintf(posj.y*fboxi.y);
//                posj.z = posi.z - posj.z;
//                posj.z -= (fbox.z) * rintf(posj.z*fboxi.z);
//
//                // Square and inverse distance 
//                float rsq = posj.x*posj.x + posj.y*posj.y + posj.z*posj.z;
//                if (rsq < 0.0000000001f) { rsq = 2.0f*rcutsq; }// no self-interaction - not needed for neighbor list
//                float ri = rsqrtf(rsq); // r^-1
//                float r = rsq*ri; // r
//
//                // Real Ewald part
//                float t = __fdividef(1.0f,1.0f + palfa*r); // r
//                float force_added = spiialfa*((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*r + 1.0f;
//                force_added *= exp2f(-alfasq*rsq);
//                force_added *= qfac_ew*posi.w*posj.w*ri;
//                ri *= ri; // r^-2
//                force_added *= ri;
//                // LJ part
//                ri *= ri; ri *= ri;  // r^-8
//                unsigned int jt = it + 0;      // particle-particle type (only 0 for now)
//                force_added += (lj1[jt] * ri * rsq - lj2[jt]) * ri; // 12-6 VdW // consider putting lj to shared memory
//
//                // Forces + virial
//                f.w += force_added * rsq;     // virial*6
//                f.x += force_added * posj.x;  // forces
//                f.y += force_added * posj.y;
//                f.z += force_added * posj.z;
//  	    }
//
//            __syncthreads();
//
//        }
//        force[i] = f;
//    }	// end if i < natoms
//}

