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

#include "simulation.h"
#include "gpu_properties.h"
#include "config.h"
#include "neighborlist.h"
#include "forces.h"
#include "integrator.h"
#include "prints.h"

// General configuration parameters
__constant__ int natoms;  // number of atoms
__constant__ int ntype;  // number of atomtypes
// NVT only (what to do wiht NPT?)
__constant__ float4 fbox;  // box dimensions
__constant__ float4 fboxi; // inverse box dimensions
// Neighborlist parameters
__constant__ int pitch_f;
// Nonbonded interactions
__constant__ float rcutsq;  // cutoff distance squared
__constant__ float reacf;  // cutoff distance squared
__constant__ float qfac = 1389.35427611722f;
__constant__ float lj1[NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES];
__constant__ float lj2[NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES];
// Ewald (erf) parameters
__constant__ float a1 = 0.254829592f;  
__constant__ float a2 = -0.284496736f; 
__constant__ float a3 = 1.421413741f;  
__constant__ float a4 = -1.453152027f; 
__constant__ float a5 = 1.061405429f;  
__constant__ float palfa;     
__constant__ float alfasq;     
__constant__ float spiialfa;     
__constant__ float wolf;     
// Bonded interactions
__constant__ float bnd1[NUM_BOND_TYPES];
__constant__ float bnd2[NUM_BOND_TYPES];
__constant__ float ang1[NUM_BOND_TYPES];
__constant__ float ang2[NUM_BOND_TYPES];

// Texture definition
texture<float4, 1, cudaReadModeElementType> posTex;
//texture<float4, 1, cudaReadModeElementType> forTex;

// Include forces files
#include "pair_forces.cu" // pair interactions reaction field and real Ewald
#include "bonded_forces.cu" // harmonic bonds and harmonic angles
#include "energy.cu" // harmonic bonds and harmonic angles

gpu_ffield* forces_init(ffield* fld, params* pars, config* conf, gpu_struct* g_prop, gpu_config* g_conf)
{
    gpu_ffield* g_fld;
    g_fld = (gpu_ffield *)malloc(sizeof(gpu_ffield));

    // config constants
    cudaMemcpyToSymbol( "natoms", &conf->natom, sizeof(int) );
    cudaMemcpyToSymbol( "fbox", &conf->box, sizeof(float4) ); // may not be constant for npt
    cudaMemcpyToSymbol( "fboxi", &conf->boxi, sizeof(float4) ); // may not be constant for npt
    cudaMemcpyToSymbol( "ntype", &fld->ntype, sizeof(int) );

    // nnlist parameters
    int apitch = g_prop->dimGrid.x*g_prop->dimBlock.x;
    cudaMemcpyToSymbol( "pitch_f", &apitch, sizeof(int) ); checkCUDAError("pitch_f"); 

    // Nonbonded interactions
    cudaMemcpyToSymbol( "rcutsq", &pars->rcutsq, sizeof(float) );
    cudaMemcpyToSymbol( "reacf", &pars->reacf, sizeof(float) ); checkCUDAError("lj0"); 
    cudaMemcpyToSymbol( lj1, fld->lj1, NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES*sizeof(float), 0 ); checkCUDAError("lj1");
    cudaMemcpyToSymbol( lj2, fld->lj2, NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES*sizeof(float), 0 ); checkCUDAError("lj2");

//  Ewald & Wolf parameters
    if (pars->lkspace > 0) {
        float palfa = pars->alfa*Pex;
        cudaMemcpyToSymbol( "palfa", &palfa, sizeof(float) ); checkCUDAError("palfa");
        float spiialfa = 2.0f/sqrtf(M_Pi)*pars->alfa;
        cudaMemcpyToSymbol( "spiialfa", &spiialfa, sizeof(float) ); checkCUDAError("spiialfa");
        float alfasq = -pars->alfa*pars->alfa/logf(2.0f);
        cudaMemcpyToSymbol( "alfasq", &alfasq, sizeof(float) ); checkCUDAError("alfasq");
        float wolf = 0.0f;
        if (pars->lkspace == 1) {
            float ff = 1.0f/(1.0f + palfa*sqrt(pars->rcutsq));
            wolf = spiialfa*sqrtf(pars->rcutsq) + ((((A5*ff + A4)*ff + A3)*ff + A2)*ff + A1)*ff; 
            wolf *= exp2f(alfasq*pars->rcutsq)/sqrt(pars->rcutsq);
        }
        cudaMemcpyToSymbol( "wolf", &wolf, sizeof(float) ); checkCUDAError("wolf");
    }
    
    // Harmonic Bonds 
    if (fld->nbonds > 0) {
        cudaMemcpyToSymbol( bnd1, fld->bnd1, NUM_BOND_TYPES*sizeof(float), 0 ); checkCUDAError("bnd1");
        cudaMemcpyToSymbol( bnd2, fld->bnd2, NUM_BOND_TYPES*sizeof(float), 0 ); checkCUDAError("bnd2");

        fld->nmaxbond = 4;
        size_t uint_array_size = sizeof(unsigned int) * g_prop->dimGrid.x*g_prop->dimBlock.x;
        fld->bndlist_max = (unsigned int *)malloc( uint_array_size );
        size_t uint2_array_size = sizeof(uint2) * g_prop->dimGrid.x*g_prop->dimBlock.x*fld->nmaxbond;
        fld->bndlist = (uint2 *)malloc( uint2_array_size );
 
        for (int i = 0; i < g_prop->dimGrid.x*g_prop->dimBlock.x; i++) { fld->bndlist_max[i] = 0; }
        for (int it = 0; it < fld->nbonds; it++) { 
            unsigned int i = fld->iatn[it].x-1;
            unsigned int j = fld->iatn[it].y-1;
            fld->bndlist[i + fld->bndlist_max[i]*apitch].x = j; 
            fld->bndlist[i + fld->bndlist_max[i]*apitch].y = fld->iatn[it].w-1; 
            fld->bndlist[j + fld->bndlist_max[j]*apitch].x = i; 
            fld->bndlist[j + fld->bndlist_max[j]*apitch].y = fld->iatn[it].w-1; 
            fld->bndlist_max[i]++;
            fld->bndlist_max[j]++;
        }
 
        cudaMalloc( (void **)&g_fld->bndlist_max, uint_array_size );
        cudaMalloc( (void **)&g_fld->bndlist, uint2_array_size );
        cudaMemcpy( g_fld->bndlist_max, fld->bndlist_max, uint_array_size, cudaMemcpyHostToDevice );
        cudaMemcpy( g_fld->bndlist, fld->bndlist, uint2_array_size, cudaMemcpyHostToDevice );
    }

    // Harmonic angles
    if (fld->nangle > 0) {
        cudaMemcpyToSymbol( ang1, fld->ang1, NUM_BOND_TYPES*sizeof(float), 0 ); checkCUDAError("ang1");
        cudaMemcpyToSymbol( ang2, fld->ang2, NUM_BOND_TYPES*sizeof(float), 0 ); checkCUDAError("ang2");
 
        fld->nmaxbond = 4;
        size_t int_array_size = sizeof(unsigned int) * g_prop->dimGrid.x*g_prop->dimBlock.x;
        fld->anglist_max = (unsigned int *)malloc( int_array_size );
        size_t uint4_array_size = sizeof(uint4) * g_prop->dimGrid.x*g_prop->dimBlock.x*fld->nmaxbond;
        fld->anglist = (uint4 *)malloc( uint4_array_size );
 
        for (int i = 0; i < g_prop->dimGrid.x*g_prop->dimBlock.x; i++) { fld->anglist_max[i] = 0; }
        for (int it = 0; it < fld->nangle; it++) { 
            unsigned int i = fld->iang[it].x-1;
            unsigned int j = fld->iang[it].y-1;
            unsigned int k = fld->iang[it].z-1;
            fld->anglist[i + fld->anglist_max[i]*apitch].x = 0; 
            fld->anglist[i + fld->anglist_max[i]*apitch].y = j; 
            fld->anglist[i + fld->anglist_max[i]*apitch].z = k; 
            fld->anglist[i + fld->anglist_max[i]*apitch].w = fld->iang[it].w-1; 
            fld->anglist[j + fld->anglist_max[j]*apitch].x = 1; 
            fld->anglist[j + fld->anglist_max[j]*apitch].y = i; 
            fld->anglist[j + fld->anglist_max[j]*apitch].z = k; 
            fld->anglist[j + fld->anglist_max[j]*apitch].w = fld->iang[it].w-1; 
            fld->anglist[k + fld->anglist_max[k]*apitch].x = 0; 
            fld->anglist[k + fld->anglist_max[k]*apitch].y = j; 
            fld->anglist[k + fld->anglist_max[k]*apitch].z = i; 
            fld->anglist[k + fld->anglist_max[k]*apitch].w = fld->iang[it].w-1; 
            fld->anglist_max[i]++;
            fld->anglist_max[j]++;
            fld->anglist_max[k]++;
        }
 
        cudaMalloc( (void **)&g_fld->anglist_max, int_array_size );
        cudaMalloc( (void **)&g_fld->anglist, uint4_array_size );
        cudaMemcpy( g_fld->anglist_max, fld->anglist_max, int_array_size, cudaMemcpyHostToDevice );
        cudaMemcpy( g_fld->anglist, fld->anglist, uint4_array_size, cudaMemcpyHostToDevice );
    }

    // Positions in texture
    cudaBindTexture(0, posTex, g_conf->pos, sizeof(float4)*conf->natom); checkCUDAError("texture"); 
//    cudaBindTexture(0, forTex, g_conf->force, sizeof(float4)*conf->natom); checkCUDAError("texture"); 

    return g_fld;
}

// run pair, bonded, angle forces
void forces(gpu_struct* g_prop, gpu_nlist* g_nlist, gpu_config* g_conf, gpu_ffield* g_fld, ffield* fld, config* conf, params* pars)
{
    // real space pair interactions
    if (pars->lkspace == 0) {  // reaction field
        pair_forces <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->force, g_nlist->nnlist, g_nlist->nnlist_max );
    } else {  // real Ewald (Wolf)
        pair_forces_ewa <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->force, g_nlist->nnlist, g_nlist->nnlist_max );
    }
//    print_ave_force_momentum(g_conf, conf);

    // harmonic bonds
    if (fld->nbonds > 0) {
        bonded_forces <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->force, g_fld->bndlist, g_fld->bndlist_max );
    }
//    print_ave_force_momentum(g_conf, conf);

    // harmonic angles
    if (fld->nangle > 0) {
        angle_forces <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->force, g_fld->anglist, g_fld->anglist_max);
    }
//    print_ave_force_momentum(g_conf, conf);
}

void forces_finish(ffield* fld, gpu_ffield* g_fld)
{
    cudaFree( g_fld->bndlist );
    cudaFree( g_fld->bndlist_max );
    cudaFree( g_fld->anglist );
    cudaFree( g_fld->anglist_max );
}

// Copy particle property arrays to device
void forces_HostToDevice (ffield* fld, gpu_ffield* g_fld )
{
}

void forces_DeviceToHost (ffield* fld, gpu_ffield* g_fld)
{
}

void energy(gpu_struct* g_prop, gpu_nlist* g_nlist, gpu_config* g_conf, gpu_ffield* g_fld, ffield* fld, config* conf, params* pars)
{
}
