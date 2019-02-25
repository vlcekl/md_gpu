/*
   File name: inte.cu
   Date:      2009/03/31 22:46
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

#include "gpu_properties.h"
#include "integrator.h"

// config constants
__constant__ float4 box;  // box dimensions
__constant__ float4 boxi; // inverse box dimensions
__constant__ int natom;  // number of atoms
// Integrator constants
__constant__ float dt;   // time step
__constant__ unsigned int enstype;  // ensemble type (0 - NVE, 1 - NVT, 2 - NPT)
__constant__ float temper_0;  // required temperature
__constant__ float tcons;  // thermostat constant
__constant__ float press_0;  // required pressure
__constant__ float pcons;  // barostat constant
__constant__ unsigned int nblocks; // number of blocks

// integrator kernels
#include "nve_leapfrog.cu"
#include "nvt_leapfrog_ber.cu"
#include "npt_leapfrog_ber.cu"
#include "nve_vverlet.cu"

// initializes variables needed by a gpu simulation
gpu_integr* integrate_init( integr* integ, gpu_struct* g_prop, config* conf  )
{
    gpu_integr* g_integ;
    g_integ = (gpu_integr *)malloc(sizeof(gpu_integr));


    // block organization
    g_prop->sum_dimGrid.x = 1; g_prop->sum_dimGrid.y = 1; g_prop->sum_dimGrid.z = 1;
    g_prop->sum_dimBlock.x = g_prop->dimGrid.x; g_prop->sum_dimBlock.y = 1; g_prop->sum_dimBlock.z = 1;
    g_prop->sum_sharedMemSize = g_prop->sharedMemSize;

    // config constants
    cudaMemcpyToSymbol( "natom", &conf->natom, sizeof(int) ); checkCUDAError("natom");
    cudaMemcpyToSymbol( "box", &conf->box, sizeof(float4) ); // may not be constant for npt
    cudaMemcpyToSymbol( "boxi", &conf->boxi, sizeof(float4) ); // may not be constant for npt
    // integrator constants
    cudaMemcpyToSymbol( "dt", &integ->dt, sizeof(float) );
    cudaMemcpyToSymbol( "enstype", &integ->enstype, sizeof(float) );
    if (integ->enstype > 0)  
    {   // NVT
        cudaMemcpyToSymbol( "temper_0", &integ->temper_0, sizeof(float) );
        cudaMemcpyToSymbol( "tcons", &integ->tcons, sizeof(float) ); // thermostat constant
        if (integ->enstype > 1)
        {   // NPT
            cudaMemcpyToSymbol( "press_0", &integ->press_0, sizeof(float) );
            cudaMemcpyToSymbol( "pcons", &integ->pcons, sizeof(float) ); // barostat constant
        }
    }
    cudaMemcpyToSymbol( "nblocks", &g_prop->dimGrid.x, sizeof(unsigned int) ); checkCUDAError("nblocks");

    integ->temper = integ->temper_0;

    size_t float_array_size = sizeof(float) * g_prop->dimGrid.x;
    integ->Ke_sum = (float *)malloc(float_array_size);
    for (int i = 0; i < g_prop->dimGrid.x; i++) { integ->Ke_sum[i] = 0.0f; }
    integ->Wi_sum = (float *)malloc(float_array_size);
    for (int i = 0; i < g_prop->dimGrid.x; i++) { integ->Wi_sum[i] = 0.0f; }
    integ->Ke_total = (float *)malloc(sizeof(float));
    integ->Wi_total = (float *)malloc(sizeof(float));
//    integ->therm = (thermodyn *)malloc(sizeof(thermodyn));

    cudaMalloc( (void **)&g_integ->Ke_sum, float_array_size );
    cudaMalloc( (void **)&g_integ->Wi_sum, float_array_size );
    cudaMalloc( (void **)&g_integ->Ke_total, sizeof(float) );
    cudaMalloc( (void **)&g_integ->Wi_total, sizeof(float) );
    cudaMalloc( (void **)&g_integ->g_therm, sizeof(gpu_thermodyn) );
//    cudaMemcpy( g_integ->g_therm, integ->therm, sizeof(gpu_thermodyn), cudaMemcpyHostToDevice );
    cudaMemcpy( g_integ->Ke_sum, integ->Ke_sum, float_array_size, cudaMemcpyHostToDevice );

    return g_integ;
}

// Copy particle property arrays to device
void integrate_HostToDevice (integr* integ, gpu_integr* g_integ, gpu_struct* g_prop)
{
    size_t float_array_size = sizeof(float) * g_prop->dimGrid.x;
    cudaMemcpy( g_integ->Ke_sum, integ->Ke_sum, float_array_size, cudaMemcpyHostToDevice );
    cudaMemcpy( g_integ->Wi_sum, integ->Wi_sum, float_array_size, cudaMemcpyHostToDevice );
    cudaMemcpy( g_integ->Ke_total, integ->Ke_total, sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( g_integ->Wi_total, integ->Wi_total, sizeof(float), cudaMemcpyHostToDevice );
}

void integrate_DeviceToHost (integr* integ, gpu_integr* g_integ, gpu_struct* g_prop)
{
    
    size_t float_array_size = sizeof(float) * g_prop->dimGrid.x;
    cudaMemcpy( integ->Ke_sum,   g_integ->Ke_sum,   float_array_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( integ->Wi_sum,   g_integ->Wi_sum,   float_array_size, cudaMemcpyDeviceToHost );
//    cudaMemcpy( integ->Ke_total, g_integ->Ke_total, float4_array_size, cudaMemcpyDeviceToHost );
//    cudaMemcpy( integ->Wi_total, g_integ->Wi_total, float4_array_size, cudaMemcpyDeviceToHost );
}

void integrate_finish(integr* integ, gpu_integr* g_integ)
{
    cudaFree( g_integ->Ke_sum );
    cudaFree( g_integ->Wi_sum );
    cudaFree( g_integ->Ke_total );
    cudaFree( g_integ->Wi_total );
    cudaFree( g_integ->g_therm );
}

// integrate() advances positions and velocities, and calculates tempareture +
// scaling factor for Berendsen thermostat
// Most of the work is done in kernels
void integrate( gpu_config* g_conf, gpu_integr* g_integ, gpu_struct* g_prop, integr* integ, config* conf)
{
    float scale;
    float tco;
//        size_t float4_array_size = sizeof(float4) * conf->natom;
//        cudaMemcpy( conf->force, g_conf->force, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "ifor1: %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w);
//        cudaMemcpy( conf->pos, g_conf->pos, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "ipos1: %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z, conf->pos[0].w);
    switch (integ->enstype)
    {
        case 0:
            leapfrog_nve <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->pos, g_conf->vel, g_conf->force );
            break;
        case 1:
            tco = 0.01;
            scale = sqrt(1.0f + integ->dt*tco*(integ->temper_0/integ->temper - 1.0f));
            leapfrog_nvt_berendsen <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->pos, g_conf->vel, g_conf->force, scale, g_integ->Ke_sum );
//            checkCUDAError("nvt");
            kinetic_energy_sum <<< g_prop->sum_dimGrid, g_prop->sum_dimBlock, g_prop->sum_sharedMemSize >>> ( g_integ->Ke_sum);
//            checkCUDAError("k_sum");
            cudaMemcpy( integ->Ke_sum, g_integ->Ke_sum, sizeof(float)*g_prop->dimGrid.x, cudaMemcpyDeviceToHost );
//            checkCUDAError("memcp");
            integ->Ke_sum[0] *= 0.5f;
            integ->temper = integ->Ke_sum[0]*2.0f/(3.0f*(float)conf->natom - 3.0f)*1000.0/8.314;
            break;
//        case 2: // NPT Leapfrog integrator with Berendsen thermostat (scale = 0 -> NVE)
//            float tco = 0.1;
//            float scale = sqrt(1.0f + dt*tco*(temper_0/integ->temper - 1.0f));
//            leapfrog_npt_berendsen <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->pos, g_conf->vel, g_conf->force, scale, g_integ->Ke_sum );
//            kinetic_energy_and_virial_sum <<< g_prop->sum_dimGrid, g_prop->sum_dimBlock, g_prop->sum_sharedMemSize >>> ( g_integ->Ke_sum, g_integ->Ke_total);
//            cudaMemcpy( integ->Ke_sum, g_integ->Ke_sum, sizeof(float)*g_prop->dimGrid.x, cudaMemcpyDeviceToHost ); checkCUDAError("k_sum");
//            integ->temper = integ->Ke_sum[0]*2.0f/(3.0f*(float)conf->natom - 3.0f)*1000.0/8.314;
//            break;
        default:
            printf("No such enseble implemented\n");
    }

//        cudaMemcpy( conf->force, g_conf->force, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "ifor2: %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w);
//        cudaMemcpy( conf->pos, g_conf->pos, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "ipos2: %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z, conf->pos[0].w);
}

void integrate_pre( gpu_config* g_conf, gpu_integr* g_integ, gpu_struct* g_prop, integr* integ, config* conf)
{
//    float scale;
//   float tco;
    switch (integ->enstype)
    {
        case 0:
            vverlet_nve_pre <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->pos, g_conf->vel, g_conf->force );
            checkCUDAError("nve_vverlet");
            break;
        case 1:
//            tco = 1.0;
//            scale = sqrt(1.0f + integ->dt*tco*(integ->temper_0/integ->temper - 1.0f));
//            leapfrog_nvt_berendsen <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->pos, g_conf->vel, g_conf->force, scale, g_integ->Ke_sum );
//            kinetic_energy_sum <<< g_prop->sum_dimGrid, g_prop->sum_dimBlock, g_prop->sum_sharedMemSize >>> ( g_integ->Ke_sum);
//            cudaMemcpy( integ->Ke_sum, g_integ->Ke_sum, sizeof(float)*g_prop->dimGrid.x, cudaMemcpyDeviceToHost );
//            integ->Ke_sum[0] *= 0.5f;
//            integ->temper = integ->Ke_sum[0]*2.0f/(3.0f*(float)conf->natom - 3.0f)*1000.0/8.314;
//            break;
        default:
            printf("No such enseble implemented\n");
    }

}

void integrate_cor( gpu_config* g_conf, gpu_integr* g_integ, gpu_struct* g_prop, integr* integ, config* conf)
{
    float scale;
    float tco;
    switch (integ->enstype)
    {
        case 0:
            vverlet_nve_cor <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->vel, g_conf->force );
            break;
        case 1:
            tco = 1.0;
            scale = sqrt(1.0f + integ->dt*tco*(integ->temper_0/integ->temper - 1.0f));
            leapfrog_nvt_berendsen <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->pos, g_conf->vel, g_conf->force, scale, g_integ->Ke_sum );
            kinetic_energy_sum <<< g_prop->sum_dimGrid, g_prop->sum_dimBlock, g_prop->sum_sharedMemSize >>> ( g_integ->Ke_sum);
            cudaMemcpy( integ->Ke_sum, g_integ->Ke_sum, sizeof(float)*g_prop->dimGrid.x, cudaMemcpyDeviceToHost );
            integ->Ke_sum[0] *= 0.5f;
            integ->temper = integ->Ke_sum[0]*2.0f/(3.0f*(float)conf->natom - 3.0f)*1000.0/8.314;
            break;
        default:
            printf("No such enseble implemented\n");
    }
}

/* end of integrator.cu */
