/*
 *  runSimulation.cu
 *  SimpleMD
 *
 *  Created by Aaron Thompson on 6/16/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
// Fastest GPU runtime with BLOCK_SIZE 32 and particles struct align(4): 2724
// Fastest CPU runtime: 28519

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>

#include "gpu_properties.h"
#include "pair_forces.h"
#include "config.h"

#pragma mark Simulation functions

gpu_config* config_init( config* conf )
{

    gpu_config* g_conf;
    g_conf = (gpu_config *)malloc(sizeof(gpu_config));
 
    size_t float4_array_size = sizeof(float4) * conf->natom;
    cudaMalloc( (void **)&g_conf->pos, float4_array_size );
    cudaMalloc( (void **)&g_conf->vel, float4_array_size );
    cudaMalloc( (void **)&g_conf->force, float4_array_size );

    config_HostToDevice(conf, g_conf);

    return g_conf;
}

void config_finish( config* conf, gpu_config* g_conf )
{
    config_DeviceToHost(conf, g_conf, 2);

    // Free allocated device memory
    cudaFree( g_conf->pos );
    cudaFree( g_conf->vel );
    cudaFree( g_conf->force );
}

// Copy particle property arrays to device
void config_HostToDevice (config* conf, gpu_config* g_conf )
{
    size_t float4_array_size = sizeof(float4) * conf->natom;
    cudaMemcpy( g_conf->pos, conf->pos, float4_array_size, cudaMemcpyHostToDevice );
    cudaMemcpy( g_conf->vel, conf->vel, float4_array_size, cudaMemcpyHostToDevice );
    cudaMemcpy( g_conf->force, conf->force, float4_array_size, cudaMemcpyHostToDevice );
}

void config_DeviceToHost (config* conf, gpu_config* g_conf, int icopy)
{
    size_t float4_array_size = sizeof(float4) * conf->natom;
    cudaMemcpy( conf->pos, g_conf->pos, float4_array_size, cudaMemcpyDeviceToHost );
    if (icopy > 0 )
    {
        cudaMemcpy( conf->vel, g_conf->vel, float4_array_size, cudaMemcpyDeviceToHost );

        if (icopy > 1)
        {
            cudaMemcpy( conf->force, g_conf->force, float4_array_size, cudaMemcpyDeviceToHost );
        }
    }
}

