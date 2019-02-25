/*
 *  runSimulation_kernel.cu
 *  SimpleMD
 *
 *  Created by Aaron Thompson on 6/30/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __DEVICE_EMULATION__
#include <stdio.h>
#endif
#include "gpu_properties.h"
#include "neighborlist.h"

#define WRAP(x,m) (((x)<m)?(x):(x-m))  // Mod without divide, works on values from 0 up to 2m

__constant__ float rskinsq;  // max. distance for nnlist rcut + skin
__constant__ float rmaxdif;  // maximum displacement before nnlist update
__constant__ float4 nbox;  // box dimensions
__constant__ float4 nboxi; // inverse box dimensions
__constant__ int nposit; // number of positions (atoms)
//__constant__ int nmaxlist; // maximum nnlist length
__constant__ int pitch; // nnlist array pitch
__constant__ unsigned int numTiles; // nnlist array pitch

extern __shared__ float4 sharedPos[];
extern __shared__ float sdata[];

gpu_nlist* neighborlist_init(config* conf, gpu_struct* g_prop, params* pars, gpu_config* g_conf, ffield* fld)
{
    gpu_nlist* g_nlist;
    g_nlist = (gpu_nlist *)malloc(sizeof(gpu_nlist));

    float arskinsq = pow(pars->skin + sqrt(pars->rcutsq), 2);
    float armaxdif = pars->skin*pars->skin/4.0f;
    int apitch = g_prop->dimGrid.x*g_prop->dimBlock.x;

    cudaMemcpyToSymbol( "rskinsq", &arskinsq, sizeof(float) );
    cudaMemcpyToSymbol( "rmaxdif", &armaxdif, sizeof(float) );
    cudaMemcpyToSymbol( "nposit", &conf->natom, sizeof(int) );
    cudaMemcpyToSymbol( "nbox", &conf->box, sizeof(float4) ); // may not be constant for npt
    cudaMemcpyToSymbol( "nboxi", &conf->boxi, sizeof(float4) ); // may not be constant for npt
//    cudaMemcpyToSymbol( "nmaxlist", &pars->nmaxlist, sizeof(int) ); checkCUDAError("nnlist"); 
    cudaMemcpyToSymbol( "pitch", &apitch, sizeof(int) ); checkCUDAError("nnlist"); 

    unsigned int numi = conf->natom / g_prop->block_size;
    cudaMemcpyToSymbol( "numTiles", &numi, sizeof(unsigned int) );

    //cudaMalloc( (void **)&g_nlist->update, sizeof(int) );
    cudaMalloc( (void **)&g_nlist->update, sizeof(int) );

    // initial 'old' positions
    size_t float4_array_size = sizeof(float4) * conf->natom;
    cudaMalloc( (void **)&g_nlist->pos_old, float4_array_size );
    cudaMemcpy( g_nlist->pos_old, conf->pos, float4_array_size, cudaMemcpyHostToDevice );

    size_t int_array_size = sizeof(unsigned int) * apitch;
    cudaMalloc( (void **)&g_nlist->nnlist_max, int_array_size );

    int_array_size = sizeof(unsigned int) * apitch * pars->nmaxlist;
    cudaMalloc( (void **)&g_nlist->nnlist, int_array_size );

    // Exclusions
    uint4* excl;
    size_t uint4_array_size = sizeof(uint4) * apitch;
    excl = (uint4 *)malloc( uint4_array_size );
    for (int i = 0; i < conf->natom; i++) {
        excl[i].x = conf->natom + 1;
        excl[i].y = conf->natom + 1;
        excl[i].z = conf->natom + 1;
        excl[i].w = conf->natom + 1;
    }
    for (int it = 0; it < fld->nbonds; it++) { 
        unsigned int i = fld->iatn[it].x-1;
        unsigned int j = fld->iatn[it].y-1;
        if (excl[i].x > conf->natom)
            excl[i].x = j;
        else if (excl[i].y > conf->natom) 
            excl[i].y = j;
        else if (excl[i].z > conf->natom) 
            excl[i].z = j;
        else if (excl[i].w > conf->natom) 
            excl[i].w = j;
        else 
            printf("Too many exclusions for particle %d!\n", i);
        if (excl[j].x > conf->natom)
            excl[j].x = i;
        else if (excl[i].y > conf->natom) 
            excl[j].y = i;
        else if (excl[i].z > conf->natom) 
            excl[j].z = i;
        else if (excl[i].w > conf->natom) 
            excl[j].w = i;
        else 
            printf("Too many exclusions for particle %d!\n", j);
    }
    cudaMalloc( (void **)&g_nlist->excl, uint4_array_size);checkCUDAError("nnmalloc");
    cudaMemcpy( g_nlist->excl, excl, uint4_array_size, cudaMemcpyHostToDevice);checkCUDAError("nnmemcp");

    // initialize neighborlist
    checkCUDAError("nsquared_zero"); 
    neighborlist_nsquared<<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>>(g_conf->pos, g_nlist->pos_old, g_nlist->nnlist, g_nlist->nnlist_max, g_nlist->excl);
    checkCUDAError("nsquared_start");

    return g_nlist;
}

void neighborlist(gpu_struct* g_prop, gpu_nlist* g_nlist, gpu_config* g_conf, config* conf)
{
//    size_t int_array_size = g_prop->dimGrid.x*g_prop->dimBlock.x*sizeof(unsigned int);
//    unsigned int* auxlist;
//    auxlist = (unsigned int*) malloc( int_array_size );
//    cudaMemcpy( auxlist, g_nlist->nnlist_max, int_array_size, cudaMemcpyDeviceToHost );
//    for (int i=0; i < conf->natom; i++)  printf( "nlist: %d  %d\n", i, auxlist[i]);

    int* update;
    update = (int *)malloc(sizeof(int));
    *update = 0;

//        size_t float4_array_size = sizeof(float4) * conf->natom;
//        cudaMemcpy( conf->force, g_conf->force, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "nfor1: %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w);
//        cudaMemcpy( conf->pos, g_conf->pos, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "npos1: %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z, conf->pos[0].w);
    // check if nnlist update needed
    cudaMemcpy( g_nlist->update, update, sizeof(int), cudaMemcpyHostToDevice );
    neighborlist_check<<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>>(g_conf->pos, g_nlist->pos_old, g_nlist->update);
    cudaMemcpy( update, g_nlist->update, sizeof(int), cudaMemcpyDeviceToHost );
//        cudaMemcpy( conf->force, g_conf->force, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "nfor2: %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w);
//        cudaMemcpy( conf->pos, g_conf->pos, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "npos2: %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z, conf->pos[0].w);
    
    if (*update == 1) 
    {
//        printf("\n");
        neighborlist_nsquared<<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>>(g_conf->pos, g_nlist->pos_old, g_nlist->nnlist, g_nlist->nnlist_max, g_nlist->excl);
//        printf( "nfor3: %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w);
//        cudaMemcpy( conf->pos, g_conf->pos, float4_array_size, cudaMemcpyDeviceToHost );
//        printf( "npos3: %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z, conf->pos[0].w);
    } else {
//        printf("#");
    }

}

void neighborlist_finish(gpu_nlist* g_nlist)
{
    cudaFree( g_nlist->update );
    cudaFree( g_nlist->excl );
    cudaFree( g_nlist->pos_old );
    cudaFree( g_nlist->nnlist );
    cudaFree( g_nlist->nnlist_max );
}

__global__ void neighborlist_check(float4 *pos, float4 *pos_old, int* update )
{

    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if (i < nposit)
    {
	float4 pos_n = pos[i];
	float4 pos_o = pos_old[i];

        float dx = pos_n.x - pos_o.x;
        float dy = pos_n.y - pos_o.y;
        float dz = pos_n.z - pos_o.z;
        dx -= nbox.x * rintf(dx*nboxi.x);
        dy -= nbox.y * rintf(dy*nboxi.y);
        dz -= nbox.z * rintf(dz*nboxi.z);

	if ((dx*dx + dy*dy + dz*dz) >= rmaxdif) *update = 1;
    }
}

__global__ void neighborlist_nsquared_tile(float4 *pos, float4 *pos_old, unsigned int *nnlist, unsigned int *nnlist_max, uint4* excl)
{
    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if( i < nposit )
    {
	float4 posi = pos[i];
        uint4 exi = excl[i]; // maximum 4 exclusions (usually enough)

        int ibin = 0;
        for(unsigned int tile = 0; tile < numTiles; tile++)
	{
            // Each thread loads a particle position from global memory into  shared // load a neighborlist?
            int start = WRAP(blockIdx.x + tile, gridDim.x) * blockDim.x;
            __syncthreads();
            sharedPos[threadIdx.x] = pos[start + threadIdx.x];
            __syncthreads();

  	    for(unsigned int counter=0; counter < blockDim.x; counter++)
  	    {
                float4 posj = sharedPos[counter];
//                ibin <<= 1;

                float dx = posi.x - posj.x;
                float dy = posi.y - posj.y;
                float dz = posi.z - posj.z;

                dx -= (nbox.x) * rintf(dx*nboxi.x);
                dy -= (nbox.y) * rintf(dy*nboxi.y);
                dz -= (nbox.z) * rintf(dz*nboxi.z);

                int j = start + counter; // in shared mem - should be global
#ifdef __DEVICE_EMULATION__
        if (i == 0) printf("NNSQ  i: %d j: %d, ibin: %d, rsq: %f, rskinsq: %f\n", i, j, ibin, dx*dx + dy*dy + dz*dz, rskinsq);
#endif
                if ((dx*dx + dy*dy + dz*dz < rskinsq) && (j != i) && (j != exi.x) && (j != exi.y) && (j != exi.z) && (j != exi.w)) {
#ifdef __DEVICE_EMULATION__
        if (i == 0) printf("NNSQ ACCEPTED\n", rskinsq - (dx*dx + dy*dy + dz*dz));
#endif
		    nnlist[i + pitch*ibin] = j;
                    ibin++; 
                }
            }
        }
        nnlist_max[i] = ibin;
        pos_old[i] = posi;
    }
    else
    {
        nnlist_max[i] = 0;
    }
}

__global__ void neighborlist_nsquared(float4 *pos, float4 *pos_old, unsigned int *nnlist, unsigned int *nnlist_max, uint4* excl)
{
    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    float4 posi = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    if ( i < nposit)
        posi = pos[i];
    float px = posi.x;
    float py = posi.y;
    float pz = posi.z;

    uint4 exi = make_uint4(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);
    if ( i < nposit)
        exi = excl[i];

    int ibin = 0;

    for (int start = 0; start < nposit; start += blockDim.x)
    {
	float4 posj = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	if (start + threadIdx.x < nposit)
            posj = pos[start + threadIdx.x];

        __syncthreads();
	    sdata[threadIdx.x] = posj.x;
	    sdata[threadIdx.x + blockDim.x] = posj.y;
	    sdata[threadIdx.x + 2*blockDim.x] = posj.z;
	    sdata[threadIdx.x + 3*blockDim.x] = posj.w; //< unused, but try to get compiler to fully coalesce reads
        __syncthreads();
        
	int end_offset= blockDim.x;
        end_offset = min(end_offset, nposit - start);

        if( i < nposit )
        {
	    for (int cur_offset = 0; cur_offset < end_offset; cur_offset++)
            {
                float dx = px - sdata[cur_offset];
                float dy = py - sdata[cur_offset + blockDim.x];
                float dz = pz - sdata[cur_offset + blockDim.x*2];
                dx -= (nbox.x) * rintf(dx*nboxi.x);
                dy -= (nbox.y) * rintf(dy*nboxi.y);
                dz -= (nbox.z) * rintf(dz*nboxi.z);

                int j = start + cur_offset;
#ifdef __DEVICE_EMULATION__
//        if ((dx*dx + dy*dy + dz*dz) < rskinsq) {
//            printf("NN: %d %d %d ", i, j, ibin);
//            printf("%d %d %d %d ", exi.x, exi.y, exi.z, exi.w);
//            printf("%f %f\n", sqrtf(dx*dx + dy*dy + dz*dz), sqrtf(rskinsq));
//        }
#endif
                if ((dx*dx + dy*dy + dz*dz < rskinsq) && (j != i) && (j != exi.x) && (j != exi.y) && (j != exi.z) && (j != exi.w)) {
                    nnlist[i + pitch*ibin] = j;
                    ibin++; 
                }

            }
        }
    }

    nnlist_max[i] = 0;
    if (i < nposit) {
        nnlist_max[i] = ibin;
        pos_old[i] = posi;
    }
}

