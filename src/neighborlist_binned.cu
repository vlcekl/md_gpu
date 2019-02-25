/*
 *  runSimulation_kernel.cu
 *  SimpleMD
 *
 *  Created by Aaron Thompson on 6/30/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "neighborlist_kernel.h"
#include "global_kernel_data.h"
#include "device_functions.h"
#ifdef __DEVICE_EMULATION__
#include <stdio.h>
#endif

// update_neighborlist_nsquared
// Updates the neighbor list using a simple n^2 iteration algorithm
// Called for each particle i=0 to NUM_PARTICLES-1
__global__ void update_neighborlist_nsquared(float4 *posD, unsigned int *neighborlistD)
{
	unsigned int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if(i < NUM_PARTICLES)
	{
		float4 posi = tex1Dfetch(posDTex, i);
		float4 posj;
		unsigned int neighborlist_count = 0;

		// Loop over all other particles
                unsigned int ibin = 0;
		for(unsigned int j = 0; j < natoms && neighborlist_count < nnlist_max; j++)
		{
			posj = tex1Dfetch(posDTex, j);
                        ibin <<= 1;

  	                float dx = boxh.x - fabs(boxh.x - fabs(posi.x - posj.x));
                  	float dy = boxh.y - fabs(boxh.y - fabs(posi.y - posj.y));
                        float dz = boxh.z - fabs(boxh.z - fabs(posi.z - posj.z));

                        if (dx*dx + dy*dy + dz*dz < rnnlsq && i != j) ibin++; 
		}
		// Save position of last neighbor in this row of the neighborlist array to position zero in the row
		neighborlistD[i*NEIGHBORLIST_CAPACITY] = neighborlist_count;
	}
}

//neighborlistD[i*nnlist_max + neighborlist_count + 1] = j;
//neighborlist_count++;

__global__ void update_neighborlist_binned(float4 *posD, unsigned int *neighborlistD)
{
	unsigned int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if(i < NUM_PARTICLES)
	{
		float4 pos1 = posD[i];
		uint3 bin_index1 = BIN_INDEX(pos1);
		unsigned int bin_array_index1 = bin_array_index

		float4 pos2;
		uint3 bin_index2;
		unsigned int bin_array_index2;
		unsigned int posD2Index;
		unsigned int num_particles_bin2;

		unsigned int neighborlist_count = 0;

		// Loop over all 27 bins containing & surrounding the particle i
		for(unsigned int bin2i=0; bin2i<1; bin2i++)
		{
			// Determine the bin index based on bin_array_index and bin2i
			bin_index2 = {	bin2i
			bin_array_index2 = bin_array_index;

			// Get the # of particles in this bin
			num_particles_bin2 = tex1Dfetch(binsDTex, bin_array_index2);

			// Loop over all particles in the bin
			for(unsigned int j=0; j<num_particles_bin2 && neighborlist_count<NEIGHBORLIST_CAPACITY; j++)
			{
				// Fetch the posD index of the jth particle in this bin
				posD2Index = tex1Dfetch(binsDTex, bin_array_index2 + j + 1);
				pos2 = posD[posD2Index];
				
				if(are_neighbors(&pos1, &pos2))
				{
					// The particles are within range of each other
					neighborlistD[i*NEIGHBORLIST_CAPACITY + neighborlist_count+1] = j;
					neighborlist_count++;
//						printf("Particles %i and %i with pos (%e %e %e) and (%e %e %e) are neighbors\n",
//							i, j, pos1.x, pos1.y, pos1.z, pos2.x, pos2.y, pos2.z);
				}
			} // End loop over all particles in the bin
		} // End loop over all particles in the current bin
		// Save position of last neighbor in this row of the neighborlist array to position zero in the row
		neighborlistD[i*NEIGHBORLIST_CAPACITY] = neighborlist_count;
	}
}

__global__ void neighborlist_bin_sort( float4 *posD, unsigned int *binsD )
{
	unsigned int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if( i < NUM_PARTICLES )
	{
		float4 pos = posD[i];

		unsigned int bin_index = bin_index_linear(&pos);
		// Add this particle to the bin, making sure to avoid a race condition by using an atomic operation
		atomicAdd(&binsD[bin_index], 1);
		unsigned int bin_particle_count = binsD[bin_index];
		// Don't overfill the bin
		if( bin_particle_count < BIN_CAPACITY ) {
			binsD[bin_index] = i;
		} else {
		}
	}
}
