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


// BIN_INDEX
// Computes the bin index that a particle at position pos should be located in
#define BIN_INDEX(pos) \
	make_uint3( \
		(pos.x + BOX_LENGTH.x/2.0f) / BIN_LENGTH.x, \
		(pos.y + BOX_LENGTH.y/2.0f) / BIN_LENGTH.y, \
		(pos.z + BOX_LENGTH.z/2.0f) / BIN_LENGTH.z \
	);


// update_neighborlist_nsquared
// Updates the neighbor list using a simple n^2 iteration algorithm
// Called for each particle i=0 to NUM_PARTICLES-1
__global__ void update_neighborlist_nsquared(float4 *posD, unsigned int *neighborlistD)
{
	unsigned int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if(i < NUM_PARTICLES)
	{
		float4 pos1 = tex1Dfetch(posDTex, i);
		float4 pos2;
		unsigned int neighborlist_count = 0;
		
		// Loop over all other particles
		for(unsigned int j=0; j<NUM_PARTICLES && neighborlist_count<NEIGHBORLIST_CAPACITY; j++)
		{
			pos2 = tex1Dfetch(posDTex, j);

			if( are_neighbors(&pos1, &pos2) && (i != j) ) {
				neighborlistD[i*NEIGHBORLIST_CAPACITY + neighborlist_count+1] = j;
				neighborlist_count++;

			}
		}
		// Save position of last neighbor in this row of the neighborlist array to position zero in the row
		neighborlistD[i*NEIGHBORLIST_CAPACITY] = neighborlist_count;
	}
}


/*
// update_neighborlist_binned
// Updates the neighbor list using our sorted bins (should already be sorted by the time this is called)
// Called for each bin x=0 to NUM_BINS.x-1, y=0 to NUM_BINS.y-1, z=0 to NUM_BINS.z-1
__global__ void update_neighborlist_binned(float4 *posD, unsigned int *neighborlistD)
{
	uint4 bin_index = make_uint4( __mul24(blockIdx.x, blockDim.x) + threadIdx.x,
								__mul24(blockIdx.y, blockDim.y) + threadIdx.y,
								__mul24(blockIdx.z, blockDim.z) + threadIdx.z,
								0 );

	if(bin_index.x<NUM_BINS.x && bin_index.y<NUM_BINS.y && bin_index.z<NUM_BINS.z)
	{
		unsigned int bin_array_index = bin_array_index_linear(bin_index);
		unsigned int num_particles_bin1 = tex1Dfetch(binsDTex, bin_array_index);
		unsigned int num_particles_bin2;
		float4 pos1;
		float4 pos2;
		uint3 bin_index2;
		unsigned int neighborlist_count = 0;
		
		// Loop over all particles in the current bin
		for(unsigned int i=0; i<num_particles_bin1; i++)
		{
			pos1 = posD[i];

			// Loop over the current and 26 surrounding bins
			for(unsigned int 
			
			// Loop over all particles in the bin
			for(unsigned int j=0; j<num_particles_bin1 && neighborlist_count<NEIGHBORLIST_CAPACITY; j++)
			{
				pos2 = posD[j];

				bin_index2 = BIN_INDEX(pos2);
				
				// Only consider particles in the boxes at or one greater than the current index in any direction
				// This means we consider only the current box, one up, one right, and/or one deep, a 2x2x2 cube
				// (Making sure to wrap around)
				if(abs(bin_index1.x % (NUM_BINS.x-1) - bin_index2.x % (NUM_BINS.x-1)) < 2 && (i != j)) {
					if(abs(bin_index1.y % (NUM_BINS.y-1) - bin_index2.y % (NUM_BINS.y-1)) < 2) {
						if(abs(bin_index1.z % (NUM_BINS.z-1) - bin_index2.z % (NUM_BINS.z-1)) < 2) {
							// The particles are within range of each other
							neighborlistD[i*NEIGHBORLIST_CAPACITY + neighborlist_count+1] = j;
							neighborlist_count++;
	//						printf("Particles %i and %i with pos (%e %e %e) and (%e %e %e) are neighbors\n",
	//							i, j, pos1.x, pos1.y, pos1.z, pos2.x, pos2.y, pos2.z);
						}
					}
				}
			} // End loop over all particles in the bin
		} // End loop over all particles in the current bin

		// Save position of last neighbor in this row of the neighborlist array to position zero in the row
		neighborlistD[i*NEIGHBORLIST_CAPACITY] = neighborlist_count;
	}
}
*/


// update_neighborlist_binned
// Updates the neighbor list using our sorted bins (should already be sorted by the time this is called)
// Called for each particle i=0 to NUM_PARTICLES-1
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


// empty_neighborlist_bins
// Empties out the bins before re-sorting
// Called for each bin i=0 to NUM_BINS-1
__global__ void empty_neighborlist_bins( unsigned int *binsD )
{
	unsigned int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if( i < NUM_BINS.x )
	{
		// The first entry of a cell of binsD contains the number of particle indices in that bin
		binsD[i*NUM_BINS.x] = 0;
	}
}


// neighborlist_bin_sort
// Sorts particles into the different bins
// Called for each particle i=0 to NUM_PARTICLES-1
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
//			printf("Putting particle %i with pos (%e, %e, %e) in bin at index = %i\n",
//				i, pos.x, pos.y, pos.z, bin_index);
			binsD[bin_index] = i;
		} else {
//			printf("Bin %i for particle %i at pos (%e, %e, %e) is overfilled! Count = %i\n",
//				bin_index, i, pos.x, pos.y, pos.z, bin_particle_count);
		}
	}
}


// are_neighbors
// Determines the force and mu value between two particles [52 FLOPS]
__device__ int are_neighbors(float4 *pos1, float4 *pos2)
{
	float3 dist;
	// Calculate distance [3 FLOPS]
	dist.x = pos2->x - pos1->x;
	dist.y = pos2->y - pos1->y;
	dist.z = pos2->z - pos1->z;

	// Factor minimum image conversion into distance calculation [12 FLOPS]
	dist.x -= BOX_LENGTH.x * rintf(dist.x*BOX_LENGTH_INV.x);
	dist.y -= BOX_LENGTH.y * rintf(dist.y*BOX_LENGTH_INV.y);
	dist.z -= BOX_LENGTH.z * rintf(dist.z*BOX_LENGTH_INV.z);

	// Calculate scalar distance [6 FLOPS]
	float dist_squared = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
	// Added hack because it seems to be calculating distances between a particle and itself sometimes
	return (dist_squared < FORCE_DISTANCE_CUTOFF*FORCE_DISTANCE_CUTOFF) ? 1 : 0;
}


// bin_array_index_linear
// Converts a 3D index into a 1D index in order to linearize the array
__device__ unsigned int bin_array_index_linear(uint3 bin_index)
{
	return  bin_index.x * NUM_BINS.y * NUM_BINS.x * BIN_CAPACITY +
						 bin_index.y * NUM_BINS.x * BIN_CAPACITY +
									  bin_index.z * BIN_CAPACITY;
}


// bin_index_linear
// Computes the bin index that a particle at position pos should be located in
// Returns the result as a position in a 1D array as opposed to 3D
#define INDEX_X(pos) (pos->x + BOX_LENGTH.x/2.0f) / BIN_LENGTH.x
#define INDEX_Y(pos) (pos->y + BOX_LENGTH.y/2.0f) / BIN_LENGTH.y
#define INDEX_Z(pos) (pos->z + BOX_LENGTH.z/2.0f) / BIN_LENGTH.z
__device__ unsigned int bin_index_linear(float4 *pos)
{
	return	INDEX_X(pos) * NUM_BINS.y * NUM_BINS.x * BIN_CAPACITY +
			INDEX_Y(pos) * NUM_BINS.x * BIN_CAPACITY +
			INDEX_Z(pos) * BIN_CAPACITY;
}
