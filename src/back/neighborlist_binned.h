/*
 *  neighborlist_kernel.h
 *  SimpleMD
 *
 *  Created by Aaron Thompson on 6/30/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _neighborlist_kernel_h_
#define _neighborlist_kernel_h_

#include "sharedSimulationItems.h"

// update_neighborlist_nsquared
// Updates the neighbor list using a simple n^2 iteration algorithm
// Called for each particle pair i=0 to NUM_PARTICLES-1, j=i to NUM_PARTICLES
__global__ void update_neighborlist_nsquared(float4 *posD, unsigned int *neighborlistD);

// update_neighborlist_binned
// Updates the neighbor list using our sorted bins
// Called for each particle pair i=0 to NUM_PARTICLES-1, j=i to NUM_PARTICLES
__global__ void update_neighborlist_binned(float4 *posD, unsigned int *neighborlistD);

// empty_neighborlist_bins
// Empties out the bins before re-sorting
// Called for each bin i=0 to NUM_BINS-1
__global__ void empty_neighborlist_bins(unsigned int *binsD);

// neighborlist_bin_sort
// Sorts particles into the different bins
// Called for each particle i=0 to NUM_PARTICLES-1
__global__ void neighborlist_bin_sort(float4 *posD, unsigned int *binsD);

// are_neighbors
// Tells if particles at 2 positions are neighbors
__device__ int are_neighbors(float4 *pos1, float4 *pos2);

// bin_array_index_linear
// Converts a 4D index into a 1D index in order to linearize the array
__device__ unsigned int bin_array_index_linear(uint3 bin_index);

// bin_index_linear
// Computes the bin index that a particle at position pos should be located in
// Returns the result as a position in a 1D array as opposed to 3D
__device__ unsigned int bin_index_linear(float4 *pos);

#endif
