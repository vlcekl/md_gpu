/*
   File name: nvt_leapfrog_ber.cu
   Date:      2009/04/03 13:48
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

#include "nvt_leapfrog_ber.h"

__global__ void leapfrog_nvt_berendsen( float4 *pos, float4 *vel, float4 *force, float scale, float* Ke_sum )
{
        extern __shared__ float Ke_shared[];
	int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if( i < natom )
	{
            float4 veli = vel[i];
            float massif = __fdividef(dt, veli.w);
            float4 f = force[i];

            veli.x = (veli.x + f.x*massif)*scale;
            veli.y = (veli.y + f.y*massif)*scale;
            veli.z = (veli.z + f.z*massif)*scale;
            vel[i] = veli;

            Ke_shared[threadIdx.x] = veli.w*(veli.x*veli.x + veli.y*veli.y + veli.z*veli.z);

            float4 posi = pos[i];
            posi.x += veli.x*dt;
            posi.y += veli.y*dt;
            posi.z += veli.z*dt;
            posi.x -= box.x*truncf(posi.x*2.0f*boxi.x);
            posi.y -= box.y*truncf(posi.y*2.0f*boxi.y);
            posi.z -= box.z*truncf(posi.z*2.0f*boxi.z);
            pos[i] = posi;
	}
        else
        {
            Ke_shared[threadIdx.x] = 0.0f;
        }

        __syncthreads();

        // b) sum (reduction) according to CUDA documentation
        for (unsigned int s = blockDim.x>>1; s > 0; s>>=1) {
            if(threadIdx.x < s) {
                Ke_shared[threadIdx.x] += Ke_shared[threadIdx.x + s];
            }
            __syncthreads();
        }

        if (threadIdx.x == 0) { Ke_sum[blockIdx.x] = Ke_shared[0]; }
}

// Sum kinetic energy in a single block according to CUDA doc + HOOMD
__global__ void kinetic_energy_sum( float* Ke_sum )
{
        extern __shared__ float Ke_shared[];
        float ksum = 0.0f;
        for (unsigned int ib = 0; ib < nblocks; ib += blockDim.x )
        {
            __syncthreads();
            if (ib + threadIdx.x < nblocks)
                Ke_shared[threadIdx.x] = Ke_sum[threadIdx.x + ib];
            else
                Ke_shared[threadIdx.x] = 0.0f;
            __syncthreads();

            for (unsigned int s = blockDim.x>>1; s > 0; s>>=1)
            {
                if(threadIdx.x < s)
                {
                    Ke_shared[threadIdx.x] += Ke_shared[threadIdx.x + s];
                }
                __syncthreads();
            }
 
            ksum += Ke_shared[0];
        }

        if (threadIdx.x == 0) { Ke_sum[0] = ksum; }
}

