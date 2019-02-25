/*
   File name: npt_integ.cu
   Date:      2009/04/03 13:46
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

#include "npt_leapfrog_ber.h"

extern __shared__ float2 KW_shared[];

__global__ void leapfrog_npt_berendsen( float4 *pos, float4 *vel, float4 *force, float scale, float2* Ke_sum )
{
        
        // Leapfrog integration + particle kinetic energy

	int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

//        float Ke_i;
//        float Wi_i;
	if( i < natom )
	{       // LEAPFROG VERLET + Beredsen thermostat
            float4 veli = vel[i];
            float massif = __fdividef(dt, veli.w);
            //float massif = dt/veli.w;
            float4 f = force[i];

            veli.x = (veli.x + f.x*massif)*scale;
            veli.y = (veli.y + f.y*massif)*scale;
            veli.z = (veli.z + f.z*massif)*scale;
            vel[i] = veli;

//            Ke_i = veli.w*(veli.x*veli.x + veli.y*veli.y + veli.z*veli.z);
//            Wi_i = f.w;
        KW_shared[threadIdx.x].x = veli.w*(veli.x*veli.x + veli.y*veli.y + veli.z*veli.z);
        KW_shared[threadIdx.x].y = f.w;

            float4 posi = pos[i];
            posi.x += veli.x*dt;
            posi.x -= box.x *truncf(posi.x*2.0f*boxi.x );
            posi.y += veli.y*dt;
            posi.y -= box.y *truncf(posi.y*2.0f*boxi.y);
            posi.z += veli.z*dt;
            posi.z -= box.z *truncf(posi.z*2.0f*boxi.z );
            pos[i] = posi;
	}
        else
        {
            KW_shared[threadIdx.x].x = 0.0f;
            KW_shared[threadIdx.x].y = 0.0f;
//            Ke_i = 0.0f;
//            Wi_i = 0.0f;
        }

        // b) sum (reduction) according to CUDA documentation

//        KW_shared[threadIdx.x].x = Ke_i;
//        KW_shared[threadIdx.x].y = Wi_i;

        __syncthreads();

        for (unsigned int s = blockDim.x>>1; s > 0; s>>=1) {
            if(threadIdx.x < s) {
                KW_shared[threadIdx.x].x += KW_shared[threadIdx.x + s].x;
                KW_shared[threadIdx.x].y += KW_shared[threadIdx.x + s].y;
            }
            __syncthreads();
        }

        if (threadIdx.x == 0)
        {
            Ke_sum[blockIdx.x] = KW_shared[0];
        }
}

// Sum kinetic energy and virial in a single block according to CUDA doc + HOOMD
// If less efficient than pure energy sum, reuse the energy kernel for virial
__global__ void kinetic_energy_and_virial_sum( float2* Ke_sum)
{
        float2 ksum = {0.0f, 0.0f};
        for (unsigned int ib = 0; ib < nblocks; ib += blockDim.x )
        {
            __syncthreads();
            if (ib + threadIdx.x < nblocks)
                KW_shared[threadIdx.x] = Ke_sum[threadIdx.x + ib];
            else
                KW_shared[threadIdx.x] = make_float2(0.0f, 0.0f);
            __syncthreads();
            for (unsigned int s = blockDim.x>>1; s > 0; s>>=1)
            {
                if(threadIdx.x < s)
                {
                    KW_shared[threadIdx.x].x += KW_shared[threadIdx.x + s].x;
                    KW_shared[threadIdx.x].y += KW_shared[threadIdx.x + s].y;
                }
                __syncthreads();
            }
 
            ksum.x += KW_shared[0].x;
            ksum.y += KW_shared[0].y;
        }

        if (threadIdx.x == 0) {
            ksum.y *= float(1.0/6.0);
            Ke_sum[0] = ksum;
        }
}

/* end of npt_integ.cu */
