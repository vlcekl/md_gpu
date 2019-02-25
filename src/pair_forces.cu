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
#include "neighborlist.h"
#include "gpu_properties.h"
#include "pair_forces.h"

// Computes the forces between one particle and all others
__global__ void pair_forces( float4* force, unsigned int* nnlist, unsigned int* nnlist_max) // maybe enough
{
    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if( i < natoms )
    {
	float4 posi = tex1Dfetch(posTex, i);
        float4 posj;

	int num_neighbors = nnlist_max[i];
        int it = __float2int_rd(posi.w);
        float qi = (posi.w - __int2float_rd(it))*10.0f - 5.0f;
        it *= ntype;

	int j;
        int next_j = nnlist[i];  // prefetch

        float4 f = {0.0f, 0.0f, 0.0f, 0.0f}; //make_float4(rcutsq, (float)i, (float)natoms, 4.0f);
//        float fdx = 0.0;
//        float fdy = 0.0;
//        float fdz = 0.0;
//        float fdw = 0.0;
//        float cx = 0.0f;
//        float cy = 0.0f;
//        float cz = 0.0f;
//        float cw = 0.0f;
  	for(int jx=0; jx < num_neighbors; jx++)
  	{
            j = next_j;
            next_j = nnlist[i + (jx+1)*pitch_f];  // prefetch

            posj = tex1Dfetch(posTex, j);

            int jt = __float2int_rd(posj.w);
            float qj = (posj.w - __int2float_rd(jt))*10.0f - 5.0f;
            jt += it;
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
//            float ff = qfac*qi*qj*(ri - reacf*rsq); // Coulomb
            float ff = qfac*qi*qj*(ri - reacf*rsq); // Coulomb

            ri *= ri;            // r^-2
            ff *= ri;
            ri *= ri; ri *= ri;  // r^-8
            ff += (lja * ri * rsq - ljc) * ri; // 12-6 VdW // consider putting lj to shared memory

            if (rsq > rcutsq) { ff = 0.0f; }

#ifdef __DEVICE_EMULATION__
        if ((i == 0 && j == 4) || (i == 4 && j == 0)) printf("Pair: %d %d %e %f %f %f\n", i, j, ff, posj.x, posj.y, posj.z);
#endif
//            float yw = ff*rsq    - cw;
//            float tw = f.w + yw;
//            cw = (tw - f.w) - yw;
//            f.w = tw;     // virial*6
//            float yx = ff*posj.x - cx;
//            float tx = f.x + yx;
//            cx = (tx - f.x) - yx;
//            f.x = tx;  // forces
//            float yy = ff*posj.y - cy;
//            float ty = f.y + yy;
//            cy = (ty - f.y) - yy;
//            f.y = ty;
//            float yz = ff*posj.z - cz;
//            float tz = f.z + yz;
//            cz = (tz - f.z) - yz;
//            f.z = tz;

            f.w += ff * rsq;     // virial*6
            f.x += ff * posj.x;  // forces
            f.y += ff * posj.y;
            f.z += ff * posj.z;

//            fdw += ff * rsq;     // virial*6
//            fdx += ff * posj.x;  // forces
//            fdy += ff * posj.y;
//            fdz += ff * posj.z;
  	}
        //    __syncthreads();
        force[i] = f;
//        force[i] = make_float4(fdx, fdy, fdz, fdw);
    }	// end if i < natoms
}

__global__ void pair_forces_ewa( float4* force, unsigned int* nnlist, unsigned int* nnlist_max) // maybe enough
{
    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if( i < natoms )
    {
	float4 posi = tex1Dfetch(posTex, i);
        float4 posj;

	int num_neighbors = nnlist_max[i];
        int it = __float2int_rd(posi.w);
        float qi = (posi.w - __int2float_rd(it))*10.0f - 5.0f;
        it *= ntype;

	int j;
        int next_j = nnlist[i];  // prefetch

        float4 f = {0.0f, 0.0f, 0.0f, 0.0f}; //make_float4(rcutsq, (float)i, (float)natoms, 4.0f);
  	for(int jx=0; jx < num_neighbors; jx++)
  	{
            j = next_j;
            next_j = nnlist[i + (jx+1)*pitch_f];  // prefetch

            posj = tex1Dfetch(posTex, j);

            int jt = __float2int_rd(posj.w);
            float qj = (posj.w - __int2float_rd(jt))*10.0f - 5.0f;
            jt += it;
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

            // real Ewald  (Wolf)
            float ff = __fdividef(1.0f,1.0f + palfa*rsq*ri); 

            ff = spiialfa*rsq*ri + ((((a5*ff + a4)*ff + a3)*ff + a2)*ff + a1)*ff;
            ff *= exp2f(alfasq*rsq)*ri;
            ff -= wolf;
            ff *= qfac*qi*qj;

#ifdef __DEVICE_EMULATION__
//        if (rsq > rcutsq-0.5 && rsq < rcutsq+0.5 && jt == 0) printf("Wolf: %d %d %d %d %f %f %f\n", i, j, it/ntype, jt-it, sqrt(rsq), wolf, ff);
#endif

            // Lennard-Jones
            ri *= ri;            // r^-2
            ff *= ri;
            ri *= ri; ri *= ri;  // r^-8
            ff += (lja * ri * rsq - ljc) * ri; // 12-6 VdW // consider putting lj to shared memory

            if (rsq > rcutsq) { ff = 0.0f; }

            f.w += ff * rsq;     // virial*6
            f.x += ff * posj.x;  // forces
            f.y += ff * posj.y;
            f.z += ff * posj.z;
  	}
        force[i] = f;
    }
}

