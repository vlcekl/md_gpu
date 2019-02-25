/*
   File name: measure.cu
   Date:      2009/03/31 18:38
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

#include "energy.h"

__global__ void energy_ewa( float4* ener, unsigned int* nnlist, unsigned int* nnlist_max) // maybe enough
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
            float lja = C1_12*lj1[jt];
            float ljc = C1_6*lj2[jt];

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
            ff = ((((a5*ff + a4)*ff + a3)*ff + a2)*ff + a1)*ff;
            ff *= exp2f(alfasq*rsq);
            ff -= wolf;

            ri *= ri; ri *= ri*ri; // r^-6

#ifdef __DEVICE_EMULATION__
//        if (rsq > rcutsq-0.5 && rsq < rcutsq+0.5 && jt == 0) printf("Wolf: %d %d %d %d %f %f %f\n", i, j, it/ntype, jt-it, sqrt(rsq), wolf, ff);
#endif
            if (rsq < rcutsq) { 
                f.x = ff*qfac*qi*qj*ri; // Coulomb
                f.y += (lja * ri - ljc) * ri; // LJ
            }
  	}
        ener[i] = f;
    }
}
