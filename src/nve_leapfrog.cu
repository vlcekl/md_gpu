/*
   File name: nve_leapfrog.cu
   Date:      2009/04/03 13:49
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

#include "nve_leapfrog.h"

// Called for each particle i=0 to natom
__global__ void leapfrog_nve( float4 *pos, float4 *vel, float4 *force )
{
	int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if( i < natom )
	{       // LEAPFROG VERLET //
            float4 veli = vel[i];
            float massif = __fdividef(dt, veli.w);
            //float massif = dt/veli.w;
            float4 f = force[i];
#ifdef __DEVICE_EMULATION__
        printf("Massif: %e %e %e %e %e %e %e %e\n", massif, f.x, veli.x, veli.y, veli.z, f.x*massif, f.y*massif, f.z*massif);
#endif

            veli.x += f.x*massif;
            veli.y += f.y*massif;
            veli.z += f.z*massif;
            vel[i] = veli;

            float4 posi = pos[i];
            posi.x += veli.x*dt;
            posi.x -= box.x *truncf(posi.x*2.0f*boxi.x );
            posi.y += veli.y*dt;
            posi.y -= box.y *truncf(posi.y*2.0f*boxi.y);
            posi.z += veli.z*dt;
            posi.z -= box.z *truncf(posi.z*2.0f*boxi.z );
            pos[i] = posi;
	}
}



/* end of nve_leapfrog.cu */
