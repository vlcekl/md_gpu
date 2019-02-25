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
#include "bonded_forces.h"

// Computes the forces between one particle and all others
__global__ void bonded_forces( float4* force, uint2* bndlist, unsigned int* bndlist_max) 
{
    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if( i < natoms )
    {

	float4 posi = tex1Dfetch(posTex, i);
        float4 posj;

	int num_neighbors = bndlist_max[i];

        uint2 next_j = bndlist[i];  // prefetch

	float4 f = force[i];
//        float4 f = {0.0f, 0.0f, 0.0f, 0.0f}; //make_float4(rcutsq, (float)i, (float)natoms, 4.0f);
  	for(int jx=0; jx < num_neighbors; jx++)
  	{
            unsigned int j = next_j.x;  // bonded particle
            unsigned int jt = next_j.y; // Bond type

            next_j = bndlist[i + (jx+1)*pitch_f];  // prefetch

            posj = tex1Dfetch(posTex, j);

            float kk = bnd1[jt];
            float r0 = bnd2[jt];

            // Distance for the minimum image
            posj.x = posi.x - posj.x;
            posj.y = posi.y - posj.y;
            posj.z = posi.z - posj.z;

            posj.x -= fbox.x * rintf(posj.x*fboxi.x);
            posj.y -= fbox.y * rintf(posj.y*fboxi.y);
            posj.z -= fbox.z * rintf(posj.z*fboxi.z);

            // Square distance 
            float rsq = posj.x*posj.x + posj.y*posj.y + posj.z*posj.z; // r^2

            // Forces + virial
            float force_added = kk*(r0*rsqrtf(rsq) - 1.0f); // harmonic forces
#ifdef __DEVICE_EMULATION__
//    if (i == 0 || i == 1) 
//        f.x = 0.0f;f.y = 0.0f; f.z = 0.0f;
#endif

            f.w += force_added * rsq;     // virial*6
            f.x += force_added * posj.x;  // forces
            f.y += force_added * posj.y;
            f.z += force_added * posj.z;
#ifdef __DEVICE_EMULATION__
    if (i == 0 || i == 1) {
        printf("Bond: %d %d %f %f %f %f %f %f %f\n", i, jt,f.x, f.y, f.z, force_added, kk, r0, rsqrtf(rsq));
        printf("Bond: %d %d %f\n", i, j, sqrt(rsq));
    }
#endif
  	}
        //    __syncthreads();
        force[i] = f; 
    }	// end if i < natoms
}

// Computes the forces between one particle and all others
__global__ void angle_forces( float4* force, uint4* anglist, unsigned int* anglist_max) 
{
    int id = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if( id < natoms )
    {
	float4 posid = tex1Dfetch(posTex, id);
        float4 posi;
        float4 posj;
        float4 posk;

	int num_neighbors = anglist_max[id];

        uint4 next_j = anglist[id];  // prefetch

#ifdef __DEVICE_EMULATION__
//    if (id == 0 || id == 1 || id == 2)
//        printf("ix: %d j: %d, k: %d, it %d\n", next_j.x, next_j.y, next_j.z, next_j.w);
#endif

	float4 f = force[id];
  	for(int jx=0; jx < num_neighbors; jx++)
  	{

            posj = tex1Dfetch(posTex, next_j.y);
            posi = posid;
            posk = tex1Dfetch(posTex, next_j.z);

            if (next_j.x == 1) { posi = posj; posj = posid; }

            float ff = ang1[next_j.w];
            float thet0 = ang2[next_j.w];

            if (jx < num_neighbors-1) 
                next_j = anglist[id + (jx+1)*pitch_f];  // prefetch

            // Distance for the minimum image
            float ax = posi.x - posj.x;
            float ay = posi.y - posj.y;
            float az = posi.z - posj.z;
            ax -= fbox.x * rintf(ax*fboxi.x);
            ay -= fbox.y * rintf(ay*fboxi.y);
            az -= fbox.z * rintf(az*fboxi.z);
            float ai = rsqrtf(ax*ax + ay*ay + az*az);

            float bx = posk.x - posj.x;
            float by = posk.y - posj.y;
            float bz = posk.z - posj.z;
            bx -= fbox.x * rintf(bx*fboxi.x);
            by -= fbox.y * rintf(by*fboxi.y);
            bz -= fbox.z * rintf(bz*fboxi.z);
            float bi = rsqrtf(bx*bx + by*by + bz*bz);

            // cos(theta)
            float ab = ax*bx + ay*by + az*bz;
            float cang = ab*ai*bi;
            if (cang <= -0.9999999f) cang = -0.9999999f;
            if (cang >= 0.9999999f) cang = 0.9999999f;

#ifdef __DEVICE_EMULATION__
//    if (id == 0 || id == 1 || id == 2) {
//        f.x = 0.0f;f.y = 0.0f; f.z = 0.0f;
//        printf("%d ff: %f thet0: %f, thet: %f\n", id, ff, thet0, acosf(cang));
//    }
#endif

            // Forces + virial
            ff *= (acosf(cang) - thet0)*rsqrtf(1.0f - cang*cang)*ai*bi; // harmonic forces
            //float ff = ang1[next_j.w]*(acosf(cang) - ang2[next_j.w])*rsqrtf(1.0f - cang*cang)*ai*bi; // harmonic forces

            ai = ab*ai*ai;
            posj.x = bx - ai*ax;
            posj.y = by - ai*ay;
            posj.z = bz - ai*az;

            if (next_j.x == 1) {
                bi = ab*bi*bi;
                posj.x = bi*bx - ax - posj.x;
                posj.y = bi*by - ay - posj.y;
                posj.z = bi*bz - az - posj.z;
            }

            f.x += ff * posj.x;  // forces
            f.y += ff * posj.y;
            f.z += ff * posj.z;
//            f.w += ff * rsq;     // virial*6
#ifdef __DEVICE_EMULATION__
//    if (id == 0 || id == 1 || id == 2) 
//        printf("Ang: %d %f %f %f %f\n", id, f.x, f.y, f.z, ff);
#endif
  	}
        //    __syncthreads();
        force[id] = f; // ADD not OVERWRITE !
    }	// end if i < natoms
}

