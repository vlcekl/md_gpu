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

#include "measure.h"

// CPU version
void pair_forces( float4 *pos, measure* meas, float2 *Ue_sum )
{
    float lxh = 0.5f * box.x;
    float lyh = 0.5f * box.y;
    float lzh = 0.5f * box.z;

    float ene_vdw = 0.0f;
    float ene_coul = 0.0f;

    for (unsigned int i = 0; i < natom ; i++) {
        // go through the neighbor list
        for (unsigned int k = 0; k < nnlistmax[i]; k++) {
            unsigned int j = nnlist[i][k]

            // minimum image distance
  	    float dx = lxh - fabs(lxh - fabs(pos[i].x - pos[j].x));
  	    float dy = lyh - fabs(lyh - fabs(pos[i].y - pos[j].y));
  	    float dz = lzh - fabs(lzh - fabs(pos[i].z - pos[j].z));

            float rsq = dx*dx + dy*dy + dz*dz;

            if (rsq < rcutsq)
            {
  	        float r = sqrt(rsq);

                // Pair distribution functions
                if (lpdf) 
                {
                    int l = (int)(r*dri);
                }

                r = 1.0f/r;
                // Coulombic energies 
                ene_coul[ijt] += qfac*pos1.w*pos2.w*(r + rsq*reacf);

                // VdW energies
                ene_vdw[ijt] += (lj[jt].x * ri * rsq - lj[jt].y) * ri; // 12-6 VdW // consider putting lj to shared memory
            }
        }
    }

    Ue_sum.x = ene_coul;
    Ue_sum.y = ene_vdw;
}

