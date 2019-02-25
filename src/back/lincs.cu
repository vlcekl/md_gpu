/*
   File name: lincs.cu
   Date:      2009/04/07 17:48
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

#include "lincs.h"

// Linear Constraint Solver

lincs(x, xp, invmass, K, nrec, atom1, atom2, length, ncc, cmax, con, Sdiag, coef)
{
    float B[K, 3];
    float A[K, cmax]; 
    float rhs[2, K];
    float sol[K];
    int k, n, a1, a2;
    float len, p;

    int i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if (i < natom) 
    {
        // get bonded atom - j

        float dx = posi.x - posj.x;
        float dy = posi.y - posj.y;
        float dz = posi.z - posj.z;
        
        dx -= fbox.x * rintf(dx*fboxi.x);
        dy -= fbox.y * rintf(dy*fboxi.y);
        dz -= fbox.z * rintf(dz*fboxi.z);
 
        float rsq = dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i]; // r^2
        float ri = rsqrt(rsq);
 
        dx *= ri;
        dy *= ri;
        dz *= ri;

        for (n = 0; n < ncc; n++)
        { 
            k = con[i, n]; 
            A[i, n] = coef[i, n]*(dx[i]*dx[k] + dy[i]*dy[k] + dz[i]*dz[k]);
        }
        a1 = atom1[i]; 
        a2 = atom2[i];
        rhs[1, i] = Sdiag[i]*(B[i, 1]*(xp[a1,1] - xp[a2,1]) +
                              B[i, 2]*(xp[a1,2] - xp[a2,2]) +
                              B[i, 3]*(xp[a1,3] - xp[a2,3]) - length[i]);
        sol[i] = rhs[1, i]; 
    }

    solve(xp, invmass, K, nrec, atom1, atom2, ncc, con, Sdiag, B, rhs, sol );

    if (i < natom) 
    {
        a1 = atom1[i]; 
        a2 = atom2[i]; 
        p = sqrt(2*length[i]*length[i] - sqr(xp[a1,1]-xp[a2,1]) -
                                         sqr(xp[a1,2]-xp[a2,2]) -
                                         sqr(xp[a1,3]-xp[a2,3]));

        rhs[1, i] = Sdiag[i]*(length[i] - p);
        sol[i] = rhs[1, i];

        solve(xp, invmass, K, nrec, atom1, atom2, ncc, con, Sdiag, B, rhs, sol );
    }
} 

solve(xp, invmass, K, nrec, atom1, atom2, ncc, con, Sdiag, B, rhs, sol)
{
    int i, n, rec, w, a1, a2;
    w = 2;
    for (rec = 0; rec < nrec; rec++)
    {
        for (i = 0; i <  K; i++) 
        {
            rhs[w,i] = 0;
            for (n = 0; n < ncc[i]; n++)
                rhs[w,i] += A[i, con[i,n]]*rhs[3 - w, con[i, n]]
            sol[i] += rhs[w, i];
        }
        w = 3 - w;
    }
    for (i = 0; i <  K; i++)
    {
        a1 = atom1[i];
        a2 = atom2[i];
        xp[a1,1] -= invmass[a1]*B[i,1]*Sdiag[i]*sol[i];
        xp[a2,1] -= invmass[a2]*B[i,1]*Sdiag[i]*sol[i];
        xp[a1,2] -= invmass[a1]*B[i,2]*Sdiag[i]*sol[i];
        xp[a2,2] -= invmass[a2]*B[i,2]*Sdiag[i]*sol[i];
        xp[a1,3] -= invmass[a1]*B[i,3]*Sdiag[i]*sol[i];
        xp[a2,3] -= invmass[a2]*B[i,3]*Sdiag[i]*sol[i];
    }
}

/* end of lincs.cu */
