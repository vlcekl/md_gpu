/*
   File name: prints.cu
   Date:      2009/04/03 16:35
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

#include <stdlib.h>
#include <stdio.h>
#include <vector_types.h>

#include "gpu_properties.h"
#include "config.h"
#include "forces.h"
#include "integrator.h"
#include "prints.h"
#include "runGPU.h"
#include "prints.h"


void print_statistics( unsigned int cur_timestep, integr* integ )
{
    printf( "%i\t%f\t%f\n", cur_timestep, integ->temper, integ->Ke_sum[0]  );
//    print_random_stuff(pars, integ, fld, conf, 1)
}

void print_random_stuff( params* pars, integr* integ, ffield* fld, config* conf,
gpu_struct* g_prop, gpu_config* g_conf, gpu_integr* g_integ, gpu_ffield* g_fld, int id )
{
    switch (id) {
        case 0:
            config_DeviceToHost(conf, g_conf, 0);
            integrate_DeviceToHost(integ, g_integ, g_prop);
            forces_DeviceToHost(fld, g_fld);

            printf( "Sample input copied back from GPU\n" );
            printf( "sim: %s %d %f %f %f %f\n", integ->ensemble, pars->timrun, pars->reacf, pars->rcutsq, integ->dt, pars->rdel );
            printf( "fld: %d %d %s %f %f %f %f\n", fld->ntype, fld->nt[0], fld->name[0], fld->mass[0], fld->charge[0], fld->lj1[0], fld->lj2[0] );
            printf( "cfg: %d %d %f %f %f %f\n", conf->natom, conf->itype[0], conf->box.x, conf->pos[0].x, conf->vel[0].x, conf->pos[20].y );
            break;
        case 1:
            config_DeviceToHost(conf, g_conf, 2);
            float4 fsum;
            float4 psum;
            fsum.x = 0.0f; fsum.y = 0.0f; fsum.z = 0.0f;
            psum.x = 0.0f; psum.y = 0.0f; psum.z = 0.0f;
            for (int i=0; i < conf->natom; i++) {
               fsum.x += conf->force[i].x;
               fsum.y += conf->force[i].y;
               fsum.z += conf->force[i].z;
               psum.x += conf->vel[i].x*conf->vel[i].w;
               psum.y += conf->vel[i].y*conf->vel[i].w;
               psum.z += conf->vel[i].z*conf->vel[i].w;
            }
            printf( "force: %f %f %f\n", fsum.x, fsum.y, fsum.z);
            printf( "moment: %f %f %f\n", psum.x, psum.y, psum.z);
            //printf( "pos: %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z, conf->pos[0].w);
            //printf( "vel: %f %f %f %f\n", conf->vel[0].x, conf->vel[0].y, conf->vel[0].z, conf->vel[0].w);
            //printf( "for: %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w);
            //printf( "pos: %f %f %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z,
            //conf->pos[conf->natom-1].x, conf->pos[conf->natom-1].y, conf->pos[conf->natom-1].z );
            //printf( "for: %f %f %f %f %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z, conf->force[0].w,
            //conf->force[conf->natom-1].x, conf->force[conf->natom-1].y, conf->force[conf->natom-1].z, conf->force[conf->natom-1].w );
            break;
        case 2:
            printf( "g_prop: %d %d %d\n", g_prop->dimGrid.x, g_prop->dimBlock.x, g_prop->sharedMemSize);
        
            printf( "integr1: %f %f %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z,
            conf->force[conf->natom-1].x, conf->force[conf->natom-1].y, conf->force[conf->natom-1].z );

            config_DeviceToHost(conf, g_conf, 0);

            printf( "integr1b: %f %f %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z,
            conf->pos[conf->natom-1].x, conf->pos[conf->natom-1].y, conf->pos[conf->natom-1].z );

            printf( "integr2: %f %f %f %f %f %f\n", conf->force[0].x, conf->force[0].y, conf->force[0].z,
            conf->force[conf->natom-1].x, conf->force[conf->natom-1].y, conf->force[conf->natom-1].z );
            printf( "integr2b: %f %f %f %f %f %f\n", conf->pos[0].x, conf->pos[0].y, conf->pos[0].z,
            conf->pos[conf->natom-1].x, conf->pos[conf->natom-1].y, conf->pos[conf->natom-1].z );
            break;
        default:
            printf( "Nazdar.\n");
    }
}

void print_ave_force_momentum(gpu_config* g_conf, config* conf)
{
    config_DeviceToHost(conf, g_conf, 2);
    float4 fsum;
    float4 psum;
    fsum.x = 0.0f; fsum.y = 0.0f; fsum.z = 0.0f;
    psum.x = 0.0f; psum.y = 0.0f; psum.z = 0.0f;
    for (int i=0; i < conf->natom; i++) {
       fsum.x += conf->force[i].x;
       fsum.y += conf->force[i].y;
       fsum.z += conf->force[i].z;
       psum.x += conf->vel[i].x*conf->vel[i].w;
       psum.y += conf->vel[i].y*conf->vel[i].w;
       psum.z += conf->vel[i].z*conf->vel[i].w;
    }
    printf( "force: %f %f %f\n", fsum.x, fsum.y, fsum.z);
    printf( "moment: %f %f %f\n", psum.x, psum.y, psum.z);
}

/* end of prints.cu */
