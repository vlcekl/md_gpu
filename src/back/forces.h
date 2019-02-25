/*
   File name: pair_forces.h
   Date:      2009/04/01 01:18
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


#ifndef __PAIR_FORCES_H__
#define __PAIR_FORCES_H__

#include <vector_types.h>
#include "simulation.h"
#include "gpu_properties.h"
#include "config.h"
#include "neighborlist.h"

typedef struct {
    float something;
    uint2* bndlist;
    unsigned int* bndlist_max;
    uint4* anglist;
    unsigned int* anglist_max;
//    float2* lj;
} gpu_ffield;

// pair forces
gpu_ffield* forces_init(ffield* fld, params* pars, config* conf, gpu_struct* g_prop, gpu_config* g_conf);
void forces(gpu_struct* g_prop, gpu_nlist* g_nlist, gpu_config* g_conf, gpu_ffield* g_fld, ffield* fld, config* conf, params* pars);
void forces_finish(ffield* fld, gpu_ffield* g_fld);
void forces_HostToDevice (ffield* fld, gpu_ffield* g_fld );
void forces_DeviceToHost (ffield* fld, gpu_ffield* g_fld);

void energy(gpu_struct* g_prop, gpu_nlist* g_nlist, gpu_config* g_conf, gpu_ffield* g_fld, ffield* fld, config* conf, params* pars);

#endif

/* end of pair_forces.h */

