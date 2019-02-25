/*
   File name: neighborlist.h
   Date:      2009/05/21 00:32
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

#include "gpu_properties.h"

#ifndef __NEIGHBORLIST_H__
#define __NEIGHBORLIST_H__

typedef struct {    
//    float rskinsq;
//    float rmaxdif;
//    float nmaxlist
    float4* pos_old;
    unsigned int* nnlist;
    unsigned int* nnlist_max;
    uint4* excl;
    int* update;
} gpu_nlist;

gpu_nlist* neighborlist_init(config* conf, gpu_struct* g_prop, params* pars, gpu_config* g_conf, ffield* fld);
void neighborlist(gpu_struct* g_prop, gpu_nlist* g_nlist, gpu_config* g_conf, config* conf);
void neighborlist_finish(gpu_nlist* g_nlist);
__global__ void neighborlist_check(float4* pos, float4* pos_old, int* update);
__global__ void neighborlist_nsquared(float4 *pos, float4 *pos_old, unsigned int *nnlist, unsigned int *nnlist_max, uint4 *excl);

#endif

/* end of neighborlist.h */
