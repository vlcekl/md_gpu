/*
   File name: config.h
   Date:      2009/03/31 16:58
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


#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "simulation.h"

// Configuration

typedef struct {
//    int *itype;  // particle types (will look different than pointer to int)
//    float4 box; // box dimensions + some npt variable
//    float4 boxi; // inverse box dimensions + some npt variable
    float4* pos; // positions + charge
    float4* vel; // Velocities + mass
    float4* force; // forces
} gpu_config;

gpu_config* config_init( config* conf );
void config_finish( config* conf, gpu_config* g_conf );
void config_HostToDevice (config* conf, gpu_config* g_conf);
void config_DeviceToHost (config* conf, gpu_config* g_conf, int icopy);

#endif

/* end of config.h */
