/*
   File name: integr.h
   Date:      2009/03/31 19:54
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


#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include <vector_types.h>
#include "simulation.h"
#include "gpu_properties.h"
#include "config.h"

typedef struct { // Not sure if it fits here
    float4 tdm;  // x-Temperature, y-Pressure, z-Kinetic energy, w-Virial
    float4 energ_nonbond; // x-vdw energy, y-Coulombic short range, z-Coul. long range, w-other
    float4 energ_bonded;  // x-bonds, y-angles, z-dihedrals, w-core/shell
} gpu_thermodyn;

typedef struct {        // Thermostat & barostat parameters for Nose-Hoover thermostats
    float2 st ;         // 'position' and 'velocity' of thermostat degree of freedom
    float2 sp ;         // 'position' and 'velocity' of barostat degree of freedom
    float* Ke_sum;
    float* Wi_sum;
    float* Ke_total;
    float* Wi_total;
    gpu_thermodyn* g_therm;
} gpu_integr;

// NVE/T Leapfrog integrator with Berendsen thermostat
gpu_integr* integrate_init( integr* integ, gpu_struct* g_prop, config* conf  );
void integrate_finish(integr* integ, gpu_integr* g_integ);
void integrate( gpu_config* g_conf, gpu_integr* g_integ, gpu_struct* g_prop, integr* integ, config* conf);
void integrate_pre( gpu_config* g_conf, gpu_integr* g_integ, gpu_struct* g_prop, integr* integ, config* conf);
void integrate_cor( gpu_config* g_conf, gpu_integr* g_integ, gpu_struct* g_prop, integr* integ, config* conf);
void integrate_HostToDevice (integr* integ, gpu_integr* g_integ, gpu_struct* g_prop);
void integrate_DeviceToHost (integr* integ, gpu_integr* g_integ, gpu_struct* g_prop);
__global__ void leapfrog_nve( float4 *pos, float4 *vel, float4 *force );
__global__ void leapfrog_nvt_berendsen( float4 *pos, float4 *vel, float4 *force, float scale, float* Ke_sum );
__global__ void leapfrog_npt_berendsen( float4 *pos, float4 *vel, float4 *force, float scale, float2* Ke_sum );
//__global__ void kinetic_energy_sum( float* Ke_sum, float* Ke_total);
__global__ void kinetic_energy_sum( float* Ke_sum);
//__global__ void kinetic_energy_and_virial_sum( float2* Ke_sum, float2* Ke_total );
__global__ void kinetic_energy_and_virial_sum( float2* Ke_sum );

#endif

/* end of integr.h */
