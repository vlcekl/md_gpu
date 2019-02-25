/*
   File name: runGPU.cu
   Date:      2009/03/31 21:45
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

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector_types.h>
#include <cuda.h>

#include "gpu_properties.h"
#include "config.h"
#include "pair_forces.h"
#include "integrator.h"
#include "measure.h"
#include "prints.h"
#include "runGPU.h"

#pragma mark Simulation functions


void SimulationRun_CUDA( params* pars, integr* integ, ffield* fld, config* conf )
{

    gpu_struct* g_prop;
    gpu_config* g_conf;
    gpu_ffield* g_fld;
    gpu_integr* g_integ;
    //gpu_thermodyn* g_therm;                   // now tdm with integrator, later will be separate

    g_prop = gpu_properties_init(conf);     // Device properties
    g_conf = config_init(conf); checkCUDAError("config");         // Initialize configuration
    g_fld = pair_forces_init(fld, pars, conf, g_prop); checkCUDAError("ffield");// Initialize pair forces
    g_integ = integrate_init(integ, g_prop, conf); checkCUDAError("integ");    // Initialize integrator
    //g_therm = thermodynamics_init(therm);        // Initialize thermodynamics
    //measure_init();
    
    print_random_stuff( pars, integ, fld, conf, g_prop, g_conf, g_integ, g_fld, 0 );
    
    // Main computation loop through every timestep
    printf( "Starting the main loop ...\n");
    for( unsigned int t = 0; t < pars->timrun; t++ )
    {

        // pair interactions 
        pair_forces <<< g_prop->dimGrid, g_prop->dimBlock, g_prop->sharedMemSize >>> ( g_conf->force, g_conf->pos );

        // integration
        integrate(g_conf, g_integ, g_prop, integ, conf);
    	//cudaThreadSynchronize();
    
    	if( t % pars->timprn == 0 )
        {
            print_statistics(t, integ); 
            print_random_stuff( pars, integ, fld, conf, g_prop, g_conf, g_integ, g_fld, 1 );
        }
    }
    
    // finish
    pair_forces_finish(fld, g_fld);
    config_finish(conf, g_conf);
    integrate_finish(integ, g_integ);
    //thermodynamics_finish(integ, g_integ);
    //measure_finish();
    gpu_properties_finish(g_prop);
}

