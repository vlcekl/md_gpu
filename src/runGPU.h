/*
   File name: run_gpu.h
   Date:      2009/04/04 00:04
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


#ifndef __RUNGPU_H__
#define __RUNGPU_H__

#include <cuda.h>
#include "integrator.h"
#include "pair_forces.h"
#include "simulation.h"

void SimulationRun_CUDA( params* pars, integr* integ, ffield* fld, config* conf );

#ifndef CUT_CHECK_ERROR
#define CUT_CHECK_ERROR(errorMessage) do {                                 \
cudaError_t err = cudaGetLastError();                                    \
if( cudaSuccess != err) {                                                \
printf("Cuda error: %s in file '%s' in line %i : %s.\n",    \
errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
exit(EXIT_FAILURE);                                                  \
}                                                                        \
err = cudaThreadSynchronize();                                           \
if( cudaSuccess != err) {                                                \
printf( "Cuda error: %s in file '%s' in line %i : %s.\n",    \
errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
exit(EXIT_FAILURE);                                                  \
} } while (0)
#endif

#endif

/* end of run_gpu.h */
