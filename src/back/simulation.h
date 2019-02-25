/*
   File name: simulation.h
   Date:      2009/04/01 02:04
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


#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <string>
#include <vector_types.h>

using namespace std;

#define LINE_LENGTH             80

#define PH_kB			1.3806503e-23		// Boltzmann constant
#define PH_Na			6.0221415e23		// Avogadro's constant
#define M_Pi			3.1415927f		// Avogadro's constant
// Erf approximation parameters
#define A1                      0.254829592f
#define A2                      -0.284496736f
#define A3                      1.421413741f
#define A4                      -1.453152027f
#define A5                      1.061405429f
#define Pex                     0.3275911f


#define NUM_PARTICLE_TYPES      8
#define NUM_BOND_TYPES      8

// config - holds a configuration (changes every time step)
typedef struct {
    int natom; // particle numbers
    int *itype;  // particle types
    float4 box; // box dimensions + some npt variable
    float4 boxi; // inverse box dimensions + some npt variable
    float4 *pos; // positions + charge
    float4 *vel; // Velocities + mass
    float4 *force; // forces
} config;

// ffield - holds force field parameters (constant over the whole run)
typedef struct {
    int ntype;
    int *nt;
    char **name;
    float *mass;
    float *charge;
    float *lj1;
    float *lj2;
    float *bnd1;
    float *bnd2;
    float *ang1;
    float *ang2;
    int tbonds;
    int tangle;
    int nrigid;
    int nbonds;
    int nangle;
    int nmaxbond;
    unsigned int *bndlist_max;
    uint2 *bndlist;
    uint4 *iatn;
    unsigned int *anglist_max;
    uint4 *anglist;
    uint4 *iang;
} ffield;

typedef struct {
    int timeq;  // equilibration timesteps
    int timrun;  // production timesteps
    int timprn;  // print interval
    float rcutsq;  // real cutoff
    float reacf;  // real cutoff
    int lkspace; // any long range coulombic contribution?
    int kmax;  // number of k-vectors
    float alfa; // Ewald's alpha
    float skin; // skin for a neighborlist
    int nmaxlist;
    float rdel; // bin size for pcfs
} params;

typedef struct { // Not sure if it fits here
    float4 tdm;  // x-Temperature, y-Pressure, z-Kinetic energy, w-Virial
    float4 energ_nonbond; // x-vdw energy, y-Coulombic short range, z-Coul. long range, w-other
    float4 energ_bonded;  // x-bonds, y-angles, z-dihedrals, w-core/shell
} thermodyn;

typedef struct {
    char* ensemble; // type of ensemble
    int enstype;
    float dt;  // timestep
    float temper_0; // required temperature in NVT
    float tcons;
    float press_0; // required temperature in NVT
    float pcons;
    float temper;
    float* Ke_sum;
    float* Wi_sum;
    float* Ke_total;
    float* Wi_total;
    thermodyn* therm;
} integr;

typedef struct {
    char* simname;  // name of the simulation (from the command line)
    params* pars;
    integr* integ;
    config* conf;
    ffield* fld;
} simulation;

simulation* SimulationInit( char *name_sim );
void SimulationRun( simulation* sim );
void SimulationFinish( simulation* sim );
ffield* FFieldInput( char *name );
void FFieldFinish( ffield* fld );
config* ConfigInput( char *name, float *mass );
void ConfigFinish( config* cfg, ffield* fld, char *name );

#endif

/* end of simulation.h */
