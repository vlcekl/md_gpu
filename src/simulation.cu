/*
   File name: simulation.cpp
   Date:      2009/04/01 02:05
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

// Very ugly, should be reorganized

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector_types.h>

#include "runGPU.h"
#include "simulation.h"

using namespace std;

simulation *SimulationInit(char *name_sim)
{
    FILE *fin;
    char *saux;

    simulation* sim;
    sim = (simulation *)malloc(sizeof(simulation));

    sim->simname = (char *)malloc(LINE_LENGTH);
    strcpy(sim->simname, name_sim);

    
    params* pars;
    pars = (params *)malloc(sizeof(params));

    integr* integ;
    integ = (integr *)malloc(sizeof(integr));

    // read data from .inp
    fin = fopen(strcat(name_sim, ".inp") , "r");
        saux = (char *) malloc(LINE_LENGTH);

        // ensemble 
        integ->ensemble = (char *)malloc(3);
        integ->enstype = 0;
        fscanf (fin, "%s %s", saux, integ->ensemble);
        if (strcmp(integ->ensemble, "nvt") == 0)
        {
            integ->enstype = 1;
            fscanf (fin, "%s %f", saux, &integ->temper_0);
        }
        if (strcmp(integ->ensemble, "npt") == 0)
        {
            integ->enstype = 2;
            fscanf (fin, "%s %f %f", saux, &integ->temper_0, &integ->press_0);
        }

        // time and time steps
        fscanf (fin, "%s %d %d %d %f", saux, &pars->timeq, &pars->timrun, &pars->timprn, &integ->dt);
        fscanf (fin, "%s %f %d", saux, &pars->rcutsq, &pars->lkspace);

        // force & energy calcuation parameters
        if (pars->lkspace == 0) // 0 ... reaction field
        {
            float eps;
            fscanf (fin,"%s %f", saux, &eps);
            pars->reacf = 2.0f*(eps - 1.0f)/(2.0f*eps + 1.0f)/pow(pars->rcutsq, 3);
        }
        else if (pars->lkspace == 1) { // 1 ... Wolf
            fscanf (fin,"%s %f", saux, &pars->alfa);
        } 
        else //2 ... Ewald
        {
            fscanf (fin,"%s %d", saux, &pars->kmax);
            pars->alfa = 3.09/pars->rcutsq; // one way to estimate good alpha (?)
        }
        pars->rcutsq *= pars->rcutsq; // square of cutoff distance

        // neighbor list skin
        fscanf (fin, "%s %f %d", saux, &pars->skin, &pars->nmaxlist);

        // PDF bin size
        fscanf (fin, "%s %f", saux, &pars->rdel);

        free(saux);
    fclose(fin);

    // read data from .fld
    strcpy(name_sim, sim->simname);
    ffield* fld;
    fld = (ffield *)malloc(sizeof(ffield));
    fld = FFieldInput(name_sim);

    // read data from .cfg
    strcpy(name_sim, sim->simname);
    config* conf;
    conf = (config *)malloc(sizeof(config));
    conf = ConfigInput(name_sim, fld->mass);

    sim->pars = pars;
    sim->integ = integ;
    sim->fld = fld;
    sim->conf = conf;

    return(sim);
}

void SimulationRun(simulation* sim)
{
    time_t start, end;
    printf( "About to run simulation on GPU ...\n" );
    time(&start);
    SimulationRun_CUDA(sim->pars, sim->integ, sim->fld, sim->conf);
    time(&end);
    printf( "Finished simulation.\n Total time: %f s.\n", difftime(end, start));
}

void SimulationFinish(simulation* sim)
{
    FFieldFinish(sim->fld);
    ConfigFinish(sim->conf, sim->fld, sim->simname);
}

// ------ Force field ------
ffield* FFieldInput(char *name)
{
    int i, j, iaux, iatom;
    float faux1, faux2;
    char *saux;
    FILE *fin;
    ffield* fld;

    fld = (ffield *)malloc(sizeof(ffield));
    fin = fopen(strcat(name, ".fld") , "r");
        saux = (char *) malloc(LINE_LENGTH);
        fscanf (fin, "%s %d", saux, &fld->ntype); // number of atom types
        int nty =  NUM_BOND_TYPES;
        // allocate atom type and ffield structures
        saux = (char *) malloc(LINE_LENGTH);
        fld->nt = (int *)malloc( nty * sizeof(int) );
        fld->name = (char **)malloc( nty * sizeof(char *) );
        for (int i = 0; i < fld->ntype; i++) fld->name[i] = (char *) malloc(2);
        fld->mass = (float *)malloc( nty * sizeof(float) );
        fld->charge = (float *)malloc( nty * sizeof(float) );
        fld->lj1 = (float *) malloc ( nty*nty * sizeof(float) );
        fld->lj2 = (float *) malloc ( nty*nty * sizeof(float) );
        for (int i = 0; i < nty*nty; i++)
        {
            fld->lj1[i] = 0.0f;
            fld->lj2[i] = 0.0f;
        }
 
        // atomic parameters: number, name, mass, charge
        iatom = 0;
        for (int it = 0; it < fld->ntype; it++) {
            fscanf (fin, "%d %d", &iaux, &fld->nt[it]);
            fscanf (fin, "%s %f %f", fld->name[it], &fld->mass[it], &fld->charge[it]);
            fld->mass[it] *= 10000.0; // g/mol*(A/fs)**2 -> kJ/mol*(m/s)
            iatom += fld->nt[it];
        }
 
        // VdW
        int nvdw;
        fscanf (fin, "%s %d", saux, &nvdw); 
        for (int it = 0; it < nvdw; it++) {
            fscanf (fin, "%d %d %s %f %f", &i, &j, saux, &faux1, &faux2);
            i--; j--; 
            if (strcmp(saux, "LJ") == 0) // LJ form
            {
           //     faux2 *= PH_kB*PH_Na/1000.0f; // K->kJ/mol if needed
                // for forces calculation - 12 and 6 factors included
                fld->lj1[i*fld->ntype + j] = 12.0f*4.0f*faux2*pow(faux1, 12);
                fld->lj1[j*fld->ntype + i] = fld->lj1[i*fld->ntype + j];

                fld->lj2[i*fld->ntype + j] = 6.0f*4.0f*faux2*pow(faux1, 6);
                fld->lj2[j*fld->ntype + i] = fld->lj2[i*fld->ntype + j];
            }
            else if (strcmp(saux, "12-6") == 0)  // 12-6 form
            {
                // for forces calculation - 12 and 6 factors included
                fld->lj1[i*fld->ntype + j] = 12.0f*faux1;
                fld->lj1[j*fld->ntype + i] = fld->lj1[i*fld->ntype + j];
                fld->lj2[i*fld->ntype + j] = 6.0f*faux2;
                fld->lj2[j*fld->ntype + i] = fld->lj2[i*fld->ntype + j];
            }
        }

        // Fixed atoms
//        fscanf (fin, "%s %d", saux, &fld->nfixed); 
//        for (int it = 0; it < fld->nfixed; it++) { ; }

        // Rigid molecules (just water)
//        fscanf (fin, "%s %d", saux, &fld->nrigid); 
//        for (int it = 0; it < fld->nrigid; it++) { ; }

        // Bond types
        fscanf (fin, "%s %d", saux, &fld->tbonds);  // bond types
        int nbo = NUM_BOND_TYPES;
        fld->bnd1 = (float *)malloc( nbo * sizeof(float) );
        fld->bnd2 = (float *)malloc( nbo * sizeof(float) );
        for (int i = 0; i < nbo*nbo; i++)
        {
            fld->bnd1[i] = 0.0f;
            fld->bnd2[i] = 0.0f;
        }
        for (int it = 0; it < fld->tbonds; it++) { 
            fscanf (fin, "%d %f %f", &iaux, &fld->bnd1[it], &fld->bnd2[it]);
            fld->bnd1[it] *= 2.0f; // 
        }

        // Angle types
        fscanf (fin, "%s %d", saux, &fld->tangle);  // bond types
        int nag = NUM_BOND_TYPES;
        fld->ang1 = (float *)malloc( nag * sizeof(float) );
        fld->ang2 = (float *)malloc( nag * sizeof(float) );
        for (int i = 0; i < nag*nag; i++)
        {
            fld->ang1[i] = 0.0f;
            fld->ang2[i] = 0.0f;
        }
        for (int it = 0; it < fld->tangle; it++) { 
            fscanf (fin, "%d %f %f", &iaux, &fld->ang1[it], &fld->ang2[it]);
            fld->ang1[it] *= 2.0f; //
            fld->ang2[it] *= M_Pi/180.0f; // deg2rad
        }

        // Bond list
        fscanf (fin, "%s %d", saux, &fld->nbonds); 
        fld->iatn = (uint4 *)malloc( fld->nbonds * sizeof(uint4) );
        for (int it = 0; it < fld->nbonds; it++) { 
            fscanf (fin, "%d %d %d %d", &iaux, &fld->iatn[it].w, &fld->iatn[it].x, &fld->iatn[it].y);
        }

        // Angle list
        fscanf (fin, "%s %d", saux, &fld->nangle); 
        fld->iang = (uint4 *)malloc( fld->nangle * sizeof(uint4) );
        for (int it = 0; it < fld->nangle; it++) { 
            fscanf (fin, "%d %d %d %d %d", &iaux, &fld->iang[it].w, &fld->iang[it].x, &fld->iang[it].y, &fld->iang[it].z);
        }

        free(saux);
    fclose(fin);

    return(fld);
}

void FFieldFinish(ffield* fld)
{

    free(fld->nt);
    free(fld->name);
    free(fld->mass);
    free(fld->charge);
    free(fld->lj1);
    free(fld->lj2);
}

// ------ Configuration ------
config* ConfigInput(char *name, float *mass)
{
    int i, iaux, cftyp;
    FILE *fin;

    config* conf;
    conf = (config *)malloc(sizeof(config));

    fin = fopen(strcat(name, ".cfg") , "r");
        fscanf (fin, "%d %d", &conf->natom, &cftyp);
 
        fscanf (fin, "%f %f %f", &conf->box.x, &conf->box.y, &conf->box.z);
        conf->boxi.x = 1.0/conf->box.x;
        conf->boxi.y = 1.0/conf->box.y;
        conf->boxi.z = 1.0/conf->box.z;
 
        // allocate configuration arrays
        conf->itype = (int *)malloc( conf->natom * sizeof(int) );
        conf->pos = (float4 *)malloc( conf->natom * sizeof(float4) );
        conf->vel = (float4 *)malloc( conf->natom * sizeof(float4) );
        conf->force = (float4 *)malloc( conf->natom * sizeof(float4) );
        for (i = 0; i < conf->natom; i++) {
            conf->pos[i].x = 0.0; conf->pos[i].y = 0.0; conf->pos[i].z = 0.0; conf->pos[i].w = 0.0;
            conf->vel[i].x = 0.0; conf->vel[i].y = 0.0; conf->vel[i].z = 0.0; conf->vel[i].w = 0.0;
            conf->force[i].x = 0.0; conf->force[i].y = 0.0; conf->force[i].z = 0.0; conf->force[i].w = 0.0;
        }
 
        for (i = 0; i < conf->natom; i++) {
            fscanf (fin, "%d %d %f", &iaux, &conf->itype[i], &conf->pos[i].w);
            conf->itype[i] = conf->itype[i] - 1;
            conf->vel[i].w = mass[conf->itype[i]];
            fscanf (fin, "%f %f %f", &conf->pos[i].x, &conf->pos[i].y, &conf->pos[i].z);
            if (cftyp > 0)
                fscanf (fin, "%f %f %f", &conf->vel[i].x, &conf->vel[i].y, &conf->vel[i].z);
        }
    fclose(fin);

    return(conf);
}

void ConfigFinish(config* conf, ffield* fld, char *name)
{
    FILE *fout;

    // write an xyz file 
    fout = fopen(strcat(name, ".xyz") , "w");
        fprintf (fout, "%d\n", conf->natom);
        fprintf (fout, "%f %f %f\n", conf->box.x, conf->box.y, conf->box.z);
        for (int i = 0; i < conf->natom; i++)
            fprintf (fout, "%s %f %f %f\n", fld->name[conf->itype[i]], conf->pos[i].x, conf->pos[i].y, conf->pos[i].z);
    fclose(fout);

    // write a configuration file 
    fout = fopen(strcat(name, ".cfg") , "w");
        fprintf (fout, "%d %d\n", conf->natom, 1);
        fprintf (fout, "%f %f %f\n", conf->box.x, conf->box.y, conf->box.z);
        for (int i = 0; i < conf->natom; i++)
        {
            fprintf (fout, "%d %d %f\n", i+1, conf->itype[i]+1, conf->pos[i].w);
            fprintf (fout, "%f %f %f\n", conf->pos[i].x, conf->pos[i].y, conf->pos[i].z);
            fprintf (fout, "%f %f %f\n", conf->vel[i].x, conf->vel[i].y, conf->vel[i].z);
        }
    fclose(fout);

    free(conf->itype);
    free(conf->pos);
    free(conf->force);
}

/* end of simulation.cpp */
