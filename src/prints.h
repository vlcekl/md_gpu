/*
   File name: prints.h
   Date:      2009/04/03 20:29
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


#ifndef __PRINTS_H__
#define __PRINTS_H__

void print_statistics( unsigned int cur_timestep, integr* integ );
void print_random_stuff( params* pars, integr* integ, ffield* fld, config* conf, gpu_struct* g_prop, gpu_config* g_conf, gpu_integr* g_integ, gpu_ffield* g_fld, int id );
void print_ave_force_momentum(gpu_config* g_conf, config* conf);

#endif

/* end of prints.h */
