/*
   File name: athom.cpp
   Date:      2009/03/31 20:24
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

#include <string>
#include <iostream>
#include "simulation.h"

using namespace std;

int main (int argc, char * argv[])
{
	if( argc != 2 ) {
		cout << "USAGE: In a directory containing files\n"
                     << "simulation.inp simulation.fld and simulation.cfg\n"
                     << "run ./athom simulation\n";
                return 1;
	}

        simulation* sim;
        sim = SimulationInit(argv[1]); // read input, initialize, and pack in 'sim' structure

	SimulationRun(sim); // run a simulation

        SimulationFinish(sim); // print output and finalize

	return 0;
}


/* end of athom.cpp */
