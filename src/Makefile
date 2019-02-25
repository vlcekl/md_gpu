#
#  File name: Makefile
#  Date:      2009/03/31 23:45
#  Author:    Lukas Vlcek
#
#  Copyright (C) 2009 Lukas Vlcek
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  in a file called COPYING along with this program; if not, write to
#  the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
#  02139, USA.
#

target  :=  athom
opt  ?= -O3
warn ?= -Wall

cudaroot  := /usr/local/cuda

#----- C++ compilation options ------
gpp         := /usr/bin/g++
ccflags     += $(opt) $(warn)
cclib_paths :=
ccinc_paths := -I $(cudaroot)/include
cclibraries :=

#----- CUDA compilation options -----

nvcc        := $(cudaroot)/bin/nvcc
cuflags     += $(opt) #cuflags     += $(opt) -Xptxas -fastimul
cuflagst    := $(opt) --ptxas-options=-v # -keep --cubin #cuflagsd    := $(opt) -ptx
cuflagsd    := -O0 -deviceemu # -keep --cubin

culib_paths := -L$(cudaroot)/lib
cuinc_paths := -I$(cudaroot)/include
culibraries := -lcuda -lcudart

lib_paths   := $(culib_paths) $(cclib_paths) $(clib_paths)
libraries   := $(culibraries) $(cclibraries) $(clibraries)

#ccsources := simulation.cpp athom.cpp
cusources := gpu_properties.cu config.cu neighborlist.cu forces.cu integrator.cu runGPU.cu prints.cu simulation.cu athom.cu
hsources := gpu_properties.h config.h neighborlist.h forces.h integrator.h runGPU.h prints.h simulation.h
msources  := Makefile ../README

objects := gpu_properties.o config.o neighborlist.o forces.o integrator.o runGPU.o prints.o simulation.o athom.o

$(target): $(objects) 
	$(nvcc) -o $@ $(lib_paths) $(libraries) $(objects)

%.o: %.cu
	$(nvcc) -c $^ $(cuflags) $(cuinc_paths) -o $@ 

#%.o: %.cpp
#	$(gpp) -c $^ $(ccflags) $(ccinc_paths) -o $@

#----- Test occupancy of registers and shared memory
objectst := gpu_propertiest.o configt.o neighborlistt.o forcest.o integratort.o runGPUt.o printst.o simulationt.o athomt.o

test: $(objectst) 
	$(nvcc) -o $@ $(lib_paths) $(libraries) $(objectst)

%t.o: %.cu
	$(nvcc) -c $^ $(cuflagst) $(cuinc_paths) -o $@ 

#----- Compile debugging version
objectsd := gpu_propertiesd.o configd.o neighborlistd.o forcesd.o integratord.o runGPUd.o printsd.o simulationd.o athomd.o

debug: $(objectsd) 
	$(nvcc) -o $@ $(lib_paths) $(libraries) $(objectsd)

%d.o: %.cu
	$(nvcc) -c $^ $(cuflagsd) $(cuinc_paths) -o $@ 

#---- packs all sources + this Makefile + README to a tarball
tar: $(ccsources) $(cusources) $(hsources) $(msources)
	tar czvf athom_n.tgz $^

#----- Clean rubbish
clean:
	rm -f *.o 

#----- Clean rubbish
cleanall:
	rm -f *.o $(target) test debug

info:
	@echo "make [athom]  ..... makes binary athom"
	@echo "make debug    ..... compiles a debugging version for cpu"
	@echo "make test     ..... prints device memory usage"
	@echo "make tar      ..... packs all sources, Makefile, and README to athom_n.tgz"
	@echo "make clean    ..... cleans rubbish, a.k.a. object files"
	@echo "make cleanall ..... cleans object files and old binaries"
	@echo "make info     ..... prints this useless info. More in ../README."

# end of Makefile
