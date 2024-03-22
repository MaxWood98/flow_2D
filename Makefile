#build arguments  #-Wall -fbounds-check -ffpe-trap=underflow,zero -fopt-info-optimized=$@_opt.dat
# buildargs = -O0 -Wall -fbounds-check -ffpe-trap=underflow,zero,invalid -fopenmp
# buildargs = -O0 -Wall -fbounds-check -ffpe-trap=inexact,denormal -fopenmp
# buildargs = -O0 -Wall -fbounds-check -ffpe-trap=denormal -fopenmp 
buildargs = -O2 -fopt-info-optimized=$@_opt.dat -fopenmp

#build settings
buildsettings = -J obj 

#object file directory
OBJDIR = obj/

#list object file names -> 1 for each .f90 in the correct compilation order
OBJS = $(addprefix $(OBJDIR), \
		flow_data.o\
		flow_general.o\
		io_utilities_module.o\
		flow_io.o\
		flow_mesh.o\
		flow_boundaries.o\
		flow_calculation.o\
		flow_iterator.o\
		)

#object patturn rule -> for every file in $(OBJDIR) of the form var.o make it from src/var.f90
$(OBJDIR)%.o : src/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

#main build procedure
build: flow_link

#linking procedure
flow_link: $(OBJS) $(addprefix $(OBJDIR), flow_main.o)
	gfortran -o flow2d $^ $(buildsettings) -I obj $(buildargs) 

#clean procedure 
clean: 
	rm obj/*.o 
	rm obj/*.mod
	rm obj/*.dat