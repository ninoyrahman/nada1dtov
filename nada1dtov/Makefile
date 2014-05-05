EXED = build/
SRCD = src/

SRCS =tmp_mdparam.f90 tmp_mdgrid.f90 tmp_mdmetric.f90 tmp_mdhydro.f90 MD_Boundary.f90 MD_HydroSubroutines.f90 PPM_RECONSTRUCTION.f90 Schwarzschild_isotropic.f90 main_BSSN_1D.f90 Grid1D.f90  Output_1D.f90 Derivatives.f90 Einstein.f90 Derivatives_func_1.f90 Derivatives_func_12.f90 Einstein_2nd.f90 Runge_kutta_2nd_BCs.f90 Gauge_dyn.f90 Einstein_2nd_trK.f90 TOV_source.f90 Einstein_2nd_deltax.f90 Einstein_2nd_capB.f90 Einstein_2nd_beta.f90 Eigenvalues.f90 Fluxespm.f90 Con2Prim.f90 Output_Max.f90 BH_ID.f90 Shock.f90

#OBJS =	tmp_mdparam.o tmp_mdgrid.o tmp_mdmetric.o tmp_mdhydro.o MD_Boundary.o MD_HydroSubroutines.o Schwarzschild_isotropic.o main_BSSN_1D.o Grid1D.o  Output_1D.o Derivatives.o Einstein.o Derivatives_func_1.o Derivatives_func_12.o Einstein_2nd.o Runge_kutta_2nd_BCs.o Gauge_dyn.o Einstein_2nd_trK.o TOV_source.o Einstein_2nd_deltax.o Einstein_2nd_capB.o Einstein_2nd_beta.o Eigenvalues.o Fluxespm.o Con2Prim.o Output_Max.o BH_ID.o Shock.o

OBJS=$(SRCS:.f90=.o)

BINOBJECTS=$(patsubst %.o,$(EXED)%.o,$(OBJS))

EXECTBL=nada1D

LIB_PATH = 

LIBS =	

F90 = gfortran

F90FLAGS = -fopenmp


all: $(EXECTBL) move


$(EXECTBL): $(OBJS)
	@echo '' 
	@echo '============ ...building the executable ============'
	@echo ''
	mkdir -p build
	mkdir -p data
	$(F90) $(F90FLAGS) $(BINOBJECTS) -o $(EXED)$@
	@echo ''
	@echo '============ ...Done! ============'

$(OBJS):%.o : $(SRCD)%.f90
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	$(F90) $(F90FLAGS) -c $< -o $(EXED)$@ 

move:
	mv *.mod $(EXED)

clean:
	rm -f $(BINOBJECTS) $(EXED)$(EXECTBL) $(EXED)*.mod

flush:
	rm -f data/*.dat

