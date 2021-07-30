FC      = ftn
FCFLAGS = -g #-fast -xMIC-AVX512 
LDFLAGS = -h omp
TARGET  = do.exe

OBJ  = mt19937.o \
       constant_parameter.o \
       lshell_setting.o \
       constants_in_the_simulations.o \
       variables.o \
       subroutines.o \
       main.o \

all: $(TARGET)

.SUFFIXES : .f90
.f90.o:
	$(FC) $(FCFLAGS) $(LDFLAGS) -c $< 

$(TARGET): $(OBJ)
	$(FC) -o $@ $(OBJ) $(FCFLAGS) $(LDFLAGS)

clean:
	rm -f *.o  *.mod *.exe fort.* *.out
