FC      = mpiifort #ftn 
FCFLAGS = #-g #-fbacktrace -fbounds-check -fno-range-check
LDFLAGS = -qopenmp #-h omp 
TARGET  = do.exe

OBJ  = mt19937.o \
		constant_parameter.o \
		lshell_setting.o \
		von_Neumann.o \

all: $(TARGET)

.SUFFIXES : .f90
.f90.o:
		$(FC) $(FCFLAGS) $(LDFLAGS) -c $< 

$(TARGET): $(OBJ)
		$(FC) -o $@ $(OBJ) $(FCFLAGS) $(LDFLAGS)

clean:
		rm -f *.o  *.mod *.exe fort.* *.out