FC      = mpiifort #ftn 
FCFLAGS = -O4 #-g #-fbacktrace -fbounds-check -fno-range-check
LDFLAGS = -qopenmp #-h omp 
TARGET  = do.exe

OBJ  = mt19937.o \
	constant_parameter.o \
	lshell_setting.o \
	von_Neumann.o \

all: $(TARGET)

.SUFFIXES : .f90

$(TARGET): $(OBJ)
	$(FC) -o $@ $(OBJ) $(FCFLAGS) $(LDFLAGS)

.f90.o:
	$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) -c $<

clean:
	rm -f *.o  *.mod *.exe fort.* *.out