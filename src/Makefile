FC=ifort

target1=nagcap
target2=nagcap_plot

src=$(wildcard *.f90)
obj=$(patsubst %.f90, %.o, $(src))

subs = get_par.o sacio.o load_data.o detrend.o ftiout.o conv.o findpeaks.o raytracing.o triangle.o align_seis.o butter.o nmtensor.o
mod = raytracing.mod sacio.mod load_data.mod
Flags=-qopenmp

all : $(target1) $(target2)

$(target1) : $(subs) fitness.o NA.o nagcap.o
	$(FC) $^ $(Flags) -o $@

$(target2) : $(subs) sort.o gen_p_data.o nagcap_plot.o
	$(FC) $^ $(Flags) -o $@

%.o : %.f90 
	$(FC) $^ $(Flags) -c 

%.mod : %.f90
	$(FC) $^ -c

.PHONY : clean
clean :
	rm *.o $(target1) $(target2)
