################# TrajXCUDA #########################
######### Suresh Kondati Natarajan ##################
##### Lehrstuhl fuer theoretische chemie ############
######## Ruhr Universitaet Bochum ###################

CC=nvcc

Fort=gfortran

fftw3=-lfftw3

libs=-I/home/snataraj/SOFTS/CUDA-7/include

intel=-D__INTEL_COMPILER

arch=-arch compute_20

sources_all=square.cu angle.cu compute_coulomb_charge.cu compute_lifetime.cu compute_tetrahedral-water.cu compute_resident_time.cu compute_adatoms.cu compute_HbondsSSP.cu compute_rdf_frames.cu compute_diffusion_z.cu select_atoms.cu compute_angledist.cu compute_Hbonds.cu compute_Hbonds2.cu compute_Hbonds1.cu compute_diffusion_lifetime.cu compute_densitygrid1.cu compute_densitygrid.cu compute_densityvelocitygrid.cu compute_coverage.cu compute_rdf.cu compute_exchange.cu  compute_Zdist.cu compute_Vdist.cu compute_diffusion.cu compute_diffusion1.cu compute_VDOS.cu read_inputs.cu density_top_avg_pos.cu density_side_avg_pos.cu compute_Ofrag.cu maximum.cu minimum.cu

kernel_all=density_calc.cu density_filter.cu velocity_calc.cu Xdist_calc.cu Vdist_calc.cu boxmat.cu distmat.cu distmatmono.cu SD_calc.cu VAF_calc.cu HBAF_calc.cu Ofragmat.cu Hbondmat.cu Hbondmatmono.cu Hbond_calc.cu covermat.cu covermatmono.cu reduce.cu reducef.cu SDreduce.cu SDreduceHB.cu SDreducef.cu SDreducefHB.cu


objects=$(sources_all:.cu=.o)

kernel_objects=$(kernel_all:.cu=.o)

exe= TrajXCUDA.cux

all: $(sources_all) $(kernel_all) $(exe)
	#$(Fort) unwrapper.f90 -o unwrapper.x
	cp TrajXCUDA.cux unwrapper.x ~/bin/

$(exe): $(objects) $(kernel_objects)
	$(CC) ${fftw3} ${arch} -o $@ $(objects) $(kernel_objects) main.cu ${libs} ${intel}

%.o : %.cu
	$(CC) $(fftw3) -c ${arch} -o $@  $< ${libs} ${intel}
	
clean:
	rm *.o *.cux *.x
	
