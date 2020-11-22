all: shallow_llf_2D.so

shallow_llf_2D.so:
	f2py -c rpn2_shallow_hlle.f90 rpt2_dummy.f90 -m shallow_llf_2D

clean:
	rm *.so	
