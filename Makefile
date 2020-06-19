#shallow_hllem_2D.so:
#	f2py -c rpn2_shallow_hllem.f90 rpt2_dummy.f90 -m shallow_hllem_2D
		
shallow_hllc_2D.so:
	f2py -c rpn2_shallow_hllc.f90 rpt2_dummy.f90 -m shallow_hllc_2D

clean:
	rm *.so
