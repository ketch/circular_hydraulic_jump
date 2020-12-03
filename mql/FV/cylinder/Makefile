all: shallow_hllemcc_2D.so shallow_es_2D.so 

shallow_hllemcc_2D.so:
	#f2py -c rpn2_sw_hllemcc.f90 rpt2_dummy.f90 -m shallow_hllemcc_2D
	f2py -c rpn2_sw_hllemcc.f90 rpt2_shallow_roe_with_efix.f90 -m shallow_hllemcc_2D

shallow_es_2D.so:
	f2py -c rpn2_shallow_es.f90 rpt2_dummy.f90 -m shallow_es_2D

clean:
	rm *.so	
