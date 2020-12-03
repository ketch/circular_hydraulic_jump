all: shallow_annular_hllemcc_2D.so shallow_annular_es_2D.so

shallow_annular_hllemcc_2D.so:
	f2py -c rpn2_sw_annular_hllemcc.f90 rpt2_dummy.f90 -m shallow_annular_hllemcc_2D
	#f2py -c rpn2_sw_annular_hllemcc.f90 rpt2_shallow_roe_with_efix_annulus.f90 -m shallow_annular_hllemcc_2D	

shallow_annular_es_2D.so:
	f2py -c rpn2_sw_annular_es.f90 rpt2_dummy.f90 -m shallow_annular_es_2D

clean:
	rm *.so	
