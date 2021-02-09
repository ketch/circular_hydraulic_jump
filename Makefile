all: shallow_hllemcc_2D.so shallow_quad_hllemcc_2D.so shallow_es_2D.so

shallow_quad_hllemcc_2D.so: rpn2_swq_hllemcc.f90 rpt2_swq.f90
	f2py -c rpn2_swq_hllemcc.f90 rpt2_swq.f90 -m shallow_quad_hllemcc_2D
#	f2py -c rpn2_swq_hllemcc.f90 rpt2_dummy.f90 -m shallow_quad_hllemcc_2D

shallow_hllemcc_2D.so: rpn2_sw_hllemcc.f90 rpt2_shallow_roe_with_efix.f90 
	f2py -c rpn2_sw_hllemcc.f90 rpt2_shallow_roe_with_efix.f90 -m shallow_hllemcc_2D

shallow_es_2D.so: rpn2_shallow_es.f90 rpt2_dummy.f90 
	f2py -c rpn2_shallow_es.f90 rpt2_dummy.f90 -m shallow_es_2D

clean:
	rm *.so	
