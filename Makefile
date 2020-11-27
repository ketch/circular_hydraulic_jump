all: shallow_hllemccRoEF_2D.so shallow_hllemcc_2D.so shallow_quad_hllemcc_2D.so shallow_es_2D.so

#shallow_hllem_2D.so:
#	f2py -c rpn2_shallow_hllem.f90 rpt2_dummy.f90 -m shallow_hllem_2D

#shallow_hllc_2D.so:
#	f2py -c rpn2_shallow_hllc.f90 rpt2_dummy.f90 -m shallow_hllc_2D

shallow_quad_hllemcc_2D.so:
	f2py -c rpn2_swq_hllemcc.f90 rpt2_dummy.f90 -m shallow_quad_hllemcc_2D

shallow_hllemcc_2D.so:
	f2py -c rpn2_sw_hllemcc.f90 rpt2_dummy.f90 -m shallow_hllemcc_2D

shallow_hllemccRoEF_2D.so:
	f2py -c rpn2swq_hllemccroef.f90 rpt2_dummy.f90 -m shallow_hllemccRoEF_2D

shallow_es_2D.so:
	f2py -c rpn2_shallow_es.f90 rpt2_dummy.f90 -m shallow_es_2D

clean:
	rm *.so	
