#all: shallow_hllemccRoEF_2D.so shallow_hllemcc_2D.so shallow_llf_2D.so shallow_ev_2D.so shallow_es_2D.so shallow_roe_with_efix_2D.so shallow_llf_dk_2D.so

all: shallow_hllemccRoEF_2D.so shallow_hllemcc_2D.so shallow_roe_exp_2D.so shallow_llf_2D.so shallow_ev_2D.so shallow_es_2D.so shallow_roe_with_efix_2D.so shallow_llf_dk_2D.so

shallow_hllemcc_2D.so:
	f2py -c rpn2_sw_hllemcc.f90 rpt2_dummy.f90 -m shallow_hllemcc_2D

shallow_hllemccRoEF_2D.so:
	f2py -c rpn2swq_hllemccroef.f90 rpt2_dummy.f90 -m shallow_hllemccRoEF_2D

shallow_roe_exp_2D.so:
	f2py -c rpn2_shallow_roe_exp.f90 rpt2_dummy.f90 -m shallow_roe_exp_2D

shallow_llf_2D.so:
	f2py -c rpn2_shallow_llf.f90 rpt2_dummy.f90 -m shallow_llf_2D

shallow_llf_dk_2D.so:
	f2py -c rpn2_shallow_llf_dk.f90 rpt2_dummy.f90 -m shallow_llf_dk_2D

shallow_ev_2D.so:
	f2py -c rpn2_shallow_ev.f90 rpt2_dummy.f90 -m shallow_ev_2D

shallow_es_2D.so:
	f2py -c rpn2_shallow_es.f90 rpt2_dummy.f90 -m shallow_es_2D

shallow_roe_with_efix_2D.so:
	f2py -c rpn2_shallow_roe_with_efix.f90 rpt2_dummy.f90 -m shallow_roe_with_efix_2D

clean:
	rm *.so	
