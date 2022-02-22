# test

import numpy as np
import aerokit.engine.propu_be_double_corps as tfinl
import aerokit.engine.Cycle_Turbofan        as tfclass

bpr=6.
opr=41.
Tt4=1620.

g      = tfclass.Gaz()
g_fuel = tfclass.Gaz(1.33,291.6,42800.e3)

for opr in [ 10., 20., 40., 50. ]:
	print "OPR: ",opr
	#tf_fan = tfclass.cycle_taux_fan_fixe(6.1,32.8,1410.,180.,g,g_fuel,1.6)
	tf_pow = tfclass.cycle_taux_fan_calcule(bpr, opr, Tt4, 230.,g,g_fuel,0.58)
	#
	#tf_fan.calculs_5_to_9()
	tf_pow.calculs_5_to_9()
	#
	F_spe, mk_spe, eta_th, eta_prop, eta, pi_f = tfinl.calculs(bpr, opr, Tt4, 230.,0.58)
	#
	print "  F spe ",F_spe, tf_pow.F_spe
	print "  eta th", eta_th, tf_pow.eta_th