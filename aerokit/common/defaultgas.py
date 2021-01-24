# default value of default gamma
_gamma = 1.4
_r     = 287.1

def set_gamma(gam):
	global _gamma
	_gamma = gam

def gamma():
	return _gamma

def set_r_ideal(r):
	global _r 
	_r = r

def r_ideal():
	return _r