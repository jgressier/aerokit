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

def save_default():
	global _saved_r, _saved_gamma
	_saved_r     = _r
	_saved_gamma = _gamma

def restore_default():
	global _saved_r, _saved_gamma, _r, _gamma
	_r     = _saved_r
	_gamma = _saved_gamma