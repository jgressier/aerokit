import matplotlib.pyplot as plt

def set_grid():
	plt.minorticks_on()
	plt.grid(which='major', linestyle='-', alpha=0.8)
	plt.grid(which='minor', linestyle=':', alpha=0.5)

def figure_theta_sigma(**kwargs):
	fig = plt.figure(**kwargs)
	set_grid()
	return fig

def figure_theta_pressure(**kwargs):
	fig = plt.figure(**kwargs)
	set_grid()
	plt.yscale('log')
	return fig

