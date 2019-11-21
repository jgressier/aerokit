import matplotlib.pyplot as plt

sty_shock = { 'color': 'red',    'linewidth': 3 }
sty_carac = { 'color': 'orange', 'linewidth': 2 }
sty_wall  = { 'color': 'black',  'linewidth': 3 }
sty_flow  = { 'color': 'green',  'linewidth': 3 }
sty_text  = { 'fontsize': 14 }

def set_grid(ax=plt):
	ax.minorticks_on()
	ax.grid(which='major', linestyle='-', alpha=0.8)
	ax.grid(which='minor', linestyle=':', alpha=0.5)

def figure_theta_sigma(**kwargs):
	fig = plt.figure(**kwargs)
	set_grid()
	return fig

def figure_theta_pressure(**kwargs):
	fig = plt.figure(**kwargs)
	set_grid()
	plt.yscale('log')
	return fig

