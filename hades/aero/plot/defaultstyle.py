import matplotlib.pyplot as plt

def set_grid():
	plt.minorticks_on()
	plt.grid(which='major', linestyle='-', alpha=0.8)
	plt.grid(which='minor', linestyle=':', alpha=0.5)