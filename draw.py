import numpy as np
from matplotlib.pyplot import *

maxt = 9000
maxt += 1
x_range = np.linspace(-1000, 1000, 100)
y_range = np.linspace(-1000, 1000, 100)
x, y = np.meshgrid(x_range, y_range)

for t in range(0, maxt, 50):
	string_t = str(t).zfill(4)
	print('drawing @ t = ' + string_t)
	#vor = np.loadtxt('vors/vors'+string_t+'.dat')#, allow_pickle=True)
	#vor = vor.reshape(100, 100)
	#contourf(x, y, vor, levels=50, vmin=0., vmax=0.001)
	#colorbar()
	#xlabel('km')
	#ylabel('km')
	#title('vorticity field @ t = ' + string_t)
	#savefig('vor'+string_t+'.png', dpi=300)
	#close()

	phi = np.loadtxt('phis/phis'+string_t+'.dat')
	phi = phi.reshape(100, 100)
	contourf(x, y, phi, levels=20, vmin=-50., vmax=50.)
	colorbar()
	xlabel('km')
	ylabel('km')
	title('streamline field @ t = ' + string_t)
	savefig('phi'+string_t+'.png')
	close()
