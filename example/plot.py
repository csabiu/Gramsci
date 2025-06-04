import matplotlib.pyplot as plt
import numpy as np

x=np.loadtxt('test.2pcf',skiprows=1)
plt.plot(x[:,1],x[:,-1],'.')
plt.show()
plt.savefig('2pcf.pdf')

x=np.loadtxt('test2.2pcf',skiprows=1)
xi=np.reshape(x[:,-1],(10,10))
plt.contour(np.log10(xi))
plt.show()
plt.savefig('2pcf2.pdf')

x=np.loadtxt('test.3pcf',skiprows=1)
plt.plot(x[:,-1],'.')
plt.show()
plt.savefig('3pcf.pdf')

