from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import csv


#data=np.genfromtxt("C:/Users/Spulido99/Documents/PhD_big/data/phac/HT_hiII14simulation_output.parameters.analysis.txt",delimiter='\t').astype('float')
#X = data[:,0]
#Y = data[:,1]
#Z = data[:,2]


Z=np.genfromtxt("C:/Users/Spulido99/Dropbox/PhD/PhaC/convergence_graph.csv",delimiter=',').astype('float')
X = np.arange(0, len(Z[0])) * 50 * 20
Y = np.arange(0, len(Z))

print X

X, Y = np.meshgrid(X, Y)
#Z = Z.flatten()

print len(Z)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#X, Y, Z = axes3d.get_test_data(0.05)


ax.plot_surface(X, Y, Z, linewidth=0.5)



#ax.plot_trisurf(X, Y, Z, cmap=cm.jet, linewidth=0)
#ax.scatter(X, Y, Z)
plt.show()

