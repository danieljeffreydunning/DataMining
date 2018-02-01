import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

with open('pyclustdata.txt', 'r') as infile:
	dim = int(infile.readline().strip())
	ndata = int(infile.readline().strip())
	data = []
	for i in range(dim * ndata):
		data.append(float(infile.readline().strip()))		
	data = np.reshape(data, (ndata, dim))
        
	clusters = []
	for i in range(ndata):
		clusters.append(int(infile.readline().strip()))

cluster_count = len(set(clusters)) // 2
sns.set_style("whitegrid")
cmap = ListedColormap(sns.color_palette("nipy_spectral", cluster_count)).colors
cmap += ListedColormap(sns.color_palette("terrain", cluster_count)).colors
np.random.shuffle(cmap)
cmap = ListedColormap(cmap)

fig, ax = plt.subplots()
ax.scatter(data[:, 0], data[:, 1], c=clusters, cmap=cmap)
fig.savefig('../plots/clusters.png')
plt.show()
