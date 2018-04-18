import matplotlib.pyplot as plt
import numpy as np

#with open('../data/timelist.txt', 'r') as infile:


brute_world, brute_eth, brute_cancer = [2.329, 26.279, 252.785, 2311.189], [19.830, 199.135, 2096.593, 21004.229], [14.991, 139.855, 1499.483, 14931.224]

kd_world_s = [199.082, 98.287, 58.819, 24.912, 12.314, 6.239, 4.045, 2.975, 4.179, 7.751]
kd_world_c = [1.500, 1.974, 2.597, 3.454, 4.750, 8.054, 12.232, 18.529, 33.025, 64.503]
kd_world_t = []
for i in range (10):
    kd_world_t.append(kd_world_s[i] + kd_world_c[i])
kd_x_world = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
kd_world_ticks = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

kd_eth_s = [492.240, 407.804, 223.579, 171.457, 94.309, 62.528, 51.331, 61.513]
kd_eth_c = [21.828, 25.252, 31.006, 37.528, 44.014, 59.024, 78.377, 120.159]
kd_eth_t = []
for i in range(8):
    kd_eth_t.append(kd_eth_s[i] + kd_eth_c[i])
kd_x_eth = [0, 1, 2, 3, 4, 5, 6, 7]
kd_eth_ticks = [64, 128, 256, 512, 1024, 2048, 4096, 8192]


kd_cancer_s = [418.6271, 446.5863, 517.8098, 618.0976]
kd_cancer_c = [13.168, 16.103, 19.533, 21.240]
kd_cancer_t = []
for i in range(4):
    kd_cancer_t.append(kd_cancer_s[i] + kd_cancer_c[i])
kd_x_cancer = [0, 1, 2, 3]
#kd_cancer_ticks = [16, 32, 64, 128]
kd_cancer_ticks = ['', '16', '', '32', '', '64', '', '128']

#ax = plt.subplot()
#ax.set_color_cycle(['blue', 'yellow', 'red'])

"""plt.xticks(kd_x_world, kd_world_ticks)
plt.xlabel('Number of Clusters')
plt.ylabel('Time (s)')
plt.title('KD-Tree on World Cities')
plt.plot(kd_world_s, label="Search time")
plt.plot(kd_world_c, label="Cluster time")
plt.plot(kd_world_t, label="Total time")
plt.legend()
plt.show()"""

"""plt.xticks(kd_x_eth, kd_eth_ticks)
plt.xlabel('Number of Clusters')
plt.ylabel('Time (s)')
plt.title('KD-Tree on Ethylene Gas Mixtures')
plt.plot(kd_eth_s, label="Search time")
plt.plot(kd_eth_c, label="Cluster time")
plt.plot(kd_eth_t, label="Total time")
plt.legend()
plt.show()"""

"""plt.xticks(kd_x_cancer, kd_cancer_ticks)
plt.xlabel('Number of Clusters')
plt.ylabel('Time (s)')
#plt.title('KD-Tree on Cancer Sample Expressions')
plt.plot(kd_cancer_s, label="Search time")
plt.plot(kd_cancer_c, label="Cluster time")
plt.plot(kd_cancer_t, label="Total time")"""

"""f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax.plot(kd_cancer_s, label="Search time")
ax.plot(kd_cancer_c, label="Cluster time")
ax.plot(kd_cancer_t, label="Total time")
ax2.plot(kd_cancer_s, label="Search time")
ax2.plot(kd_cancer_c, label="Cluster time")
ax2.plot(kd_cancer_t, label="Total time")

# zoom-in / limit the view to different portions of the data
ax.set_ylim(400, 650)  # outliers only
ax2.set_ylim(0, 30)  # most of the data

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax2.set_xticklabels(kd_cancer_ticks)

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xlabel('Number of Clusters')
ax.set_ylabel('Time (das)')
ax2.set_ylabel('Time (s)')
ax.set_title('KD-Tree on Cancer Sample Expressions')
ax.legend(loc="upper left")

#plt.legend()
plt.show()"""

# * * * * * * * * * * * * * * * * 
# Next algorithm
# * * * * * * * * * * * * * * * *

km_world_s = [1725.613, 1154.111, 1004.281, 580.226, 117.626]
km_world_c = [25.026, 143.965, 214.689, 359.110, 1094.297]
km_world_t = []
for i in range (5):
    km_world_t.append(km_world_s[i] + km_world_c[i])
km_x_world = [0, 1, 2, 3, 4]
km_world_ticks = [5, 10, 15, 25, 50]

km_eth_s = [6075.511, 4216.063, 3017.246]
km_eth_c = [2841.531, 4486.410, 8416.390]
km_eth_t = []
for i in range (3):
    km_eth_t.append(km_eth_s[i] + km_eth_c[i])
km_x_eth = [0, 1, 2]
km_eth_ticks = [10, 15, 25]

km_cancer_s = [10778, 10078.791, 5912.269, 6741.223]
km_cancer_c = [433.889, 640.497, 1973.352, 2035.409]
km_cancer_t = []
for i in range (4):
    km_cancer_t.append(km_cancer_s[i] + km_cancer_c[i])
km_x_cancer = [0, 1, 2, 3]
km_cancer_ticks = [15, 26, 50, 72]

"""plt.xticks(km_x_world, km_world_ticks)
plt.xlabel('Number of Clusters')
plt.ylabel('Time (s)')
plt.title('KMeans on World Cities')
plt.plot(km_world_s, label="Search time")
plt.plot(km_world_c, label="Cluster time")
plt.plot(km_world_t, label="Total time")
plt.legend()
plt.show()"""

"""plt.xticks(km_x_eth, km_eth_ticks)
plt.xlabel('Number of Clusters')
plt.ylabel('Time (s)')
plt.title('KMeans on Ethylene Gas Mixtures')
plt.plot(km_eth_s, label="Search time")
plt.plot(km_eth_c, label="Cluster time")
plt.plot(km_eth_t, label="Total time")
plt.legend()
plt.show()"""

"""plt.xticks(km_x_cancer, km_cancer_ticks)
plt.xlabel('Number of Clusters')
plt.ylabel('Time (s)')
plt.title('KMeans on Cancer Sample Expressions')
plt.plot(km_cancer_s, label="Search time")
plt.plot(km_cancer_c, label="Cluster time")
plt.plot(km_cancer_t, label="Total time")
plt.legend()
plt.show()"""

# * * * * * * * * * * * * * * * * 
# Next algorithm
# * * * * * * * * * * * * * * * *

lsh_world_s = [20.8605, 37.2615, 17.9966, 32.5491, 11.3417]
lsh_world_c = [2.442, 2.432, 3.315, 2.457, 2.506]
lsh_world_t = []
for i in range (5):
    lsh_world_t.append(lsh_world_s[i] + lsh_world_c[i])
lsh_x_world = [0, 1, 2, 3, 4]
lsh_world_ticks = ['', '<3, 500>', '', '<3, 1000>', '', '<4, 500>', '', '<3, 840>', '', '<3, 100>',]

lsh_eth_s = [162.6397, 102.19, 406.2016, 498.0375, 60.9259, 472.9079]
lsh_eth_c = [6.816, 8.386, 8.562, 7.140, 6.446, 6.447]
lsh_eth_t = []
for i in range (6):
    lsh_eth_t.append(lsh_eth_s[i] + lsh_eth_c[i])
lsh_x_eth = [0, 1, 2, 3, 4, 5]
lsh_eth_ticks = ['', '<3, 5000>', '<3, 1000>', '<4, 5000>', '<3, 2735>', '<3, 10000>', '<3, 500000>']

lsh_cancer_s = [13.726, 13.750, 14.789, 14.188]
lsh_cancer_c = [2.776, 2.679, 2.871, 2.843]
lsh_cancer_t = []
for i in range (4):
    lsh_cancer_t.append(lsh_cancer_s[i] + lsh_cancer_c[i])
lsh_x_cancer = [0, 1, 2, 3]
lsh_cancer_ticks = ['', '<3, 17>', '', '<3, 12>', '', '<3,5>', '', '<3, 8>']

"""plt.xticks(lsh_x_world, lsh_world_ticks)
plt.xlabel('M-W')
plt.ylabel('Time (s)')
plt.title('LSH on World Cities')
plt.plot(lsh_world_s, label="Search time")
plt.plot(lsh_world_c, label="Cluster time")
plt.plot(lsh_world_t, label="Total time")
plt.legend()
plt.show()"""

"""f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax.plot(lsh_world_s, label="Search time")
ax.plot(lsh_world_c, label="Cluster time")
ax.plot(lsh_world_t, label="Total time")
ax2.plot(lsh_world_s, label="Search time")
ax2.plot(lsh_world_c, label="Cluster time")
ax2.plot(lsh_world_t, label="Total time")

# zoom-in / limit the view to different portions of the data
ax.set_ylim(10, 65)  # outliers only
ax2.set_ylim(0, 5)  # most of the data

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax2.set_xticklabels(lsh_world_ticks)

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xlabel('<M value, W value>')
ax.set_ylabel('Time (das)')
ax2.set_ylabel('Time (s)')
ax.set_title('LSH on World Cities')
ax.legend(loc="upper right")

plt.show()"""


"""plt.xticks(lsh_x_eth, lsh_eth_ticks)
plt.xlabel('M-W')
plt.ylabel('Time (s)')
plt.title('LSH on Ethylene Gas Mixtures')
plt.plot(lsh_eth_s, label="Search time")
plt.plot(lsh_eth_c, label="Cluster time")
plt.plot(lsh_eth_t, label="Total time")
plt.legend()
plt.show()"""

"""f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax.plot(lsh_eth_s, label="Search time")
ax.plot(lsh_eth_c, label="Cluster time")
ax.plot(lsh_eth_t, label="Total time")
ax2.plot(lsh_eth_s, label="Search time")
ax2.plot(lsh_eth_c, label="Cluster time")
ax2.plot(lsh_eth_t, label="Total time")

# zoom-in / limit the view to different portions of the data
ax.set_ylim(50, 550)  # outliers only
ax2.set_ylim(0, 10)  # most of the data

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax2.set_xticklabels(lsh_eth_ticks)

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xlabel('<M value, W value>')
ax.set_ylabel('Time (das)')
ax2.set_ylabel('Time (s)')
ax.set_title('LSH on Ethylene Gas Mixtures')
ax.legend(loc="upper left")

plt.show()"""


"""plt.xticks(lsh_x_cancer, lsh_cancer_ticks)
plt.xlabel('M-W')
plt.ylabel('Time (s)')
plt.title('LSH on Cancer Sample Expressions')
plt.plot(lsh_cancer_s, label="Search time")
plt.plot(lsh_cancer_c, label="Cluster time")
plt.plot(lsh_cancer_t, label="Total time")
plt.legend()
plt.show()"""

"""f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax.plot(lsh_cancer_s, label="Search time")
ax.plot(lsh_cancer_c, label="Cluster time")
ax.plot(lsh_cancer_t, label="Total time")
ax2.plot(lsh_cancer_s, label="Search time")
ax2.plot(lsh_cancer_c, label="Cluster time")
ax2.plot(lsh_cancer_t, label="Total time")

# zoom-in / limit the view to different portions of the data
ax.set_ylim(10, 30)  # outliers only
ax2.set_ylim(0, 5)  # most of the data

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax2.set_xticklabels(lsh_cancer_ticks)

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xlabel('<M value, W value (in millions)>')
ax.set_ylabel('Time (ks)')
ax2.set_ylabel('Time (s)')
ax.set_title('LSH on Cancer Sample Expressions')
ax.legend(loc="upper left")

plt.show()"""

# * * * * * * * * * * * * * * * * 
# Best Graphs
# * * * * * * * * * * * * * * * *

brute = [2311.189, 21004.229, 14931.224]
kd_best_c = [18.529, 78.377, 13.168]
kd_best_s = [2.975, 51.331, 4186.271]
kd_best_t = [16.277, 121.552, 4199.439]
km_best_c = [1094.297, 8416.390, 1973.352]
km_best_c_simp = [109, 842, 197]
km_best_s = [117.626, 3017.246, 5912.269]
km_best_t = [939.336, 8702.473, 7885.621]
lsh_best_c = [2.506, 6.446, 2.776]
lsh_best_s = [113.407, 609.259, 13726]
lsh_best_t = [115.839, 615.705, 13728.679]

#alsh_c = [2.050, 6.076, 2.919]
#alsh_s = [0.515, 1.929, 4.769]
alsh_t = [2.565, 8.005, 7.688]

best_x = [0, 1, 2]
#best_x_ticks = ['', 'cities', '', '', '', 'ethylene', '', '', '', 'cancer']
best_x_ticks = ['cities', 'ethylene', 'cancer']

"""f, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.set_xticklabels(best_x_ticks)
#plt.xticks(best_x, best_x_ticks)
#plt.xlabel('Data Set')
ax1.set_ylabel('Time (s)')
ax1.plot(kd_best_s, label="KD-Tree")
ax1.plot(km_best_s, label="KMeans")
ax1.plot(lsh_best_s, label="LSH")
ax1.legend(loc="upper left")
ax1.set_title('Best Search Time')

ax2.set_xticklabels(best_x_ticks)
#ax2.set_xticklabels(best_x, best_x_ticks)
ax2.set_xlabel('Data Set')
ax2.set_ylabel('Time (s)')
ax2.plot(kd_best_c, label="KD-Tree")
ax2.plot(km_best_c_simp, label="KMeans (das)")
ax2.plot(lsh_best_c, label="LSH")
ax2.legend(loc="upper left")
ax2.set_title('Corresponding Clustering Time')

plt.show()"""

plt.xticks(best_x, best_x_ticks)
plt.xlabel('Data Set')
plt.ylabel('Time (s)')
plt.title('Total Time Taken')
plt.plot(kd_best_t, label="KD-Tree")
plt.plot(km_best_t, label="KMeans")
plt.plot(lsh_best_t, label="LSH")
plt.plot(brute, label="Exhaustive Search")
plt.plot(alsh_t, label="Approximate LSH")
plt.legend()
plt.show()
