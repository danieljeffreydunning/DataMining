import matplotlib.pyplot as plt
import numpy as np

#with open('../data/timelist.txt', 'r') as infile:


brute_world, brute_eth, brute_cancer = [2.329, 26.279, 252.785, 2311.189], [19.830, 199.135, 2096.593, 21004.229], [14.991, 139.855, 1499.483, 14931.224]

kd32_world, kd128_world, kd512_world, kd2048_world, kd8192_world = [0.083, 1.012, 9.883, 98.287], [0.017, 0.253, 2.428, 24.912], [0.005, 0.069, 0.672, 6.239], [0.003, 0.031, 0.299, 2.975], [0.008, 0.079, 0.764, 7.751]
km5_world, km10_world, km15_world, km25_world, km50_world = [1.674, 17.368, 179.306, 1725.613], [1.200, 12.511, 116.591, 1154.111], [1.158, 11.056, 105.885, 1004.281], [0.531, 6.216, 59.363, 580.226], [0.174, 1.762, 16.799, 172.626]

