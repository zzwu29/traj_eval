# convert gnss .flt result file to the traj_eval input format
# zongzhou.wu   Dec 28, 2023

# ----- input format example -----
#Seconds of Week      X-ECEF           Y-ECEF          Z-ECEF    Vx-ECEF    Vy-ECEF    Vz-ECEF     X-RMS     Y-RMS     Z-RMS    Vx-RMS    Vy-RMS    Vz-RMS  NSat  PDOP   sigma0  AmbStatus      Ratio         BL  Quality
#            (s)          (m)             (m)             (m)      (m/s)      (m/s)      (m/s)       (m)       (m)       (m)     (m/s)     (m/s)     (m/s)   (#)   (#)      (m)                             (m)        
# ----- input format example -----


# ----- output format example -----
# 1. result
# t(ns),seq,t(ns),pos(m m m),quat(x y z w),0,0,0,0,0,0,
# 1609060336768428032,0,1609060336768428032,-0.000043,-0.000466,-0.000104,0.000033,0.000106,0.000061,1.000000,0.0,0.0,0.0,0.0,0.0,0.0,
# 2. groundtruth
# %time,field.header.seq,field.header.stamp,field.pose.position.x,field.pose.position.y,field.pose.position.z,field.pose.orientation.x,field.pose.orientation.y,field.pose.orientation.z,field.pose.orientation.w
# t(ns),seq,t(ns),pos(m m m),quat(x y z w)
# 1609060363469351911,2734,1609060363469351911,0.000352647179334,-0.000128287925993,-0.000155817087489,0.0,0.0,0.0,1.0
# ----- output format example -----

import numpy as np
import matplotlib.pyplot as plt
from math import *

input_path = "D:/Datasets/HKisland_GNSS01/result_rtk/UAV0-GEC-rtk.flt"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/predict_odom_method_2.csv"


result = np.loadtxt(input_path,usecols=[0,1,2,3])
week = 2285
# https://zhuanlan.zhihu.com/p/383645111
err = 315964800 # start time difference

# http://www.leapsecond.com/java/gpsclock.htm
leapsec = 18 #GPS time was zero at 0h 6-Jan-1980 and since it is not perturbed by leap seconds GPS is now ahead of UTC by 18 seconds.

t_gps_sec = week*604800+result[:,0] + err - leapsec
result[:,0] = t_gps_sec

row_size = np.shape(result)[0]

output_result = np.array(
    [
        result[:,0] * 1e9, 1 + np.arange(row_size), result[:,0] * 1e9,
        result[:,1], result[:,2], result[:,3],
        np.zeros(row_size), np.zeros(row_size), np.zeros(row_size), np.ones(row_size),
        np.zeros(row_size), np.zeros(row_size), np.zeros(row_size),
        np.zeros(row_size), np.zeros(row_size), np.zeros(row_size)
    ]
).T

np.savetxt(output_path, output_result, delimiter = ',', fmt = '%.6f')

