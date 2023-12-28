# convert integrated navigation .ins result file to the traj_eval input format
# zongzhou.wu   Dec 28, 2023

# ----- input format example -----
# Seconds of Week            X-ECEF            Y-ECEF            Z-ECEF        VX        VY        VZ     Pitch      Roll       Yaw      GyroBiasX      GyroBiasY      GyroBiasZ      AcceBiasX      AcceBiasY      AcceBiasZ    MeasType   Nsat   PDOP   AmbStatus         
#             (s)               (m)               (m)               (m)     (m/s)     (m/s)     (m/s)     (deg)     (deg)     (deg)        (deg/h)        (deg/h)        (deg/h)           (mg)           (mg)           (mg)                  #      #            
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

def a2qua(att=[]): #att to xyzw quat
    pitch, roll, yaw = att[0]/2.0, att[1]/2.0, att[2]/2.0
    sp, sr, sy = sin(pitch), sin(roll), sin(yaw)
    cp, cr, cy = cos(pitch), cos(roll), cos(yaw)
    
    q_xyzw = np.array([
        sp*cr*cy - cp*sr*sy, #x
        cp*sr*cy + sp*cr*sy, #y
        cp*cr*sy + sp*sr*cy, #z
        cp*cr*cy - sp*sr*sy  #w
    ])
    
    return q_xyzw

def a2qua_batch(att=[]): #att to xyzw quat
    q_xyzw = np.zeros([np.shape(att)[0],4])
    for i in range(np.shape(att)[0]):
        [q_xyzw[i,0],q_xyzw[i,1],q_xyzw[i,2],q_xyzw[i,3]]=a2qua(att[i,:])
    return q_xyzw


input_path = "D:/Datasets/HKisland_GNSS01/result_ipn/UAV0-GEC-rtkins.ins"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/predict_odom_method_3.csv"

input_path = "D:/Datasets/HKisland_GNSS01/result_ipn/UAV0-GEC-rtkinsvis.ins"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/predict_odom_method_4.csv"


result = np.loadtxt(input_path,usecols=[0,1,2,3,7,8,9])
result[:,4:7] = result[:,4:7]/180*pi

week = 2285
# https://zhuanlan.zhihu.com/p/383645111
err = 315964800 # start time difference

# http://www.leapsecond.com/java/gpsclock.htm
leapsec = 18 #GPS time was zero at 0h 6-Jan-1980 and since it is not perturbed by leap seconds GPS is now ahead of UTC by 18 seconds.

t_gps_sec = week*604800+result[:,0] + err - leapsec
result[:,0] = t_gps_sec

row_size = np.shape(result)[0]

q_xyzw = a2qua_batch(result[:,4:7])

output_result = np.array(
    [
        result[:,0] * 1e9, 1 + np.arange(row_size), result[:,0] * 1e9,
        result[:,1], result[:,2], result[:,3],
        # np.zeros(row_size), np.zeros(row_size), np.zeros(row_size), np.ones(row_size),
        q_xyzw[:,0],q_xyzw[:,1],q_xyzw[:,2],q_xyzw[:,3],
        np.zeros(row_size), np.zeros(row_size), np.zeros(row_size),
        np.zeros(row_size), np.zeros(row_size), np.zeros(row_size)
    ]
).T

np.savetxt(output_path, output_result, delimiter = ',', fmt = '%.6f')

