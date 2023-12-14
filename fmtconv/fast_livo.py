# convert fast-livo result file to the traj_eval input format
# zongzhou.wu   Dec 14, 2023

# ----- input format example -----
# fout_out << setw(20) << std::fixed<<std::setprecision(6) << LidarMeasures.last_update_time  << " " << euler_cur.transpose()*57.3 << " " << state.pos_end.transpose() << " " << state.vel_end.transpose() \
#          <<" "<<state.bias_g.transpose()<<" "<<state.bias_a.transpose()<<" "<<state.gravity.transpose()<<" "<<feats_undistort->points.size()<<endl;\

# t(s) angle(deg deg deg) pos(m m m) vel(m/s m/s m/s) bg() ba() g(m/s2 m/s2 m/s2) fea_num
# 1698131923.699656 -0.018043  0.301528  0.067456 -0.000270  0.000074  0.000832  0.005261 -0.001667  0.007317 -0.001835 -0.001222 -0.004142 0.000000 0.000000 0.000000  9.797947 -0.005488 -0.486102 5552
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



input_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/lio/mat_out.txt"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/predict_odom_method_1.csv"

input_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/livo/mat_out.txt"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/predict_odom_method_2.csv"

input_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/mat_out.txt"
input_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/mat_pre.txt"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/predict_odom_method_3.csv"

# follow the x y z order
# fast-livo output rpy
result = np.loadtxt(input_path,usecols=[0,5,6,4,2,1,3])

result[:,4:7] = result[:,4:7]/180*pi

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

# check the pos xyz order
# pos = result[:,1:4]
# fig, ax = plt.subplots()
# ax.plot(pos[:,0], pos[:,1], linewidth = 2)
# plt.show()
