# convert navsatfix file to the traj_eval input format
# zongzhou.wu   Dec 14, 2023

# rostopic echo -b file_name.bag -p NavSatFix_topic > gt.csv


# ----- input format example -----
# %time,field.header.seq,field.header.stamp,field.header.frame_id,field.status.status,field.status.service,field.latitude,field.longitude,field.altitude,field.position_covariance0,field.position_covariance1,field.position_covariance2,field.position_covariance3,field.position_covariance4,field.position_covariance5,field.position_covariance6,field.position_covariance7,field.position_covariance8,field.position_covariance_type
# 1698131763008878188,869,1698131763008531176,,0,0,22.20785535829121,114.26042499370304,5.619030475616455,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0
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

glv_a = 6378137;   # the Earth's semi-major axis
glv_f  = 1/298.257;# flattening
glv_wie = 7.2921151467e-5; 
glv_b = (1- glv_f)* glv_a;  # semi-minor axis
glv_e = sqrt(2* glv_f- glv_f**2);  
glv_e2 =  glv_e**2; # 1st eccentricity
glv_ep = sqrt( glv_a**2- glv_b**2)/ glv_b; 
glv_ep2 =  glv_ep**2; # 2nd eccentricity

def blh2xyz(blh=[]):
    b, l, h=blh[0], blh[1], blh[2]

    sb,cb, sl, cl = sin(b), cos(b), sin(l), cos(l)
    N = glv_a/sqrt(1-glv_e2*sb**2)
    x = (N+h)*cb*cl
    y = (N+h)*cb*sl
    z = (N*(1-glv_e2)+h)*sb
    return np.array([x,y,z])

def blh2xyz_batch(blh=[]):
    x_all,y_all,z_all=[],[],[]
    for i in range(np.shape(blh)[0]):
        x, y, z = blh2xyz(blh[i,:])
        x_all.append(x)
        y_all.append(y)
        z_all.append(z)

    return np.array([x_all,y_all,z_all]).T


def xyz2blh(xyz=[]):
    x, y, z=xyz[0], xyz[1], xyz[2]

    bell = glv_a*(1.0-1.0/glv_f)                          
    ss = sqrt(x*x+y*y)   
    zps   = z/ss
    theta = atan( (z*glv_a) / (ss*glv_b) )
    sin3  = sin(theta) * sin(theta) * sin(theta)
    cos3  = cos(theta) * cos(theta) * cos(theta)
    
    #Closed formula
    b = atan((z + glv_ep2 * glv_b * sin3) / (ss - glv_e2 * glv_a * cos3))
    l = atan2(y,x)
    nn = glv_a/sqrt(1.0-glv_e2*sin(b)*sin(l))
    h = ss/cos(b) - nn

    i=0
    while i<=100:
        nn = glv_a/sqrt(1.0-glv_e2*sin(b)*sin(b))
        hOld = h
        phiOld = b
        h = ss/cos(b)-nn
        b = atan(zps/(1.0-glv_e2*nn/(nn+h)))
        if abs(phiOld-b) <= 1.0e-11 and abs(hOld-h) <= 1.0e-5:
            # always convert longitude to 0-360
            if l < 0.0 :
                l += 2 * pi
                break

        i+=1

    return np.array([b,l,h])

def xyz2enu(xyz=[],xyz_ref=[]):

    [b,l,h]=xyz2blh(xyz_ref)

    r=[xyz[0]-xyz_ref[0], xyz[1]-xyz_ref[1], xyz[2]-xyz_ref[2]]

    sinPhi = sin(b) # lat
    cosPhi = cos(b)
    sinlam = sin(l) # lon
    coslam = cos(l)

    n = -sinPhi * coslam * r[0] - sinPhi * sinlam * r[1] + cosPhi * r[2]
    e = -sinlam * r[0] + coslam * r[1]
    u = +cosPhi * coslam * r[0] + cosPhi * sinlam * r[1] + sinPhi * r[2]

    return np.array([e,n,u])

def xyz2enu_batch_fixed(xyz=[],xyz_ref=[]): # w.r.t a fixed point
    enu=np.zeros(np.shape(xyz))

    for idx in range(xyz.shape[0]):
        enu_idx = xyz2enu(xyz[idx,:],xyz_ref)
        enu[idx, 0], enu[idx, 1], enu[idx, 2] = enu_idx[0], enu_idx[1], enu_idx[2]
    
    return enu


# sensor_msgs/NavSatFix

input_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/gnss/gt.csv"
output_path = "C:/Users/39364/Documents/MATLAB/traj_eval/projects/hkisland_gnss01/result_dataset_1/leica_pose.csv"


# follow the x y z order
result = np.loadtxt(input_path,usecols=[0,6,7,8],delimiter=",",skiprows=1)
result[:,[1,2]] = result[:,[1,2]]/180*pi

xyz = blh2xyz_batch(result[:,1:4])
enu = xyz2enu_batch_fixed(xyz,xyz[0,:])

row_size = np.shape(result)[0]

output_result = np.array(
    [
        result[:,0], 1 + np.arange(row_size), result[:,0],
        enu[:,0], enu[:,1], enu[:,2],
        np.zeros(row_size), np.zeros(row_size), np.zeros(row_size), np.ones(row_size)
    ]
).T

np.savetxt(output_path, output_result, delimiter = ',', fmt = '%.6f')


# check the pos xyz order
# fig, ax = plt.subplots()
# ax.plot(enu[:,0], enu[:,1], linewidth = 2)
# plt.show()