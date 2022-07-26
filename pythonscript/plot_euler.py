import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def Euler_trans(euler_vector, axis = 'z'):
    euler_vector = euler_vector / 180 * np.pi

    SPH=np.sin(euler_vector[0]); CPH=np.cos(euler_vector[0])
    STH=np.sin(euler_vector[1]); CTH=np.cos(euler_vector[1])
    STM=np.sin(euler_vector[2]); CTM=np.cos(euler_vector[2])

    if axis=='x':
        return np.array([CTM*CPH-SPH*STM*CTH,-STM*CPH-SPH*CTM*CTH,SPH*STH])
    if axis=='y':
        return np.array([CTM*SPH+CPH*STM*CTH,-SPH*STM+CPH*CTM*CTH,-STH*CPH])
    if axis=='z':
        return np.array([STH*STM,CTM*STH,CTH])
    else: return np.array([0,0,1])

def trans_to_xy(axis_vec):
    axis_temp = axis_vec - np.array([0,0,-1])
    axis_temp = axis_temp / abs(axis_temp[2])
    return (axis_temp - np.array([0,0,1]))[:2]

def calc_ipf(axis_vec):
    axis_vec = axis_vec/np.linalg.norm(axis_vec)
    sorted_axis = np.sort(np.abs(axis_vec),)
    axis_ipf = np.array([sorted_axis[1],sorted_axis[0],sorted_axis[2]])
    return axis_ipf

def ipf_curve_points():
    axis = np.hstack((np.ones((31,1)),np.arange(0,1.001+1/30,1/30)[:31].reshape(31,1),np.ones((31,1))))
    pts = []
    for iaxis in axis:
        pts.append(trans_to_xy(calc_ipf(iaxis)))
    return np.array(pts)

path = os.getcwd() + "/euler_angle_grain.csv"
euler_dataf = pd.read_csv(path)
interval = len(euler_dataf)/100
slice = np.arange(0,len(euler_dataf),interval)
euler_dataf = euler_dataf.iloc[slice]
euler_dataf.reset_index()
dots = []
ipfdots = []
for ieuler in range(len(euler_dataf)):
    iaxis = Euler_trans(euler_dataf.iloc[ieuler])
    if ieuler == 0 or ieuler == len(euler_dataf)-1:
        print("Axis:\n")
        print(iaxis)
    ixy = trans_to_xy(iaxis)
    irexy = trans_to_xy(calc_ipf(iaxis))
    dots.append(ixy)
    ipfdots.append(irexy)
dots = np.array(dots)
ipfdots = np.array(ipfdots)

circlep = np.arange(0,2.2*np.pi,0.1*np.pi)
circlex = np.sin(circlep)
circley = np.cos(circlep)

pl, axes = plt.subplots(figsize=(10,10))
axes.plot(circlex,circley,'b-')
axes.plot(dots[:,0],dots[:,1],'r.')
axes.plot(dots[0,0],dots[0,1],'ro',markersize=15)
axes.plot(dots[-1,0],dots[-1,1],'rx',markersize=15)
pl.savefig(os.getcwd() + "/z_axis.png")
plt.close(pl)

pl, axes = plt.subplots(figsize=(10,7.07))
ipfcurve = ipf_curve_points()
axes.plot(ipfcurve[:,0],ipfcurve[:,1],'b-')
axes.plot([0,1.4141/3.4141],[0,0],'b-')
axes.plot([0,1/2.732],[0,1/2.732],'b-')

axes.plot(ipfdots[:,0],ipfdots[:,1],'r.')
axes.plot(ipfdots[0,0],ipfdots[0,1],'ro',markersize=15)
axes.plot(ipfdots[-1,0],ipfdots[-1,1],'rx',markersize=15)
pl.savefig(os.getcwd() + "/z_axis_rev.png")
plt.close(pl)
