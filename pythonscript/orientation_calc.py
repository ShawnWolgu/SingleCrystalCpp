import numpy as np
import os, sys

def to_euler(ori_matrix):
    TH = np.arccos(ori_matrix[2,2])
    if abs(np.cos(TH)) >= 0.99:
        TM = 0
        PH = np.arctan2(ori_matrix[0,1],ori_matrix[0,0])
    else:
        STH = np.sin(TH) + 1e-20
        TM = np.arctan2(ori_matrix[0,2]/STH,ori_matrix[1,2]/STH)
        PH = np.arctan2(ori_matrix[2,0]/STH,-ori_matrix[2,1]/STH)
    euler = np.array([PH,TH,TM]) * 180 / np.pi
    return euler

def to_matrix(euler):
    euler = euler / 180 * np.pi
    PH = euler[0]
    TH = euler[1]
    TM = euler[2]

    Mout = np.identity(3)
    SPH=np.sin(PH)
    CPH=np.cos(PH)
    STH=np.sin(TH)
    CTH=np.cos(TH)
    STM=np.sin(TM)
    CTM=np.cos(TM)
    Mout[0,0]=CTM*CPH-SPH*STM*CTH
    Mout[1,0]=-STM*CPH-SPH*CTM*CTH
    Mout[2,0]=SPH*STH
    Mout[0,1]=CTM*SPH+CPH*STM*CTH
    Mout[1,1]=-SPH*STM+CPH*CTM*CTH
    Mout[2,1]=-STH*CPH
    Mout[0,2]=STH*STM
    Mout[1,2]=CTM*STH
    Mout[2,2]=CTH

    return Mout;

def give_matrix(axis):
    axis = axis/np.linalg.norm(axis)
    axis_y = np.array([0,0,1])
    axis_x = np.cross(axis_y,axis)
    axis_y = np.cross(axis,axis_x)
    return np.hstack((axis_x.reshape(3,1),axis_y.reshape(3,1),axis.reshape(3,1)))

def check_xy(new_matrix):
    v_x = np.sort(np.abs(new_matrix[:,0]))
    v_y = np.sort(np.abs(new_matrix[:,1]))
    if np.linalg.norm(v_y-v_x) < 0.0001:
        return True
    else:
        return False

def print_ori(axis):
    print(axis)
    orient_mat = give_matrix(axis)
    print(orient_mat)
    euler = to_euler(orient_mat)
    print(euler)
    for euler_z in np.arange(0,90,0.01):
        euler[0] = euler_z
        new_matrix = to_matrix(euler)
        if check_xy(new_matrix):
            print("Rotation Matrix:")
            print(new_matrix)
            print("Euler Angle:")
            print(to_euler(new_matrix))
            print("Matrix Check:")
            print(to_matrix(to_euler(new_matrix)))
            break
        else:
            continue

axis = [float(i) for i in sys.argv[1:]]
axis = np.array(axis)
print_ori(axis)

