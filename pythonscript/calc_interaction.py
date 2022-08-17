import numpy as np

slip_sys = []
slip_sys.append([np.array([1,1,1]),np.array([0,1,-1])])
slip_sys.append([np.array([1,1,1]),np.array([1,0,-1])])
slip_sys.append([np.array([1,1,1]),np.array([1,-1,0])])
slip_sys.append([np.array([-1,1,1]),np.array([0,1,-1])])
slip_sys.append([np.array([-1,1,1]),np.array([1,0,1])])
slip_sys.append([np.array([-1,1,1]),np.array([1,1,0])])
slip_sys.append([np.array([-1,-1,1]),np.array([0,1,1])])
slip_sys.append([np.array([-1,-1,1]),np.array([1,0,1])])
slip_sys.append([np.array([-1,-1,1]),np.array([1,-1,0])])
slip_sys.append([np.array([1,-1,1]),np.array([0,1,1])])
slip_sys.append([np.array([1,-1,1]),np.array([1,0,-1])])
slip_sys.append([np.array([1,-1,1]),np.array([1,1,0])])

matrix = np.zeros((12,12))
for i,isys in enumerate(slip_sys):
    for j,jsys in enumerate(slip_sys):
        t_vec = np.cross(jsys[0],jsys[1])
        t_vec = t_vec/np.linalg.norm(t_vec)
        cos_mn = np.dot((isys[0]/np.linalg.norm(isys[0])),t_vec)
        cos_b = np.dot(isys[1],jsys[1])/np.linalg.norm(isys[1])/np.linalg.norm(jsys[1])
        matrix[i,j] = np.sqrt(1-cos_mn**2)#*cos_b

disl = np.ones((12,1))
for jjj in range(30):
    for idisl in [4,5,6,7,9,11]:
        disl[idisl] += (0.5 * np.sqrt(disl[idisl]) - 0.01*disl[idisl]) * 2.5 / (np.matmul(matrix,disl)[idisl]) 
    print(np.transpose(disl))
print(matrix)
