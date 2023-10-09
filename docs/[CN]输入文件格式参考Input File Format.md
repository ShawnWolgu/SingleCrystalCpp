# 输入文件格式参考

# 单晶参数文件

### Voce Hardening Law

```
# Elastic Modulus (MPa)
  1.683e5   1.221e5   1.221e5     0.0     0.0     0.0
  1.221e5   1.683e5   1.221e5     0.0     0.0     0.0
  1.221e5   1.221e5   1.683e5     0.0     0.0     0.0
    0.0     0.0     0.0     0.757e5   0.0e5     0.0e5
    0.0     0.0     0.0     0.0e5     0.757e5   0.0e5
    0.0     0.0     0.0     0.0e5     0.0e5     0.757e5
# Hardening Law (0: Voce Hardening, 1: Dislocation Velocity Model)
0
# Lattice Vectors (Angstrom)
3.61  0      0
0     3.61   0
0     0      3.61
# Euler
0 0 0
# Slip System: FCC{111}<1-10>
12
 1  1  1        0  0.5 -0.5    
 1  1  1        0.5  0 -0.5    
 1  1  1        0.5 -0.5  0    
-1  1  1        0  0.5 -0.5    
-1  1  1        0.5  0  0.5    
-1  1  1        0.5  0.5  0    
 1 -1  1        0  0.5  0.5    
 1 -1  1        0.5  0 -0.5    
 1 -1  1        0.5  0.5  0   
 1  1 -1        0  0.5  0.5    
 1  1 -1        0.5  0  0.5   1.0 
 1  1 -1        0.5 -0.5  0    
 40 12 70 20
#t0, t1, h0, h1
 1 1 1 1 1
#N HL CJ GJ SJ
```

### Dislocation Velocity Model

```
# Elastic Modulus (MPa)
  1.683e5   1.221e5   1.221e5     0.0     0.0     0.0
  1.221e5   1.683e5   1.221e5     0.0     0.0     0.0
  1.221e5   1.221e5   1.683e5     0.0     0.0     0.0
    0.0     0.0     0.0     0.757e5   0.0e5     0.0e5
    0.0     0.0     0.0     0.0e5     0.757e5   0.0e5
    0.0     0.0     0.0     0.0e5     0.0e5     0.757e5
# Hardening Law (0: Voce Hardening, 1: Dislocation Velocity Model)
1
# Lattice Vectors (Angstrom)
3.61  0      0
0     3.61   0
0     0      3.61
# Euler
0 0 0
# Slip System: FCC{111}<1-10>
12
 1  1  1        0  0.5 -0.5    
 1  1  1        0.5  0 -0.5    
 1  1  1        0.5 -0.5  0    
-1  1  1        0  0.5 -0.5    
-1  1  1        0.5  0  0.5    
-1  1  1        0.5  0.5  0    
 1 -1  1        0  0.5  0.5    
 1 -1  1        0.5  0 -0.5    
 1 -1  1        0.5  0.5  0   
 1  1 -1        0  0.5  0.5    
 1  1 -1        0.5  0  0.5    
 1  1 -1        0.5 -0.5  0    
1e10	1.0	1e14	1.2 	0.3	0.80	3200	2.5	0.045 0.2 100.	1.2e9	250.	1e9	17.5
#rho_init, c_l, ref_freq, Q0, tau_c, xi, v_sat, c_d, c_h, c_nuc,tau_nuc, c_multi,D,ref_rate, c/g
6	6	6	13.	16.5
#N HL CJ GJ SJ

#[velocity parameters] 
# 1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
# 6. saturated speed, 7. drag coefficient
#[hardening parameters] 
# 8. forest hardening coefficient
#[DD evolution parameters] 
# 0. initial dislocation density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
# 12. drag stress D, 13. reference strain rate, 14. c/g
```

- 参数说明：共15个参数，按照顺序依次为：
    
    [dislocation model] 
    0. initial dislocation density, 
    [velocity parameters] 
    1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent, 6. saturated speed, 7. drag coefficient
    [hardening parameters] 
    8. forest hardening coefficient
    [DD evolution parameters] 
    9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient, 12. drag stress D, 13. reference strain rate, 14. c/g
    
# 边界条件设置文件与例子

## 文件范例

```
1e-8	1	0.2	298

-5	0.0	0.0	
0.0	-5	0.0
0.0	0.0	10	

0	0	0	
1	0	0	
1	1	1	

0	0	0	
0	0	0	
0	0	0	

1	1	1	
0	1	1	
0	0	0
```

- 第一行
    - timestep: 默认时间步长
    - initial_scale: 初始时间步长比例$(\le1)$
    - max_strain: 最大应变（种植条件）
    - temperature: 温度（K）
- 矩阵1～2：速度梯度阵，速度梯度控制阵
    - 矩阵2中元素为1，代表矩阵1中对应元素的速度梯度分量已知
    - 矩阵2中元素为0，代表矩阵1中对应元素的速度梯度分量未知
    - 矩阵2中元素为d，代表矩阵1中对应元素的应变率分量已知
    - 矩阵2中元素为w，代表矩阵1中对应元素的旋率分量已知
- 矩阵3～4：应力增量阵，应力增量控制阵
    - 矩阵4中元素为1，代表矩阵3中对应元素的应力增量已知
    - 注意在当前版本中为应力增量而非应力变化率
    - 要求所有边界条件限制条件总数=9

## 速度梯度控制阵说明

原有的iuflag仅仅使用0 1代表对应元素是否被设定，现在想要改为：

- 1代表设定对应元素已知
- d代表设定对应非对角元素的对称分量已知
- w代表设定对应非对角元素的非对称分量已知

这样可以确定对于一个特定的L速率梯度下，对于某一非对角元素，只可能选择0，d，w, 1 四种情况，其中：

- 如果某元素设定为d，那么对应的转置元素如果是1、w则会自动计算为对应元素的L取值，如果为d且给定数值不一致则直接报错退出；
- 如果某元素设定为w，那么对应的转置元素如果是1、d则会自动计算为对应元素的L取值，如果为w且给定数值不为相反数则直接报错退出；
- 如果某元素设定为0，处理方法按照之前不变；
- 如果某元素设定为1，则参考上述情况处理；
- 最终判断总限定条件个数：对于对角元素只能设置1，0，d，设置w会报错，d视为1，每一个1和d对应一个限定条件；对于非对角元素，与其转置元素共同考虑，即：
    - （0，0）限定条件+0
    - （d，d），（w，w），（1，0），（d，0），（w，0）限定条件+1
    - （d，w），（1，1），（d，1），（w，1） 限定条件+2
    
    总限定条件（包括应力设定）需要为9
    
- 基本原则是尽量转换为对L分量的控制，也就是对于上述的限定条件+2的若干情况，会直接修改边界条件为L分量控制，然后按照之前的情况处理。但是对于上述的非设置10而仅提供限定条件+1的情况，处理方式如下：
    - 记录设定的d或者w分量
    - 在求解中使用的L-dw矩阵中，将对应的转换行设置为仅对称元素为1。举例而言：
    
对于原先的转换矩阵，为：
$$
\begin{pmatrix} d \\\\ w \\\\ \Delta \sigma \end{pmatrix}\_{15\times 1} 
= \begin{pmatrix} I\_{3\times 3} & 0 & 0 & 0\_{3\times 6}\\\\ 
0 & 0.5I\_{3\times 3}&0.5I\_{3\times 3}&0\_{3\times 6}\\\\ 
0 & 0.5I\_{3\times 3}&-0.5I\_{3\times 3}&0\_{3\times 6}\\\\ 
0&0&0&I\_{6\times 6}\end{pmatrix}\_{15\times 15}\begin{pmatrix}L\\\\ \Delta\sigma\end{pmatrix}\_{15\times 1}
$$
如果仅仅设置了d或者w的分量，那么右边列向量就可以写成L和d，w分量的形式，而对应的中间转换矩阵的对应行仅保留对称元素为1
    

### 边界条件控制实例

1. $\dot\varepsilon_{33} = 10/s, L_{21}=0,L_{31}=L_{32}=0$ ：沿3轴拉伸，3面法向不变
    
    ```
    0.0	0.0	0.0	
    0.0	0.0	0.0
    0.0	0.0	10	
    
    0	0	0	
    1	0	0	
    1	1	1	
    
    0	0	0	
    0	0	0	
    0	0	0	
    
    1	1	1	
    0	1	1	
    0	0	0
    ```
    
2. $\dot\varepsilon_{33} = 10/s, L_{13}=L_{23}=0,w_{12}=0$：沿3轴拉伸，3轴方向不变
    
    ```
    0.0	0.0	0.0	
    0.0	0.0	0.0
    0.0	0.0	10	
    
    0	w	1	
    0	0	1	
    0	0	1	
    
    0	0	0	
    0	0	0	
    0	0	0	
    
    1	1	1	
    0	1	1	
    0	0	0
    ```
    
# 配置文件

做了最基础的配置文件设计：指定单晶数据文件以及加载文件的位置。

```
# Single Crystal File
SingleX.txt
# Load File
Load.txt
```

主要是为了可以多个case调用同一个文件。
但是更方便的用法为直接命令行指定输入文件：
```
SXCpp -x path/to/crystalfile -l path/to/loadfile
```
