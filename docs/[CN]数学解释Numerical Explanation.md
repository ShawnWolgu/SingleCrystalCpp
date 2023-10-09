# 晶体塑性本构模型框架

# 晶体的变形

变形梯度张量 $\mathbf{F}$ 与速度梯度 $\bf{l}$ ：

$$
\bf{F}=\frac{\partial\bf{{x}}}{\partial\bf{{X}}}
$$

$$
\bf{l}=\frac{\partial \bf{v}}{\partial \bf{X}}
$$

两者之间的关系为：

$$
\dot{\mathbf{F}}=\mathbf{l}\mathbf{F}\\\mathbf{F}=\mathbf{F}^e\mathbf{F}^p
$$

$$
\mathbf{l}=\dot{\mathbf{F}}\mathbf{F}^{-1}=\dot{\mathbf{F}^e}(\mathbf{F}^e)^{-1}+\mathbf{F}^e\dot{\mathbf{F}^p}(\mathbf{F}^p)^{-1}(\mathbf{F}^e)^{-1}\\\mathbf{l}^e=\dot{\mathbf{F}^e}(\mathbf{F}^e)^{-1}\\\mathbf{l}^p=\mathbf{F}^e\dot{\mathbf{F}^p}(\mathbf{F}^p)^{-1}(\mathbf{F}^e)^{-1}
$$

对每一个速度梯度，都有其对称部分（应变率张量）与反对称部分（旋率张量）：

$$
\mathbf{l}=\mathbf{d}+\mathbf{w}\\
\mathbf{l}^e=\mathbf{d}^e+\mathbf{w}^e\\
\mathbf{l}^p=\mathbf{d}^p+\mathbf{w}^p
$$

在每一步求解时，晶粒的变形相关张量（应变，取向矩阵）更新：

$$
\mathbf{\varepsilon}_{t+\Delta t}=\mathbf{\varepsilon}_t+\mathbf{d}\Delta t \\ 
\mathbf{M}_{t+\Delta t}=\mathbf{M}_{t}(\mathbf{I}+\mathbf{w}^e\Delta t)^T
$$

这里晶粒的取向矩阵为

$$
\mathbf{M}=\left(\begin{matrix}u &m&h\\v&n&k\\w&o&l\end{matrix}\right)
$$

代表参考坐标系的x，y，z三轴分别对应晶粒坐标系下的$(uvw),(mno),(hkl)$

对于弹性变形梯度，其polar decomposition得到旋转矩阵与变形矩阵：

$$
\mathbf{F}^e=\mathbf{RU}
$$

则有：

$$
\mathbf{M}=\mathbf{M}_0\mathbf{R}^T
$$

# 弹塑性本构关系

$$
\dot{\mathbf{\sigma}}-\mathbf{w}^e\mathbf{\sigma}+\mathbf{\sigma}\mathbf{w}^e+\mathbf{\sigma} tr(\mathbf{d}^e)=\mathbb{C}:\mathbf{d}^e
$$

$$
\dot\sigma=\mathbf{C}(d-d_p)-\sigma tr(d)+w_e\sigma-\sigma w_e\\=\mathbf{C}d-\sigma tr(d) + w\sigma - \sigma w - \mathbf{C}d_p - w_p\sigma +  \sigma w_p
$$

### 收敛求解：牛顿迭代法

$$
\begin{align*}f(X)&=(\mathbf{C}d-\sigma tr(d) + \mathbf{\Sigma} w - \mathbf{C}d_p -   \mathbf{\Sigma} w_p - \dot{\sigma})dt\\&=(\mathbf{C}d-\sigma tr(d) + \mathbf{\Sigma} w - \mathbf{C}d_p -   \mathbf{\Sigma} w_p)dt - \Delta\sigma\end{align*}\\X=(w,d,\Delta\sigma)
$$

式中：

$$
\mathbf{\Sigma}=\begin{pmatrix}0&\sigma_5&\sigma_6\\\sigma_4&0&-\sigma_6\\-\sigma_4&-\sigma_5&0\\\frac12(\sigma_3-\sigma_2)&-\frac12\sigma_6&-\frac12\sigma_5\\-\frac12\sigma_6&\frac12(\sigma_3-\sigma_1)&\frac12\sigma_4\\\frac12\sigma_5&\frac12\sigma_4&\frac12(\sigma_2-\sigma_1
)\end{pmatrix}
$$

则梯度：

$$
\frac{\partial f(X)}{\partial X} = \mathbf{C'} \partial d+\mathbf{\Sigma}\partial w +\frac{\partial f(X)}{\partial{\Delta\sigma}}\partial{\Delta\sigma}
$$

$$
\mathbf{C'}=\mathbf{C}-\begin{pmatrix}\sigma&\sigma&\sigma&0&0&0\end{pmatrix}
$$

为了表述与code对应，将张量写为6维形式：

$$
\frac{\partial f(X)}{\partial \Delta \sigma}=-I-C\frac{\partial\tilde{d_p}}{\partial \Delta\sigma}-\Sigma N^{-1}\frac{\partial\tilde{w_p^*}}{\partial \Delta\sigma}\\=-I-C(M_{cpl}\frac{\partial\tilde{d_p}^c}{\partial \Delta\sigma^c}M_{cpl}^T)-\Sigma N^{-1}(M_{cpl}\frac{\partial\tilde{d_p^{*c}}}{\partial \Delta\sigma^c}M_{cpl}^T)_{4:6}
$$

式中：

$$
\tilde{d_p} = [\begin{matrix}d_1&d_2&d_3&2d_{23}&2d_{13}&2d
_{12}\end{matrix}] \\\tilde{w_p^*}=[\begin{matrix}0&0&0&2w_{23}&2w_{13}&2w
_{12}\end{matrix}]
$$

c上标代表在晶格坐标系，$M_{cpl}$为对应的6阶柔度张量旋转阵。

$$
S = [\begin{matrix}                                                                        s_1n_1&s_2n_2&s_3n_3&(s_2n_3+s_3n_2)&(s_1n_3+s_3n_1)&(s_1n_2+s_2n_1)\end{matrix}]^T\\A = [\begin{matrix}                                                                        0&0&0&(s_2n_3-s_3n_2)&(s_1n_3-s_3n_1)&(s_1n_2-s_2n_1)\end{matrix}]^T\\\frac{\partial{d_p^c}}{\partial \Delta\sigma^c}=SS^T\frac{\partial\dot\gamma dt}{\partial\tau}\\\frac{\partial{w_p^c}}{\partial \Delta\sigma^c}=AS^T \frac{\partial\dot\gamma dt}{\partial\tau}
$$

对于不同的边界条件设定，有

$$
\begin{pmatrix}d\\w\\\Delta\sigma\end{pmatrix}_{15*1}=\begin{pmatrix}I_{3*3}&0&0&0_{3*6}\\0&0.5I_{3*3}&0.5I_{3*3}&0_{3*6}\\0&0.5I_{3*3}&-0.5I_{3*3}&0_{3*6}\\0&0&0&I_{6*6}\end{pmatrix}_{15*15}\begin{pmatrix}L\\\Delta\sigma\end{pmatrix}_{15*1}
$$

得到目标函数$f(X)=0$与梯度$\partial f(X)/\partial X$之后，即可根据牛顿迭代法求$X$，即

$$
X_{n+1}= X_n + \frac{f(X)}{\partial f(X)/\partial X}
$$

### 未收敛的情况：牛顿下山法

$$
X_{n+1}= X_n + c\cdot\frac{f(X)}{\partial f(X)/\partial X}
$$

- 若$|f_{n+1}|/|f_n|>1$，令$c=1/2c$，并重算本步；若$|Y_{n+1}|/|Y_n|\le1$，令$c=1$；

# 塑性变形机制

## 滑移

塑性速度梯度的计算

$$
\mathbf{l}^p=\Sigma_\alpha\dot\gamma^\alpha \mathbf{s}^\alpha \mathbf{n}^{\alpha T} \\
\mathbf{s}^\alpha=\mathbf{F}^e\mathbf{s}_0^\alpha\\
\mathbf{n}^\alpha=(\mathbf{F}^{e-1})^{T}\mathbf{n}_0^\alpha\\
$$

位错速度与分切应力

$$
v=L/(t_w + t_r)
$$

目前对于式中的等待时间$t_w$与运动时间$t_r$的形式并没有绝对的统一形式，这里考虑到两种力的作用：强钉扎力与弱钉扎力。强钉扎力不受温度影响，弱钉扎力受到温度影响。则

$$
\begin{align*}&t_w=\frac{1}{\nu_0}\exp(\frac{Q_a}{k_BT}),\\&Q_a=Q_0[1-\text{sgn}(|\tau|-\tau_f)(\frac{||\tau|-\tau_f|}{\tau_c})^\xi]\\&t_r=L \left(v_{s}\left[\sqrt{1+\left(\frac{v_{s}}{v_m}\right)^{2}}-{\frac{v_s}{v_m}}\right]\right)^{-1},\\&
v_m=\frac{2b}{B_0}|\tau-\tau_f|，
B_{0}=\frac{c_{d} K_{\mathrm{B}} T}{v_{s} b^{\alpha^{2}}}\\\end{align*}
$$

式中$\tau_f$为强钉扎力，而$\tau_c$为弱钉扎力。对于不同的材料体系，两者的组成是不同的，例如一些FCC结构中可以考虑Peierls Barrier为弱作用，而dislocation junction为强作用。上市对应的梯度形式为：

$$
\frac{\partial t_w}{\partial\tau}=-\frac{1}{\nu_0kT}\exp(\frac{Q_a}{kT})\frac{\partial Q_a}{\partial\tau};\frac{\partial Q_a}{\partial\tau}=-\text{sgn}(\tau)\xi Q_0(\frac{||\tau|-\tau_c|}{\tau_o})^{\xi-1}\\\frac{\partial t_r}{\partial \tau}=-\text{sgn}(\tau)\frac{2bL}{v_m^2B_0}\frac{v_s^2}{v_m^2}(1-\frac{v_s}{v_m\sqrt{1+\left(\frac{v_s}{v_m}\right)^2}})
$$

### 一些尝试过又废弃的形式

下面的形式都是只考虑弱作用或者强作用，而弱形式与强形式未合理对应的形式

$$
t_w=[\nu_D\frac {b^2 l}{\Omega_k}\exp(- Q_s/k_BT)]^{-1}\\ Q_s=(Q_0) [1-(|\tau| / \tau_c)^\alpha],Q_0>2k_BT\\  \Omega_k=\frac{Q_0}{\tau_P},\tau_c=\tau_P+\tau_b
$$

$$
t_{r}=\lambda \left(v_{s}\left[\sqrt{1+\left(\frac{v_{s}}{v_m}\right)^{2}}-{\frac{v_s}{v_m}}\right]\right)^{-1} \\
v_m=\frac{2b}{B_0}(|\tau|-\tau_c)+v_c，
B_{0}=\frac{c_{d} K_{\mathrm{B}} T}{v_{s} b^{\alpha^{2}}}
$$

$$
\frac{\partial t_w}{\partial\tau} = -\frac{sgn(\tau)\alpha Q_0^2}{\nu_D b^2 l \tau_P\tau_ck_BT}\exp(\frac{Q_s}{k_B T})(|\tau|/\tau_c)^{\alpha-1}
$$

$$
\frac{\partial t_r}{\partial\tau}=-sgn(\tau)\frac{2b\lambda}{v_m^2B_0}(\frac{v_s^2}{v_m^2}(1-\frac{v_s}{v_m}(1+(\frac{v_s}{v_m})^2)^{-1/2})
$$

或者：

$$
t_w=[\nu_D\frac {b^2 l}{\Omega_k}\exp(- Q_s/k_BT)]^{-1}\\ Q_s=(Q_0) [1-(\frac{|\tau|-\tau_b} {\tau_P})^\alpha],Q_0>2k_BT\\ \Omega_k=\frac{Q_0}{\tau_P},\tau_c=\tau_P+\tau_b
$$

对应的梯度改变为： 

$$
\frac{\partial t_w}{\partial\tau} = -\frac{\text{sgn}(\tau)\alpha Q_0^2}{\nu_D b^2 l \tau_P^2k_BT}\exp(\frac{Q_s}{k_B T})(\frac{|\tau|-\tau_b}{\tau_P})^{\alpha-1}
$$

# 硬化机制

$**\rho_h$代表了其他滑移系位错对当前滑移系的阻碍作用，表示排他性；$\rho_J$代表了当前滑移系位错水平对交互作用的贡献，允许当前滑移系在位错水平较低时降低交互作用的阻碍性。**

$$
\tau_c=\tau_P+c_bbG\sqrt{\rho_h+\rho_J}\\\rho_J=\sum_{\beta\neq\alpha}{h^{\alpha\beta}\sqrt{\rho^\beta_*\rho^\alpha_*}},\rho_*=\rho-\rho_0\\\rho_h=\sum h^{\alpha\beta}\rho^\beta_e,\rho_e=\rho+min(\rho^\beta,\rho^\gamma),\alpha,\beta,\gamma \text{ are coplanar slips.}
$$

为了判断滑移系属于哪种情况，且避免一开始就输入一个巨大的矩阵，还要与滑移系一一对应，还是通过几何条件判断确定$f^{\alpha\beta}$的分类与取值：

1. 判断位错反应是否可以自发产生（Frank定理）$b^\alpha\cdot b^\beta\ne0$ ：是，继续，否则 H
2. 判断伯氏矢量是否平行 $|\cos<b^\alpha,b^\beta>|=1$ ： 是， N， 否则继续
3. 判断是否共面 $|\cos<n^\alpha,n^\beta>|=1$ ： 是： C， 否则继续
4. 判断新生位错是否可以在原滑移系滑移 $n^\alpha\cdot(b^\alpha+b^\beta)=0 || n^\beta\cdot(b^\alpha+b^\beta)=0$：是 G， 否 S

对于上面五种交互模式的说明：

1. N，No junction，两个滑移系上位错反应生成的新位错伯氏矢量与原滑移系位错平行
2. H，Hirth lock，生成的新位错伯氏矢量不符合Frank定理（新位错能量比反应前更低）
3. C，Coplanar junction，生成的新位错伯氏矢量与原有位错共面
4. G，Glissile junction，生成的新位错伯氏矢量符合能量定理，且可以在原滑移面上滑移
5. S，Sessile junction，生成的新位错伯氏矢量符合能量定理，但是不在原位错所在的两个滑移面上，无法滑移。

上面的五种情况下，取协同硬化参数修正系数为$a_N,a_{HL},a_{CJ},a_{GJ},a_{SJ}$，则有：$a_{SJ}\gt a_{GJ} \gt a_{CJ} \gt a_{HL} \gt a_{N}$。

### 饱和位错密度相关

材料中饱和位错密度是率相关的，根据[(Beyerlein & Tomé, 2008)](https://doi.org/10.1016/j.ijplas.2007.07.017) 的工作，记为：

$$
\frac1{\sqrt{\rho_{sat}}}=\frac{k_2}{k_1}=\frac{c_h b}{g}(1-\frac{kT}{Db^3}\ln(\frac{|\dot\varepsilon|}{\dot\varepsilon_0}))
$$

但这一形式还有细节需要讨论。具体位错密度的更新方法也众说纷纭，本code采用以下形式：

$$
d \rho=\left(k_{\rho, n u c} \frac{\left|\tau-\tau_{c, n u c}\right|}{G b^2}+\frac{k_{\rho, m u l}}{\bar{L}}\right)\left[1-\left(\frac{c_h b}{g}\left(1-\frac{k T}{D b^3} \log \frac{|\dot{\varepsilon}|}{\dot{\varepsilon}_0}\right)\right)^2 \rho\right] d \gamma
$$

上式中包含位错非均匀形核项$k_{\rho, n u c} \frac{\left|\tau-\tau_{c, n u c}\right|}{G b^2}$，增殖项$\frac{k_{\rho, m u l}}{\bar{L}}$，并控制最终的饱和位错密度水平。
