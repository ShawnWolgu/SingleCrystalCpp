# SingleCrystalCpp
用C++写了一个简单的单晶晶体塑性模型

## 目前实现的功能
* 支持混合边界条件
* Voce硬化模型
* 位错硬化模型
* 位错速度模型

## 编译
* 下载cmake
* 下载Eigen3，改变CMakeLists.txt中Eigen3库的位置
* 用cmake编译所有文件
* 完成

## 输入文件
* 滑移系以及塑性模型参数
* 加载与边界条件


# SingleCrystalCpp
A simple crystal plasticity model for single crystal, developed in C++.

## Available Features
* Hybrid boundary conditions
* Voce hardening
* Dislocation density hardening
* Dislocation velocity model for hardening

## Compiling
* Download cmake if not yet
* Download Eigen3 and change the CMakeLists.txt to modify the address of Eigen3 repo
* Compile the files using cmake
* Finish

## Input Files
* Slip systems and the model parameters
* Voce Hardening
* Loading and boundary conditions
