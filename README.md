# SingleCrystalCpp
用C++写了一个简单的单晶晶体塑性模型

## 目前实现的功能
* 支持混合边界条件
* Voce硬化模型
* 位错速度模型

## 编译
* 下载cmake
* 下载Eigen3 tclap，准备CMakeLists.txt
* 用cmake编译所有文件
* 完成

## 输入文件
* 滑移系以及塑性模型参数
* 加载与边界条件

## 使用说明
[输入参数格式](/docs/%5BCN%5D%E8%BE%93%E5%85%A5%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E5%8F%82%E8%80%83Input%20File%20Format.md)

[模型数学解释](/docs/%5BCN%5D%E6%95%B0%E5%AD%A6%E8%A7%A3%E9%87%8ANumerical%20Explanation.md)

## 引用
Sun X, Zhou K, Liu C, et al. A Crystal Plasticity Based Strain Rate Dependent Model across An Ultra-wide Range[J]. International Journal of Plasticity, 2024: 104056.

```bibtex
@article{sun2024crystal,
  title={A Crystal Plasticity Based Strain Rate Dependent Model across An Ultra-wide Range},
  author={Sun, Xiaochuan and Zhou, Kecheng and Liu, Chuhao and others},
  journal={International Journal of Plasticity},
  pages={104056},
  year={2024},
  publisher={Elsevier}
}
```

# SingleCrystalCpp
A simple crystal plasticity model for single crystal, developed in C++.

## Available Features
* Hybrid boundary conditions
* Voce hardening
* Dislocation velocity model for hardening

## Compiling
* Download cmake if not yet
* Download Eigen3 and tclap, prepare your CMakeLists.txt file
* Compile the files using cmake
* Finish

## Input Files
* Slip systems and the model parameters
* Loading and boundary conditions

## Usage
[Input File Format (in Chinese)](/docs/%5BCN%5D%E8%BE%93%E5%85%A5%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E5%8F%82%E8%80%83Input%20File%20Format.md) 

[Model Explanation (in Chinese)](/docs/%5BCN%5D%E6%95%B0%E5%AD%A6%E8%A7%A3%E9%87%8ANumerical%20Explanation.md)

## Citation
Sun X, Zhou K, Liu C, et al. A Crystal Plasticity Based Strain Rate Dependent Model across An Ultra-wide Range[J]. International Journal of Plasticity, 2024: 104056.

```bibtex
@article{sun2024crystal,
  title={A Crystal Plasticity Based Strain Rate Dependent Model across An Ultra-wide Range},
  author={Sun, Xiaochuan and Zhou, Kecheng and Liu, Chuhao and others},
  journal={International Journal of Plasticity},
  pages={104056},
  year={2024},
  publisher={Elsevier}
}
```
