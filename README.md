# SU2 网格读取器（ParaView Python 插件）🧩

这是一个轻量的 ParaView Python Reader 插件，用于读取 **SU2 网格文件（`.su2`）** 并进行可视化。  
分享目的很简单：**我能用，也分享给你用**；**不提供支持，不承担任何责任** ⚠️

该 Reader 提供两个输出端口：
- **Port 0：体网格（Volume Mesh）**：`vtkUnstructuredGrid` 
- **Port 1：边界网格（按 Marker 分组）**：`vtkMultiBlockDataSet`（每个 block 对应一个 SU2 marker）

**已测试环境 ✅**
- Ubuntu 24.04 LTS + ParaView 6.0.1  
- Windows 11 + ParaView 5.13.3  

---

## 安装与加载 🔧

1. 打开 ParaView  
2. `Tools → Manage Plugins...`  
3. `Load New...` 选择 `SU2Reader.py`  
4. （建议）勾选 **Auto Load**，以后启动自动加载  

---

## 使用方法 ▶️

1. `File → Open` 选择 `*.su2` 文件  
2. 点击 `Apply`  

打开后你会得到两种输出可视化 👀：
- **Port 0：体网格**（整体网格）
- **Port 1：边界网格（按 marker 分组）**（MultiBlock，可按 marker 单独查看/提取）

---

## 免责声明（重要）⚠️

本项目仅用于分享与参考，按 **“现状（AS IS）”** 提供。  
作者不提供任何形式的支持或担保，也不对使用本项目产生的任何问题、损失或纠纷承担责任。  
使用即表示你理解并接受上述条款。

---

## 开源协议 📄

MIT License（见仓库根目录 `LICENSE` 文件）。
