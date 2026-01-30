# MALDI-TOF MS 数据处理平台 - 安装与使用指南

## 📋 系统要求

### 1. Python环境
- Python 3.8 或更高版本
- pip 包管理器

### 2. R环境
- R 4.0 或更高版本
- Rscript 命令行工具

### 3. R包依赖
需要安装以下R包：
```r
install.packages("MALDIquant")
install.packages("MALDIquantForeign")
install.packages("readxl")
```

---

## 🚀 快速开始

### 步骤1: 安装R环境和包

#### Windows系统:
1. 从 https://cran.r-project.org/ 下载并安装R
2. 打开R控制台，运行：
```r
install.packages("MALDIquant")
install.packages("MALDIquantForeign")
install.packages("readxl")
```

#### macOS系统:
```bash
# 使用Homebrew安装R
brew install r

# 打开R并安装包
R
> install.packages("MALDIquant")
> install.packages("MALDIquantForeign")
> install.packages("readxl")
```

#### Linux系统:
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install r-base

# 安装R包
R -e "install.packages(c('MALDIquant', 'MALDIquantForeign', 'readxl'), repos='https://cran.rstudio.com/')"
```

### 步骤2: 安装Python依赖

```bash
# 创建虚拟环境（推荐）
python -m venv venv

# 激活虚拟环境
# Windows:
venv\Scripts\activate
# macOS/Linux:
source venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

### 步骤3: 启动应用

```bash
streamlit run maldi_app.py
```

应用将在浏览器中自动打开（默认地址：http://localhost:8501）

---

## 📖 使用说明

### 1. 准备数据

#### 训练集：
- **TXT文件**: MALDI-TOF质谱原始数据文件（多个）
- **Excel文件**: 样本标签文件，必须包含以下两列：
  - `file`: TXT文件名（包含扩展名，如 "sample1.txt"）
  - `group`: 样本分组标签（如 "Resistant", "Sensitive"）

**Excel示例**:
```
file            group
sample1.txt     Resistant
sample2.txt     Resistant
sample3.txt     Sensitive
sample4.txt     Sensitive
```

#### 验证集（可选）：
- 仅需TXT文件，无需Excel标签

### 2. 操作流程

#### 步骤1: 上传数据
1. 切换到「数据上传」标签页
2. 左侧上传训练集TXT文件（可多选）
3. 左侧上传Excel标签文件
4. 右侧上传验证集TXT文件（可选）

#### 步骤2: 配置参数
在侧边栏选择：
- **自动参数估计**（推荐）: 系统根据数据特征自动选择最佳参数
- **手动设置**: 
  - 半峰宽 (halfWindowSize): 10-200，默认90
  - 信噪比阈值 (SNR): 1-10，默认2
  - 对齐容差 (tolerance): 0.001-0.02，默认0.008

#### 步骤3: 开始处理
1. 切换到「处理与结果」标签页
2. 点击「开始处理」按钮
3. 等待处理完成（通常1-5分钟）

#### 步骤4: 查看结果
- 查看处理摘要（样本数、峰数等）
- 查看使用的处理参数
- 下载CSV结果文件

#### 步骤5: 数据可视化
切换到「数据可视化」标签页：
- 查看平均质谱图
- 查看峰强度热图
- 预览数据表格

---

## 📊 输出文件说明

### peak_intensity_train.csv
训练集峰强度矩阵：
- 行：样本组（平均谱）
- 列：m/z特征（格式：mz_1234）
- 值：峰强度

### peak_intensity_validation.csv
验证集峰强度矩阵：
- 行：单个样本
- 列：m/z特征（与训练集相同）
- 值：峰强度

### 使用结果
可将CSV文件导入到：
- Python/R进行机器学习分析
- Excel/SPSS进行统计分析
- 其他数据分析软件

---

## 🔧 常见问题

### Q1: 提示"未找到R环境"
**解决方案**:
1. 确认已安装R并配置环境变量
2. 在命令行测试：`Rscript --version`
3. 如果失败，手动添加R到PATH环境变量

### Q2: R包安装失败
**解决方案**:
```r
# 更换CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("MALDIquant")
```

### Q3: 处理时间过长
**可能原因**:
- 文件过大（>100个样本）
- 参数设置不当（halfWindowSize过大）
- 系统性能限制

**解决方案**:
- 分批处理
- 使用自动参数估计
- 增加超时时间（修改代码中的timeout参数）

### Q4: Excel文件列名不匹配
**解决方案**:
确保Excel文件包含以下两列（大小写敏感）:
- `file`
- `group`

### Q5: TXT文件格式不兼容
**支持格式**:
- 两列数据：m/z 和 intensity
- 制表符或空格分隔
- 可包含表头

---

## 🎯 参数调优建议

### 自动参数 vs 手动参数

| 场景 | 推荐方式 |
|------|---------|
| 标准MALDI-TOF数据 | 自动参数 ✅ |
| 数据质量较差 | 手动调整SNR ⬇️ |
| 峰过多/过少 | 手动调整halfWindowSize |
| 对齐效果差 | 手动增加tolerance |

### 参数说明

**halfWindowSize (半峰宽)**:
- 作用: 控制峰检测的灵敏度
- 数值小 → 检测更多小峰
- 数值大 → 只检测主要峰
- 典型值: 50-150

**SNR (信噪比阈值)**:
- 作用: 过滤噪声峰
- 数值小 → 保留更多峰（包括噪声）
- 数值大 → 只保留强峰
- 典型值: 2-5

**tolerance (对齐容差)**:
- 作用: m/z对齐的宽容度
- 数值小 → 严格对齐
- 数值大 → 宽松对齐
- 典型值: 0.002-0.01

---

## 📝 更新日志

### v1.0.0 (2024)
- ✅ 基础数据处理功能
- ✅ 自动参数估计
- ✅ 可视化界面
- ✅ 批量处理支持
- ✅ 训练集/验证集分离处理

---

## 📧 技术支持

如有问题，请：
1. 检查本文档的「常见问题」部分
2. 确认R和Python环境配置正确
3. 查看处理日志中的错误信息

---

## 📄 许可证

本工具基于开源软件构建：
- MALDIquant: GPL-3
- Streamlit: Apache 2.0
- Plotly: MIT

---

**祝您使用愉快！🎉**
