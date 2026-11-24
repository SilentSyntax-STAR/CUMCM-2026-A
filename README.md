# CUMCM 2026 · Problem A

### 🏆 2026 年“高教社杯”全国大学生数学建模竞赛 · A 题

**基于热湿耦合传递的药材分阶段烘干过程**

*Modeling the Staged Drying of Medicinal Materials via Coupled Heat and Moisture Transfer*

**竞赛成绩 / Award：** 省级一等奖

**Tech Stack：** Python · MWORKS.Syslab / Julia · FVM · BDF · PCHIP · Moving Boundary

> **从物理机制建模，到跨语言独立复算。**
>
> *From Physics-Based Modeling to Independent Cross-Language Numerical Verification.*

## 🔬 核心技术亮点 | Technical Highlights

### 01 · Python × Julia Independent Verification

**跨语言、跨数值方法的独立复算**

使用 Python 实现基于有限体积法（Finite Volume Method, FVM）的主求解器，并采用 MWORKS.Syslab / Julia，通过二阶有限差分与 MOL–ode15s 构建独立验证程序。

两套程序分别进行空间离散与数值求解，不共享中间计算状态，从而检查求解结果对程序实现方式的依赖程度。

**Cross-Implementation Consistency：**

- Q3：临界烘干时间相对差异 **0.0830%**
- Q4：临界烘干时间相对差异 **0.0393%**

### 02 · Moving Boundary Modeling

**考虑药材尺寸收缩的移动边界模型**

基于实测半径数据，利用 PCHIP 插值与随体坐标映射，将随时间收缩的物理计算域转换为固定计算域。

*Lagrangian coordinate transformation for transient shrinking-boundary simulation.*

### 03 · Robust Numerical Methods

**隐式数值求解与收敛性检验**

采用 BDF 隐式积分、全隐式时间推进与 Picard 迭代处理非线性热湿传递问题，并通过网格加密、时间步长收敛性和独立程序复算评估数值一致性。

### 04 · Sensitivity Analysis

**关键参数敏感性分析**

针对有效水分扩散系数、表面对流传质系数和换热系数开展参数扰动分析，识别长期干燥过程中的主要控制因素。

结果表明，在论文考察的扰动范围内，烘干时长对**有效水分扩散系数**最敏感。

### 05 · Controlled Comparison

**通过控制变量分析分离物性与几何效应**

在相同附录 4 物性关系下，比较固定半径与实测收缩半径模型：

| Model | Predicted Drying Time |
|---|---:|
| Fixed Radius | 129.87 h |
| Shrinking Radius | 51.09 h |
| Relative Reduction | **60.66%** |

该结果为题设参数下的模型预测，而非实际工业试验测得的效率提升。

## 📊 建模路线 | Modeling Pipeline

| Stage | Modeling Objective | Numerical Approach |
|---|---|---|
| Q1 | 预热阶段热湿场预测 | FVM + BDF |
| Q2 | 状态相关物性与热湿耦合 | Variable-Property Model |
| Q3 | 全域干燥终止事件判定 | Implicit Scheme + Picard |
| Q4 | 药材收缩与移动边界修正 | PCHIP + Lagrangian Mapping |

**Baseline Model → Variable Properties → Drying Endpoint → Moving Boundary**

## 📁 项目材料 | Project Materials

- `论文.pdf`：最终竞赛论文
- `A题题目要求/`：赛题及相关材料
- `支撑材料/code/`：Python 主求解器、Julia 独立复算及敏感性分析代码
- `支撑材料/result1.xlsx` 至 `result4.xlsx`：四问计算结果
- `支撑材料/README.txt`：程序运行与结果复现说明

## 📌 项目说明 | Notes

本仓库记录 2026 年全国大学生数学建模竞赛 A 题的建模过程、数值模拟、独立复算及论文成果。

*This repository documents our physics-based modeling, numerical simulations, independent computational verification, and reproducible results for CUMCM 2026 Problem A.*

跨语言复算用于验证数值实现的一致性，不等同于独立实验验证。仓库目前保持私有，待官方成绩公布并确认材料可公开后，再考虑正式发布。
