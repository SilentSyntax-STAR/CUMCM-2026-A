2026年高教社杯全国大学生数学建模竞赛 A题
支撑材料说明（README）

一、材料用途
本支撑材料用于复现论文《基于热湿耦合传递的药材分阶段烘干过程建模与仿真》中问题1—问题4的主要数值结果，并提供 MWORKS.Syslab/Julia 独立复算程序。

正式提交结果文件为：
1. result1.xlsx —— 问题1：0~1800 s 内温度场与水分浓度场，每隔1 s、径向每隔0.1 cm输出；
2. result2.xlsx —— 问题2：0~3 h 内温度场与水分浓度场，每隔1 s、径向每隔0.1 cm输出；
3. result3.xlsx —— 问题3：全过程水分浓度，每隔60 s、径向每隔0.1 cm输出；
4. result4.xlsx —— 问题4：考虑尺寸收缩后的全过程水分浓度，每隔60 s、径向每隔0.1 cm输出，并单列药材表面水分浓度。

所有正式结果均按赛题要求保留四位小数。

二、建议的目录结构
支撑材料/
├─ README.txt
├─ 附件1.xlsx
├─ 附件2.xlsx
├─ result1.xlsx
├─ result2.xlsx
├─ result3.xlsx
├─ result4.xlsx
└─ code/
   ├─ q1_solver.py
   ├─ q2_solver.py
   ├─ q3_solver.py
   ├─ q4_solver.py
   ├─ validation_q1.jl
   ├─ validation_q2.jl
   ├─ validation_q3.jl
   └─ validation_q4.jl

说明：附件1.xlsx、附件2.xlsx为赛题原始附件的副本，仅用于保证程序可独立复现；其内容不属于本队自主采集数据。

三、Python 主程序运行环境
Python 3.x

所需第三方库：
- numpy
- pandas
- scipy
- matplotlib
- openpyxl

可使用以下命令安装依赖：
pip install numpy pandas scipy matplotlib openpyxl

四、Python 主程序运行方法
请先将当前工作目录切换到“支撑材料”根目录，再运行 code 目录中的程序。程序采用相对路径读取附件和结果模板，因此工作目录应保持为支撑材料根目录。

依次运行：
python code/q1_solver.py
python code/q2_solver.py
python code/q3_solver.py
python code/q4_solver.py

各程序与论文问题的对应关系如下。

1. q1_solver.py
数值方法：一维径向有限体积法 + SciPy BDF 隐式积分。
输入：附件1.xlsx。
主要输出：result1_N2561.xlsx，以及正文所用图表文件。
正式提交的 result1.xlsx 与该主计算结果对应，只是按赛题要求统一文件名。

2. q2_solver.py
数值方法：变物性一维径向有限体积法 + SciPy BDF 隐式积分。
输入：附件1.xlsx。
主要输出：result2_out.xlsx，以及正文所用图表文件。
正式提交的 result2.xlsx 与该主计算结果对应，只是按赛题要求统一文件名。

3. q3_solver.py
数值方法：径向有限体积法 + 全隐式时间推进 + Picard 迭代 + 带状线性方程求解。
输入：附件1.xlsx、result3.xlsx（作为赛题规定格式的结果表模板）。
主要输出：result3_final.xlsx、Q3_table5_final.xlsx、Q3_Cmax_threshold_final.png。
正式提交的 result3.xlsx 对应最终全过程水分浓度结果。

4. q4_solver.py
数值方法：随体坐标下的移动边界有限体积法 + 全隐式时间推进 + Picard 迭代；半径采用 PCHIP 插值。
输入：附件1.xlsx、附件2.xlsx、result4.xlsx（作为赛题规定格式的结果表模板）。
主要输出：result4_final.xlsx、Q4_table_6h.xlsx、Q4_Cmax_threshold_final.png、Q4_radius_final.png。
正式提交的 result4.xlsx 对应最终全过程水分浓度结果。

五、关键结果核对
运行结束后，可用以下数值进行快速核对。

问题1：
- t=1800 s 时，中心温度约为 33.5753 ℃，表面温度约为 36.7856 ℃；
- t=1800 s 时，中心水分浓度约为 2.5500 kg/kg，表面水分浓度约为 1.5102 kg/kg。

问题2：
- t=3 h 时，中心温度约为 49.8495 ℃，表面温度约为 49.9664 ℃；
- t=3 h 时，中心水分浓度约为 1.7662 kg/kg，表面水分浓度约为 1.0081 kg/kg。

问题3：
- 连续时间模型的理论临界烘干时长约为 57.4929 h；
- 按60 s结果输出间隔，首个全域严格满足 Cmax<0.15 的记录时刻为 57.5000 h。

问题4：
- 连续时间模型的理论临界烘干时长约为 51.0939 h；
- 按60 s结果输出间隔，首个全域严格满足 Cmax<0.15 的记录时刻为 51.1000 h；
- 收缩后半径最终稳定在约 1.200 cm。

六、MWORKS.Syslab/Julia 独立复算
validation_q1.jl—validation_q4.jl 用于论文中的跨程序独立复算检验。

这些程序面向 MWORKS.Syslab/Julia 环境，使用：
- TyBase
- TyMath
- TyPlot
- Statistics
- Printf
- ode15s

建议同样以“支撑材料”根目录作为当前工作目录，再在 MWORKS.Syslab 中运行：
code/validation_q1.jl
code/validation_q2.jl
code/validation_q3.jl
code/validation_q4.jl

其中 Python 正式结果仅用于 Julia 独立求解完成后的误差比较，不参与 Julia 数值求解过程。

validation_q1.jl、validation_q2.jl 会生成相应的一致性统计文件与对比图；
validation_q3.jl、validation_q4.jl 会分别在脚本目录下建立 Q3_Julia_Output、Q4_Julia_Output 文件夹保存独立复算结果。

七、result4.xlsx 空白单元格说明
问题4中药材半径随时间减小。对于某一时刻，如果固定物理位置（如 r=1.5 cm 或 r=2.0 cm）已经位于收缩后的药材外部，则该单元格保持空白；“药材表面”列单独给出该时刻真实移动表面处的水分浓度。空白单元格不是缺失计算结果，而表示对应空间位置已不属于药材内部。

八、文件关系说明
论文正文中的表格数值均由上述正式结果文件或相应主程序直接提取；图形由同一组数值结果生成。
正式提交时应以 result1.xlsx、result2.xlsx、result3.xlsx、result4.xlsx 为赛题要求的四个结果文件，不应以带有 _final、_out、_N2561 等开发阶段后缀的文件替代正式文件名。

九、复现说明
由于不同软件版本、线性代数库及刚性积分器内部容差控制可能产生极小的末位差异，复算时应重点核对论文给出的四位小数结果、理论临界时长及收敛/一致性量级，而不要求中间浮点数逐位完全一致。
