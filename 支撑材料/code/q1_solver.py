import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.sparse import diags


# ============================================================
# 0. 计时
# ============================================================

start_clock = time.time()


# ============================================================
# 1. 读取附件1
# ============================================================

env = pd.read_excel("附件1.xlsx")

time_env = env["时间"].to_numpy(dtype=float)
T_air_data = env["温度"].to_numpy(dtype=float)
C_air_data = env["水分浓度"].to_numpy(dtype=float)


def T_air(t):
    """
    烘房温度：
    对附件1数据进行分段线性插值
    """
    return np.interp(
        t,
        time_env,
        T_air_data
    )


def C_air(t):
    """
    烘房水分浓度：
    对附件1数据进行分段线性插值
    """
    return np.interp(
        t,
        time_env,
        C_air_data
    )


# ============================================================
# 2. Q1基本物理参数
# ============================================================

R = 0.02               # 药材半径 / m

# ------------------------------------------------------------
# 正式计算采用细网格
#
# N = 2561
#
# dr = 0.02 / 2560
#    = 7.8125e-6 m
#    = 0.00078125 cm
# ------------------------------------------------------------

N = 2561

r = np.linspace(
    0.0,
    R,
    N
)

dr = r[1] - r[0]


rho = 820.0            # kg/m^3
cp = 2600.0            # J/(kg·K)
k = 0.36               # W/(m·K)

h = 25.0               # W/(m^2·K)
hm = 8e-7              # m/s

T0 = 28.0              # °C
C0 = 2.55              # kg/kg


print("=" * 70)
print("Q1 最终计算参数")
print("=" * 70)

print(f"N  = {N}")
print(f"dr = {dr:.12e} m")
print(f"dr = {dr * 100:.9f} cm")
print()


# ============================================================
# 3. 水分扩散系数
# ============================================================

def diffusion_coefficient(C):
    """
    Q1有效水分扩散系数

    D(C) = 7e-9 * exp(-0.89/C)

    单位：m^2/s
    """

    C_safe = np.maximum(
        C,
        1e-10
    )

    return (
        7e-9
        *
        np.exp(
            -0.89 / C_safe
        )
    )


# ============================================================
# 4. 调和平均
# ============================================================

def harmonic_mean(a, b):
    """
    相邻控制体界面的扩散系数采用调和平均
    """

    return (
        2.0
        *
        a
        *
        b
        /
        (
            a
            +
            b
            +
            1e-30
        )
    )


# ============================================================
# 5. 有限体积控制体几何
# ============================================================

# 节点控制体左边界
r_left = np.maximum(
    r - dr / 2.0,
    0.0
)

# 节点控制体右边界
r_right = np.minimum(
    r + dr / 2.0,
    R
)

# ------------------------------------------------------------
# 忽略公共的 πL 因子
#
# Vi = 1/2 (r_right^2 - r_left^2)
# ------------------------------------------------------------

V = (
    0.5
    *
    (
        r_right ** 2
        -
        r_left ** 2
    )
)

# ------------------------------------------------------------
# 相邻节点之间的界面半径
#
# 长度 N-1
# ------------------------------------------------------------

r_face = (
    0.5
    *
    (
        r[:-1]
        +
        r[1:]
    )
)


# ============================================================
# 6. BDF Jacobian 稀疏结构
#
# 每个节点只和：
# i-1、i、i+1
# 有关，因此 Jacobian 为三对角结构
# ============================================================

jac_pattern = diags(
    diagonals=[
        np.ones(N - 1),
        np.ones(N),
        np.ones(N - 1)
    ],
    offsets=[
        -1,
        0,
        1
    ],
    shape=(
        N,
        N
    ),
    format="csc"
)


# ============================================================
# 7. 温度方程右端
# ============================================================

def rhs_temperature(t, T):
    """
    rho cp dT/dt
    =
    1/r d/dr (k r dT/dr)
    """

    Ta = T_air(t)

    dTdt = np.empty_like(T)

    # --------------------------------------------------------
    # 内部界面导热通量
    #
    # F = r_face * k * dT/dr
    #
    # 长度 N-1
    # --------------------------------------------------------

    flux_face = (
        r_face
        *
        k
        *
        (
            T[1:]
            -
            T[:-1]
        )
        /
        dr
    )


    # --------------------------------------------------------
    # 中心节点 i=0
    #
    # 左侧中心对称通量 = 0
    # --------------------------------------------------------

    dTdt[0] = (
        flux_face[0]
        /
        (
            rho
            *
            cp
            *
            V[0]
        )
    )


    # --------------------------------------------------------
    # 内部节点
    # --------------------------------------------------------

    dTdt[1:-1] = (
        flux_face[1:]
        -
        flux_face[:-1]
    ) / (
        rho
        *
        cp
        *
        V[1:-1]
    )


    # --------------------------------------------------------
    # 表面节点
    #
    # k dT/dr = h(Ta-Ts)
    #
    # 外表面通量：
    #
    # F_R = R h (Ta - Ts)
    # --------------------------------------------------------

    flux_surface = (
        R
        *
        h
        *
        (
            Ta
            -
            T[-1]
        )
    )


    dTdt[-1] = (
        flux_surface
        -
        flux_face[-1]
    ) / (
        rho
        *
        cp
        *
        V[-1]
    )


    return dTdt


# ============================================================
# 8. 水分方程右端
# ============================================================

def rhs_moisture(t, C):
    """
    dC/dt
    =
    1/r d/dr
    (
        r D(C) dC/dr
    )
    """

    Ca = C_air(t)

    dCdt = np.empty_like(C)


    # --------------------------------------------------------
    # 节点扩散系数
    # --------------------------------------------------------

    D_node = diffusion_coefficient(
        C
    )


    # --------------------------------------------------------
    # 界面扩散系数：调和平均
    # --------------------------------------------------------

    D_face = harmonic_mean(
        D_node[:-1],
        D_node[1:]
    )


    # --------------------------------------------------------
    # 内部界面水分通量
    # --------------------------------------------------------

    flux_face = (
        r_face
        *
        D_face
        *
        (
            C[1:]
            -
            C[:-1]
        )
        /
        dr
    )


    # --------------------------------------------------------
    # 中心节点
    #
    # 中心对称，左侧通量为0
    # --------------------------------------------------------

    dCdt[0] = (
        flux_face[0]
        /
        V[0]
    )


    # --------------------------------------------------------
    # 内部节点
    # --------------------------------------------------------

    dCdt[1:-1] = (
        flux_face[1:]
        -
        flux_face[:-1]
    ) / (
        V[1:-1]
    )


    # --------------------------------------------------------
    # 表面 Robin 边界
    #
    # D dC/dr = hm(Ca-Cs)
    #
    # 外表面通量：
    #
    # F_R = R hm (Ca-Cs)
    # --------------------------------------------------------

    flux_surface = (
        R
        *
        hm
        *
        (
            Ca
            -
            C[-1]
        )
    )


    dCdt[-1] = (
        flux_surface
        -
        flux_face[-1]
    ) / (
        V[-1]
    )


    return dCdt


# ============================================================
# 9. 初值
# ============================================================

T_init = np.full(
    N,
    T0,
    dtype=float
)

C_init = np.full(
    N,
    C0,
    dtype=float
)


# ============================================================
# 10. 输出时间
#
# 0,1,2,...,1800 s
# ============================================================

t_eval = np.arange(
    0.0,
    1801.0,
    1.0
)


# ============================================================
# 11. 单独求温度场
# ============================================================

print("=" * 70)
print("开始求解 Q1 温度场...")
print("=" * 70)

clock_T = time.time()


sol_T = solve_ivp(
    fun=rhs_temperature,

    t_span=(
        0.0,
        1800.0
    ),

    y0=T_init,

    method="BDF",

    t_eval=t_eval,

    rtol=1e-10,
    atol=1e-12,
    max_step=10.0,

    jac_sparsity=jac_pattern
)


if not sol_T.success:

    raise RuntimeError(
        "温度场求解失败："
        +
        sol_T.message
    )


print(
    f"温度场求解完成，耗时 "
    f"{time.time() - clock_T:.3f} s"
)

print()


# ============================================================
# 12. 单独求水分场
# ============================================================

print("=" * 70)
print("开始求解 Q1 水分场...")
print("=" * 70)

clock_C = time.time()


sol_C = solve_ivp(
    fun=rhs_moisture,

    t_span=(
        0.0,
        1800.0
    ),

    y0=C_init,

    method="BDF",

    t_eval=t_eval,

    rtol=1e-10,
    atol=1e-12,
    max_step=10.0,

    jac_sparsity=jac_pattern
)


if not sol_C.success:

    raise RuntimeError(
        "水分场求解失败："
        +
        sol_C.message
    )


print(
    f"水分场求解完成，耗时 "
    f"{time.time() - clock_C:.3f} s"
)

print()


# ============================================================
# 13. 转换结果矩阵
#
# 行 = 时间
# 列 = 空间节点
# ============================================================

T_all = sol_T.y.T
C_all = sol_C.y.T


# ============================================================
# 14. 提取题目要求：
#
# r = 0,0.1,...,2.0 cm
# ============================================================

output_r_cm = np.arange(
    0.0,
    2.0001,
    0.1
)


# ------------------------------------------------------------
# 根据实际坐标寻找最近计算节点
#
# 对N=2561来说：
# 0.1 cm = 128个细网格
# ------------------------------------------------------------

output_r_m = (
    output_r_cm
    /
    100.0
)


indices = np.rint(
    output_r_m
    /
    dr
).astype(int)


# 防止浮点问题
indices = np.clip(
    indices,
    0,
    N - 1
)


# 检查真实节点是否正好对应
actual_output_r_cm = (
    r[indices]
    *
    100.0
)


print("=" * 70)
print("输出空间节点检查")
print("=" * 70)

print(
    np.round(
        actual_output_r_cm,
        10
    )
)

print()


if not np.allclose(
    actual_output_r_cm,
    output_r_cm,
    atol=1e-10
):

    raise RuntimeError(
        "0.1 cm 输出节点未与计算网格严格重合"
    )


T_out = T_all[
    :,
    indices
]

C_out = C_all[
    :,
    indices
]
# ============================================================
# 高精度结果：仅用于 Python-Syslab 独立复算比较
# ============================================================

df_T_validation = pd.DataFrame(
    T_out,
    columns=np.round(output_r_cm, 1)
)

df_C_validation = pd.DataFrame(
    C_out,
    columns=np.round(output_r_cm, 1)
)

df_T_validation.insert(
    0,
    "时间",
    np.arange(0, 1801)
)

df_C_validation.insert(
    0,
    "时间",
    np.arange(0, 1801)
)

with pd.ExcelWriter(
    "result1_validation_full.xlsx",
    engine="openpyxl"
) as writer:

    df_T_validation.to_excel(
        writer,
        sheet_name="温度",
        index=False
    )

    df_C_validation.to_excel(
        writer,
        sheet_name="水分浓度",
        index=False
    )

print("result1_validation_full.xlsx 已生成")

# ============================================================
# 15. 保存Excel
# ============================================================

df_T = pd.DataFrame(
    T_out,
    columns=np.round(
        output_r_cm,
        1
    )
)

df_C = pd.DataFrame(
    C_out,
    columns=np.round(
        output_r_cm,
        1
    )
)


df_T.insert(
    0,
    "时间",
    np.arange(
        0,
        1801
    )
)

df_C.insert(
    0,
    "时间",
    np.arange(
        0,
        1801
    )
)

# ------------------------------------------------------------
# 题目要求保留四位小数
# ------------------------------------------------------------

df_T.iloc[
    :,
    1:
] = (
    df_T.iloc[
        :,
        1:
    ]
    .round(4)
)


df_C.iloc[
    :,
    1:
] = (
    df_C.iloc[
        :,
        1:
    ]
    .round(4)
)


OUTPUT_EXCEL = (
    "result1.xlsx"
)


with pd.ExcelWriter(
    OUTPUT_EXCEL,
    engine="openpyxl"
) as writer:

    df_T.to_excel(
        writer,
        sheet_name="温度",
        index=False
    )

    df_C.to_excel(
        writer,
        sheet_name="水分浓度",
        index=False
    )


print(
    OUTPUT_EXCEL,
    "已生成"
)

print()


# ============================================================
# 16. 打印论文表1、表2
# ============================================================

report_times = np.array(
    [
        100,
        300,
        600,
        900,
        1200,
        1500,
        1800
    ],
    dtype=int
)


report_r_cm = np.array(
    [
        0.0,
        0.5,
        1.0,
        1.5,
        2.0
    ]
)


report_r_m = (
    report_r_cm
    /
    100.0
)


report_indices = np.rint(
    report_r_m
    /
    dr
).astype(int)


# ------------------------------------------------------------
# 因为t_eval每秒一个点，时间值本身就是数组行号
# ------------------------------------------------------------

table_T = pd.DataFrame(
    T_all[
        report_times
    ][
        :,
        report_indices
    ],

    index=report_times,

    columns=report_r_cm
).round(4)


table_C = pd.DataFrame(
    C_all[
        report_times
    ][
        :,
        report_indices
    ],

    index=report_times,

    columns=report_r_cm
).round(4)


print("=" * 70)
print("表1：温度 / °C")
print("=" * 70)

print(
    table_T.to_string()
)

print()


print("=" * 70)
print("表2：水分浓度 / kg/kg")
print("=" * 70)

print(
    table_C.to_string()
)

print()


# 同时保存两个论文小表
with pd.ExcelWriter(
    "Q1_table1_table2.xlsx",
    engine="openpyxl"
) as writer:

    table_T.to_excel(
        writer,
        sheet_name="表1_温度"
    )

    table_C.to_excel(
        writer,
        sheet_name="表2_水分浓度"
    )


# ============================================================
# 17. 重点数值检查
# ============================================================

print("=" * 70)
print("关键位置检查")
print("=" * 70)


idx_surface = N - 1


print(
    "C(100 s, 2.0 cm) =",
    f"{C_all[100, idx_surface]:.10f}"
)

print(
    "C(300 s, 2.0 cm) =",
    f"{C_all[300, idx_surface]:.10f}"
)

print(
    "C(600 s, 2.0 cm) =",
    f"{C_all[600, idx_surface]:.10f}"
)

print(
    "C(1800 s, 2.0 cm) =",
    f"{C_all[1800, idx_surface]:.10f}"
)


print()


# ============================================================
# 18. 基础物理检查
# ============================================================

print("=" * 70)
print("物理检查")
print("=" * 70)


print(
    "初始温度范围：",
    T_all[0].min(),
    T_all[0].max()
)

print(
    "1800 s 温度范围：",
    T_all[-1].min(),
    T_all[-1].max()
)

print(
    "初始水分浓度范围：",
    C_all[0].min(),
    C_all[0].max()
)

print(
    "1800 s 水分浓度范围：",
    C_all[-1].min(),
    C_all[-1].max()
)

print()


# ============================================================
# 19. 绘制论文温度时空场
# ============================================================

plot_space_idx = np.linspace(
    0,
    N - 1,
    401
).astype(int)


r_plot_cm = (
    r[
        plot_space_idx
    ]
    *
    100.0
)


# ------------------------------------------------------------
# 温度场
# ------------------------------------------------------------

plt.figure(
    figsize=(
        8.5,
        5.5
    )
)


mesh_T = plt.pcolormesh(
    sol_T.t,
    r_plot_cm,
    T_all[
        :,
        plot_space_idx
    ].T,

    shading="auto"
)


plt.colorbar(
    mesh_T,
    label="Temperature / °C"
)

plt.xlabel(
    "Time / s"
)

plt.ylabel(
    "Distance from center / cm"
)

plt.tight_layout()

plt.savefig(
    "fig_q1_temperature_field_N2561.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "fig_q1_temperature_field_N2561.png",
    dpi=400,
    bbox_inches="tight"
)

plt.close()


# ------------------------------------------------------------
# 水分浓度场
# ------------------------------------------------------------

plt.figure(
    figsize=(
        8.5,
        5.5
    )
)


mesh_C = plt.pcolormesh(
    sol_C.t,
    r_plot_cm,
    C_all[
        :,
        plot_space_idx
    ].T,

    shading="auto"
)


plt.colorbar(
    mesh_C,
    label="Moisture concentration / kg/kg"
)

plt.xlabel(
    "Time / s"
)

plt.ylabel(
    "Distance from center / cm"
)

plt.tight_layout()

plt.savefig(
    "fig_q1_moisture_field_N2561.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "fig_q1_moisture_field_N2561.png",
    dpi=400,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 20. 五个位置的曲线
# ============================================================

plt.figure(
    figsize=(
        8,
        5
    )
)


for r_cm in [
    0,
    0.5,
    1.0,
    1.5,
    2.0
]:

    idx = int(
        round(
            (r_cm / 100.0)
            /
            dr
        )
    )

    plt.plot(
        sol_T.t,
        T_all[:, idx],
        label=f"r={r_cm} cm"
    )


plt.xlabel(
    "Time / s"
)

plt.ylabel(
    "Temperature / °C"
)

plt.legend()

plt.grid(
    alpha=0.25
)

plt.tight_layout()

plt.savefig(
    "Q1_temperature_N2561.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


plt.figure(
    figsize=(
        8,
        5
    )
)


for r_cm in [
    0,
    0.5,
    1.0,
    1.5,
    2.0
]:

    idx = int(
        round(
            (r_cm / 100.0)
            /
            dr
        )
    )

    plt.plot(
        sol_C.t,
        C_all[:, idx],
        label=f"r={r_cm} cm"
    )


plt.xlabel(
    "Time / s"
)

plt.ylabel(
    "Moisture concentration / kg/kg"
)

plt.legend()

plt.grid(
    alpha=0.25
)

plt.tight_layout()

plt.savefig(
    "Q1_moisture_N2561.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 21. 最终输出
# ============================================================

total_runtime = (
    time.time()
    -
    start_clock
)


print("=" * 70)
print("Q1 N=2561 全部完成")
print("=" * 70)

print(
    f"总运行时间 = "
    f"{total_runtime:.3f} s"
)

print()

print("生成文件：")
print("1. result1.xlsx")
print("2. Q1_table1_table2.xlsx")
print("3. fig_q1_temperature_field_N2561.pdf")
print("4. fig_q1_moisture_field_N2561.pdf")
print("5. Q1_temperature_N2561.png")
print("6. Q1_moisture_N2561.png")

print()