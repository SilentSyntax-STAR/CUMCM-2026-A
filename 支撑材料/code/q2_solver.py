import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


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
    附件1数据范围内采用线性插值
    """
    return np.interp(t, time_env, T_air_data)


def C_air(t):
    """
    烘房水分浓度：
    附件1数据范围内采用线性插值
    """
    return np.interp(t, time_env, C_air_data)


# ============================================================
# 2. 基本参数
# ============================================================

R = 0.02               # m

# 正式计算采用的空间步长
dr = 0.00025           # 0.025 cm

N = int(round(R / dr)) + 1

r = np.linspace(
    0.0,
    R,
    N
)

# 表面对流系数
h = 25.0               # W/(m^2 K)
hm = 8e-7              # m/s

# 初始状态
T0 = 28.0              # °C
C0 = 2.55              # kg/kg


# ============================================================
# 3. Q2材料物性
# ============================================================

def rho_func(C):

    return 650.0 + 128.0 * C


def cp_func(C):

    C_safe = np.maximum(C, 1e-10)

    return (
        1450.0
        +
        2736.0 * C_safe / (C_safe + 1.0)
    )


def k_func(C):

    C_safe = np.maximum(C, 1e-10)

    return (
        0.21
        +
        0.38 * C_safe / (C_safe + 1.0)
    )


def D_func(C, T_celsius):

    C_safe = np.maximum(
        C,
        1e-10
    )

    T_kelvin = (
        T_celsius
        + 273.15
    )

    # 防止摄氏度误传
    if np.min(T_kelvin) < 250:
        raise ValueError(
            "温度单位异常：D公式必须使用Kelvin"
        )

    return (
        2.4e-3
        *
        np.exp(
            -0.45 / C_safe
        )
        *
        np.exp(
            -3850.0 / T_kelvin
        )
    )


def harmonic_mean(a, b):

    return (
        2.0 * a * b
        /
        (a + b + 1e-30)
    )


# ============================================================
# 4. 有限体积控制体几何
# ============================================================

r_left = np.maximum(
    r - dr / 2,
    0.0
)

r_right = np.minimum(
    r + dr / 2,
    R
)

V = 0.5 * (
    r_right**2
    -
    r_left**2
)


# ============================================================
# 5. ODE右端
# ============================================================

def rhs(t, y):

    T = y[:N]
    C = y[N:]

    dTdt = np.zeros(N)
    dCdt = np.zeros(N)

    Ta = T_air(t)
    Ca = C_air(t)

    # --------------------------------------------------------
    # 当前状态下更新物性
    # --------------------------------------------------------

    rho_node = rho_func(C)
    cp_node = cp_func(C)
    k_node = k_func(C)

    D_node = D_func(
        C,
        T
    )

    # --------------------------------------------------------
    # 面值：采用调和平均
    # --------------------------------------------------------

    k_face = harmonic_mean(
        k_node[:-1],
        k_node[1:]
    )

    D_face = harmonic_mean(
        D_node[:-1],
        D_node[1:]
    )


    # ========================================================
    # 温度场
    # ========================================================

    for i in range(N):

        # 左界面
        if i == 0:

            flux_left_T = 0.0

        else:

            rf = (
                r[i - 1] + r[i]
            ) / 2

            flux_left_T = (
                rf
                *
                k_face[i - 1]
                *
                (T[i] - T[i - 1])
                / dr
            )

        # 右界面
        if i == N - 1:

            # 表面第三类边界
            flux_right_T = (
                R
                *
                h
                *
                (Ta - T[i])
            )

        else:

            rf = (
                r[i] + r[i + 1]
            ) / 2

            flux_right_T = (
                rf
                *
                k_face[i]
                *
                (T[i + 1] - T[i])
                / dr
            )

        dTdt[i] = (
            flux_right_T
            -
            flux_left_T
        ) / (
            rho_node[i]
            *
            cp_node[i]
            *
            V[i]
        )


    # ========================================================
    # 水分场
    # ========================================================

    for i in range(N):

        # 左界面
        if i == 0:

            flux_left_C = 0.0

        else:

            rf = (
                r[i - 1] + r[i]
            ) / 2

            flux_left_C = (
                rf
                *
                D_face[i - 1]
                *
                (C[i] - C[i - 1])
                / dr
            )

        # 右界面
        if i == N - 1:

            flux_right_C = (
                R
                *
                hm
                *
                (Ca - C[i])
            )

        else:

            rf = (
                r[i] + r[i + 1]
            ) / 2

            flux_right_C = (
                rf
                *
                D_face[i]
                *
                (C[i + 1] - C[i])
                / dr
            )

        dCdt[i] = (
            flux_right_C
            -
            flux_left_C
        ) / V[i]


    return np.concatenate(
        [dTdt, dCdt]
    )


# ============================================================
# 6. 初值
# ============================================================

T_init = np.full(
    N,
    T0
)

C_init = np.full(
    N,
    C0
)

y0 = np.concatenate(
    [T_init, C_init]
)


# ============================================================
# 7. 时间
# ============================================================

# 3 h = 10800 s

t_end = 3 * 3600

# 题目要求每1s输出
t_eval = np.arange(
    0,
    t_end + 1,
    1,
    dtype=float
)


# ============================================================
# 8. BDF求解
# ============================================================

print("开始求解 Q2 ...")

sol = solve_ivp(
    rhs,
    t_span=(0.0, t_end),
    y0=y0,
    method="BDF",
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-9,
    max_step=60.0
)

if not sol.success:

    raise RuntimeError(
        sol.message
    )

print("Q2求解完成")


# ============================================================
# 9. 拆分结果
# ============================================================

T_all = (
    sol.y[:N, :]
    .T
)

C_all = (
    sol.y[N:, :]
    .T
)


# ============================================================
# 10. 输出0.1 cm间距
# ============================================================

output_r_cm = np.arange(
    0.0,
    2.0001,
    0.1
)

indices = [
    int(
        round(
            (rcm / 100.0)
            / dr
        )
    )
    for rcm in output_r_cm
]

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
    np.arange(0, t_end + 1)
)

df_C_validation.insert(
    0,
    "时间",
    np.arange(0, t_end + 1)
)

with pd.ExcelWriter(
    "result2_validation_full.xlsx",
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

print("result2_validation_full.xlsx 已生成")

# ============================================================
# 11. 保存result2
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
        t_end + 1
    )
)

df_C.insert(
    0,
    "时间",
    np.arange(
        0,
        t_end + 1
    )
)

df_T.iloc[:, 1:] = (
    df_T.iloc[:, 1:]
    .round(4)
)

df_C.iloc[:, 1:] = (
    df_C.iloc[:, 1:]
    .round(4)
)

with pd.ExcelWriter(
    "result2.xlsx",
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
    "result2.xlsx 已生成"
)


# ============================================================
# 12. 打印论文表3、表4
# ============================================================

# 0.5、1、1.5、2、2.5、3 h
report_hours = np.array([
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0
])

report_times = (
    report_hours
    * 3600
).astype(int)

report_r_cm = [
    0,
    0.5,
    1.0,
    1.5,
    2.0
]

report_indices = [
    int(
        round(
            (x / 100.0)
            / dr
        )
    )
    for x in report_r_cm
]


table_T = pd.DataFrame(
    T_all[
        report_times
    ][:, report_indices],
    index=report_hours,
    columns=report_r_cm
).round(4)


table_C = pd.DataFrame(
    C_all[
        report_times
    ][:, report_indices],
    index=report_hours,
    columns=report_r_cm
).round(4)


print("\n==========================")
print("表3：3小时内药材温度")
print("==========================")

print(table_T)


print("\n==========================")
print("表4：3小时内水分浓度")
print("==========================")

print(table_C)


# ============================================================
# 13. 基础物理检查
# ============================================================

print("\n==========================")
print("物理检查")
print("==========================")

print(
    "3h最低温度：",
    T_all[-1].min()
)

print(
    "3h最高温度：",
    T_all[-1].max()
)

print(
    "3h最低含水率：",
    C_all[-1].min()
)

print(
    "3h最高含水率：",
    C_all[-1].max()
)

print(
    "是否存在负含水率：",
    np.any(C_all < 0)
)

print(
    "是否存在NaN：",
    np.any(np.isnan(T_all))
    or
    np.any(np.isnan(C_all))
)


# ============================================================
# 14. 画温度曲线
# ============================================================

plt.figure()

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
            / dr
        )
    )

    plt.plot(
        sol.t / 3600,
        T_all[:, idx],
        label=f"r={r_cm} cm"
    )

plt.xlabel(
    "Time / h"
)

plt.ylabel(
    "Temperature / °C"
)

plt.legend()

plt.grid()

plt.tight_layout()

plt.savefig(
    "Q2_temperature.png",
    dpi=300
)

plt.close()


# ============================================================
# 15. 画水分曲线
# ============================================================

plt.figure()

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
            / dr
        )
    )

    plt.plot(
        sol.t / 3600,
        C_all[:, idx],
        label=f"r={r_cm} cm"
    )

plt.xlabel(
    "Time / h"
)

plt.ylabel(
    "Moisture concentration / kg/kg"
)

plt.legend()

plt.grid()

plt.tight_layout()

plt.savefig(
    "Q2_moisture.png",
    dpi=300
)

plt.close()


# ============================================================
# 16. 热湿耦合参数变化曲线
# ============================================================

# 取中心和表面，展示扩散系数变化
D_center = np.array([
    D_func(
        np.array([C_all[i, 0]]),
        np.array([T_all[i, 0]])
    )[0]
    for i in range(len(sol.t))
])

D_surface = np.array([
    D_func(
        np.array([C_all[i, -1]]),
        np.array([T_all[i, -1]])
    )[0]
    for i in range(len(sol.t))
])


plt.figure()

plt.plot(
    sol.t / 3600,
    D_center,
    label="center"
)

plt.plot(
    sol.t / 3600,
    D_surface,
    label="surface"
)

plt.xlabel(
    "Time / h"
)

plt.ylabel(
    "D / m²/s"
)

plt.legend()

plt.grid()

plt.tight_layout()

plt.savefig(
    "Q2_diffusion_coefficient.png",
    dpi=300
)

plt.close()