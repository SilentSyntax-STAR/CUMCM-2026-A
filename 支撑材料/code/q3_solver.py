import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.linalg import solve_banded
from openpyxl import load_workbook


# ============================================================
# 1. 读取附件1
# ============================================================

env = pd.read_excel("附件1.xlsx")

time_env = env["时间"].to_numpy(dtype=float)
T_env = env["温度"].to_numpy(dtype=float)
C_env = env["水分浓度"].to_numpy(dtype=float)

t_env_max = time_env[-1]


def T_air(t):
    """
    附件1范围内：线性插值
    超出附件1范围：50℃
    """
    if t <= t_env_max:
        return np.interp(t, time_env, T_env)

    return 50.0


def C_air(t):
    """
    附件1范围内：线性插值
    超出附件1范围：0.05 kg/kg
    """
    if t <= t_env_max:
        return np.interp(t, time_env, C_env)

    return 0.05


# ============================================================
# 2. 基本参数
# ============================================================

R = 0.02                   # m

h = 25.0                   # W/(m^2 K)
hm = 8e-7                  # m/s

T0 = 28.0                  # °C
C0 = 2.55                  # kg/kg


# ============================================================
# 3. 最终数值参数
#
# 空间网格：
# 0.0015625 cm
#
# 换成m：
# 0.0015625 / 100 = 0.000015625 m
# ============================================================

dr = 0.000015625           # m
dt = 7.5                   # s

N = int(round(R / dr)) + 1

r = np.linspace(
    0.0,
    R,
    N
)


print("================================")
print("Q3 最终计算参数")
print("================================")

print(
    f"dr = {dr * 100:.7f} cm"
)

print(
    f"dt = {dt:.1f} s"
)

print(
    f"N = {N}"
)


# ============================================================
# 4. 附录3材料物性
# ============================================================

def rho_func(C):

    return (
        650.0
        +
        128.0 * C
    )


def cp_func(C):

    return (
        1450.0
        +
        2736.0
        * C
        / (C + 1.0)
    )


def k_func(C):

    return (
        0.21
        +
        0.38
        * C
        / (C + 1.0)
    )


def D_func(C, T_celsius):

    C_safe = np.maximum(
        C,
        1e-12
    )

    T_kelvin = (
        T_celsius
        +
        273.15
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
        2.0
        * a
        * b
        /
        (a + b + 1e-30)
    )


# ============================================================
# 5. FVM几何
# ============================================================

r_left = np.maximum(
    r - dr / 2.0,
    0.0
)

r_right = np.minimum(
    r + dr / 2.0,
    R
)

V = 0.5 * (
    r_right**2
    -
    r_left**2
)

rf = 0.5 * (
    r[:-1]
    +
    r[1:]
)


# ============================================================
# 6. 三对角矩阵求解器
# ============================================================

def solve_tridiagonal(
    lower,
    diag,
    upper,
    rhs
):

    ab = np.zeros(
        (3, len(diag))
    )

    ab[0, 1:] = upper
    ab[1, :] = diag
    ab[2, :-1] = lower

    return solve_banded(
        (1, 1),
        ab,
        rhs,
        check_finite=False
    )


# ============================================================
# 7. 单个时间步
# ============================================================

def advance(
    T_old,
    C_old,
    t_new,
    dt_now
):

    T_guess = T_old.copy()
    C_guess = C_old.copy()

    # 三次Picard迭代
    for _ in range(3):

        Ta = T_air(t_new)
        Ca = C_air(t_new)


        # ====================================================
        # 温度方程
        # ====================================================

        rho = rho_func(C_guess)
        cp = cp_func(C_guess)

        k_node = k_func(
            C_guess
        )

        k_face = harmonic_mean(
            k_node[:-1],
            k_node[1:]
        )

        GT = (
            rf
            *
            k_face
            /
            dr
        )

        mass_T = (
            rho
            *
            cp
            *
            V
        )


        diag_T = np.ones(N)

        lower_T = np.zeros(
            N - 1
        )

        upper_T = np.zeros(
            N - 1
        )

        rhs_T = T_old.copy()


        # ----------------------------------------------------
        # 中心
        # ----------------------------------------------------

        diag_T[0] += (
            dt_now
            *
            GT[0]
            /
            mass_T[0]
        )

        upper_T[0] = (
            -dt_now
            *
            GT[0]
            /
            mass_T[0]
        )


        # ----------------------------------------------------
        # 内部节点
        # ----------------------------------------------------

        G_left = GT[:-1]
        G_right = GT[1:]

        M = mass_T[1:-1]

        lower_T[:-1] = (
            -dt_now
            *
            G_left
            /
            M
        )

        upper_T[1:] = (
            -dt_now
            *
            G_right
            /
            M
        )

        diag_T[1:-1] += (
            dt_now
            *
            (G_left + G_right)
            /
            M
        )


        # ----------------------------------------------------
        # 表面Robin边界
        # ----------------------------------------------------

        HT = R * h

        lower_T[-1] = (
            -dt_now
            *
            GT[-1]
            /
            mass_T[-1]
        )

        diag_T[-1] += (
            dt_now
            *
            (GT[-1] + HT)
            /
            mass_T[-1]
        )

        rhs_T[-1] += (
            dt_now
            *
            HT
            *
            Ta
            /
            mass_T[-1]
        )


        T_new = solve_tridiagonal(
            lower_T,
            diag_T,
            upper_T,
            rhs_T
        )


        # ====================================================
        # 水分方程
        # ====================================================

        D_node = D_func(
            C_guess,
            T_new
        )

        D_face = harmonic_mean(
            D_node[:-1],
            D_node[1:]
        )

        GC = (
            rf
            *
            D_face
            /
            dr
        )


        diag_C = np.ones(N)

        lower_C = np.zeros(
            N - 1
        )

        upper_C = np.zeros(
            N - 1
        )

        rhs_C = C_old.copy()


        # ----------------------------------------------------
        # 中心
        # ----------------------------------------------------

        diag_C[0] += (
            dt_now
            *
            GC[0]
            /
            V[0]
        )

        upper_C[0] = (
            -dt_now
            *
            GC[0]
            /
            V[0]
        )


        # ----------------------------------------------------
        # 内部
        # ----------------------------------------------------

        G_left = GC[:-1]
        G_right = GC[1:]

        Vi = V[1:-1]

        lower_C[:-1] = (
            -dt_now
            *
            G_left
            /
            Vi
        )

        upper_C[1:] = (
            -dt_now
            *
            G_right
            /
            Vi
        )

        diag_C[1:-1] += (
            dt_now
            *
            (G_left + G_right)
            /
            Vi
        )


        # ----------------------------------------------------
        # 表面Robin边界
        # ----------------------------------------------------

        HC = R * hm

        lower_C[-1] = (
            -dt_now
            *
            GC[-1]
            /
            V[-1]
        )

        diag_C[-1] += (
            dt_now
            *
            (GC[-1] + HC)
            /
            V[-1]
        )

        rhs_C[-1] += (
            dt_now
            *
            HC
            *
            Ca
            /
            V[-1]
        )


        C_new = solve_tridiagonal(
            lower_C,
            diag_C,
            upper_C,
            rhs_C
        )


        T_guess = T_new
        C_guess = C_new


    return T_new, C_new


# ============================================================
# 8. 输出位置
#
# 0,0.1,...,2 cm
#
# 0.1 / 0.0015625 = 64
# 因此全部正好落在计算节点上
# ============================================================

output_r_cm = np.arange(
    0.0,
    2.0001,
    0.1
)

output_indices = np.array([
    int(
        round(
            (x / 100.0)
            /
            dr
        )
    )
    for x in output_r_cm
])


# ============================================================
# 9. 开始时间推进
# ============================================================

T = np.full(
    N,
    T0
)

C = np.full(
    N,
    C0
)


t = 0.0


# 保存每60 s结果
save_times = []
save_C = []


# 同时保存用于画图的Cmax
plot_times = [0.0]
plot_Cmax = [
    np.max(C)
]


# 初始状态不写入result3，因为官方模板从60s开始
next_save_time = 60.0


tcrit = None
t_dry = None


print()
print("开始求解 Q3 最终模型...")


# 最大模拟时长设为10 d
t_max = (
    10
    *
    24
    *
    3600
)


while t < t_max:

    T_old = T.copy()
    C_old = C.copy()

    t_old = t

    Cmax_old = np.max(
        C_old
    )


    # 下一步
    t_new = (
        t
        +
        dt
    )


    T, C = advance(
        T_old,
        C_old,
        t_new,
        dt
    )


    t = t_new


    Cmax_new = np.max(
        C
    )


    # ========================================================
    # 物理检查
    # ========================================================

    if np.any(
        ~np.isfinite(C)
    ):

        raise RuntimeError(
            "出现 NaN 或 Inf"
        )


    if np.min(C) < -1e-8:

        raise RuntimeError(
            "出现负含水率"
        )


    # ========================================================
    # 每60s保存一次
    # ========================================================

    if abs(
        t - next_save_time
    ) < 1e-8:

        save_times.append(
            int(
                round(t)
            )
        )

        save_C.append(
            C[
                output_indices
            ].copy()
        )

        next_save_time += 60.0


    # 每5分钟记录一次Cmax画图
    if (
        int(round(t))
        %
        300
        ==
        0
    ):

        plot_times.append(
            t
        )

        plot_Cmax.append(
            Cmax_new
        )


    # ========================================================
    # 搜索理论临界时刻
    # ========================================================

    if (
        tcrit is None
        and
        Cmax_old > 0.15
        and
        Cmax_new <= 0.15
    ):

        fraction = (
            Cmax_old
            -
            0.15
        ) / (
            Cmax_old
            -
            Cmax_new
        )

        tcrit = (
            t_old
            +
            fraction
            *
            dt
        )


        # 严格低于0.15：
        # 向上取到下一分钟

        t_dry = (
            np.ceil(
                tcrit
                /
                60.0
            )
            *
            60.0
        )


        # 如果理论临界值恰好落在整分钟
        # 则再向后一分钟
        if abs(
            t_dry
            -
            tcrit
        ) < 1e-8:

            t_dry += 60.0


        print()
        print(
            "检测到理论临界时间："
        )

        print(
            f"{tcrit:.3f} s"
        )

        print(
            f"{tcrit / 3600:.8f} h"
        )

        print(
            "实际结束时间："
        )

        print(
            f"{t_dry:.0f} s"
        )

        print(
            f"{t_dry / 3600:.4f} h"
        )


    # ========================================================
    # 到实际结束时刻停止
    # ========================================================

    if (
        t_dry is not None
        and
        t >= t_dry
        -
        1e-8
    ):

        break


if tcrit is None:

    raise RuntimeError(
        "10天内仍未达到烘干要求"
    )


# ============================================================
# 10. 最终检查
# ============================================================

imax = np.argmax(
    C
)

Cmax_final = np.max(
    C
)


print()
print("================================")
print("Q3 最终结果")
print("================================")

print(
    f"理论临界时间 = "
    f"{tcrit / 3600:.8f} h"
)

print(
    f"实际结束时间 = "
    f"{t_dry / 3600:.4f} h"
)

print(
    f"实际结束时间 = "
    f"{int(t_dry)} s"
)

print(
    f"最终Cmax = "
    f"{Cmax_final:.10f}"
)

print(
    f"最大值位置 = "
    f"{r[imax] * 100:.7f} cm"
)


if not (
    Cmax_final < 0.15
):

    raise RuntimeError(
        "最终时刻没有严格满足 Cmax < 0.15"
    )


# ============================================================
# 11. 整理完整结果
# ============================================================

save_times = np.array(
    save_times,
    dtype=int
)

save_C = np.array(
    save_C
)


# 只保留到最终烘干时间
mask = (
    save_times
    <=
    int(t_dry)
)

save_times = (
    save_times[
        mask
    ]
)

save_C = (
    save_C[
        mask
    ]
)


# ============================================================
# 12. 直接填写官方 result3.xlsx
# ============================================================

wb = load_workbook("result3.xlsx")

print("result3工作表：", wb.sheetnames)

ws = wb[wb.sheetnames[0]]


# 官方表头
ws.cell(
    row=1,
    column=1,
    value="时间\\到药材中心的距离"
)


# 距离表头
for j, r_cm in enumerate(
    output_r_cm
):

    ws.cell(
        row=1,
        column=j + 2,
        value=float(
            round(
                r_cm,
                1
            )
        )
    )


# 数据
for i, time_s in enumerate(
    save_times
):

    row_excel = (
        i + 2
    )

    ws.cell(
        row=row_excel,
        column=1,
        value=int(
            time_s
        )
    )

    for j in range(
        len(output_r_cm)
    ):

        ws.cell(
            row=row_excel,
            column=j + 2,
            value=round(
                float(
                    save_C[i, j]
                ),
                4
            )
        )


wb.save(
    "result3.xlsx"
)


print()
print(
    "result3.xlsx 已生成"
)


# ============================================================
# 13. 论文表5
# ============================================================

regular_hours = list(
    np.arange(
        6.0,
        t_dry / 3600,
        6.0
    )
)

report_times_h = regular_hours + [
    t_dry / 3600
]

report_labels = [
    f"{h:g}"
    for h in regular_hours
] + [
    "烘干结束时间"
]


report_r_cm = np.array([
    0.0,
    0.5,
    1.0,
    1.5,
    2.0
])

report_cols = np.array([
    int(
        round(
            rcm / 0.1
        )
    )
    for rcm in report_r_cm
])


table_rows = []


for hour in report_times_h:

    ts = int(
        round(
            hour * 3600
        )
    )

    idx = (
        ts // 60
        -
        1
    )

    values = save_C[
        idx,
        report_cols
    ]

    table_rows.append(
        values
    )


table5 = pd.DataFrame(
    table_rows,
    index=report_labels,
    columns=report_r_cm
).round(4)

table5.index.name = "时间/h"


print()
print("================================")
print("表5：水分浓度")
print("================================")

print(
    table5
)


table5.to_excel(
    "Q3_table5_final.xlsx"
)


# ============================================================
# 14. 烘干结束时刻关键位置
# ============================================================

final_output_C = save_C[-1]

final_key_C = final_output_C[
    report_cols
]


print()
print("================================")
print("烘干结束时水分浓度")
print("================================")


for rr, cc in zip(
    report_r_cm,
    final_key_C
):

    print(
        f"r={rr:.1f} cm : "
        f"C={cc:.6f}"
    )


# ============================================================
# 15. 阈值图
# ============================================================

plot_times = np.array(
    plot_times
)

plot_Cmax = np.array(
    plot_Cmax
)


plt.figure(
    figsize=(8, 5)
)

plt.plot(
    plot_times / 3600,
    plot_Cmax,
    label="max C(r,t)"
)

plt.axhline(
    y=0.15,
    linestyle="--",
    label="threshold = 0.15"
)

plt.axvline(
    x=tcrit / 3600,
    linestyle="--",
    label=f"critical time = {tcrit / 3600:.3f} h"
)

plt.xlabel(
    "Time / h"
)

plt.ylabel(
    "Maximum moisture concentration / kg/kg"
)

plt.legend()

plt.grid()

plt.tight_layout()

plt.savefig(
    "Q3_Cmax_threshold_final.png",
    dpi=300
)

plt.close()