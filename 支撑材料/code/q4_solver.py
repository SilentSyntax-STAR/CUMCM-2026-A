import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.linalg import solve_banded
from scipy.interpolate import PchipInterpolator
from openpyxl import load_workbook


# ============================================================
# 1. 读取附件1：烘房环境
# ============================================================

env = pd.read_excel("附件1.xlsx")

time_env = env["时间"].to_numpy(dtype=float)
T_env = env["温度"].to_numpy(dtype=float)
C_env = env["水分浓度"].to_numpy(dtype=float)

t_env_max = np.max(time_env)


def T_air(t):
    if t <= t_env_max:
        return np.interp(t, time_env, T_env)

    # 附件1结束后的边界假设
    return 50.0


def C_air(t):
    if t <= t_env_max:
        return np.interp(t, time_env, C_env)

    # 附件1结束后的边界假设
    return 0.05


# ============================================================
# 2. 读取附件2：半径变化
# ============================================================

radius_df = pd.read_excel("附件2.xlsx")

time_R = radius_df["时间"].to_numpy(dtype=float)
R_cm_data = radius_df["半径"].to_numpy(dtype=float)

# cm -> m
R_data = R_cm_data / 100.0

# 保证按时间排序
order = np.argsort(time_R)

time_R = time_R[order]
R_data = R_data[order]

t_R_max = np.max(time_R)


# PCHIP，避免普通高阶插值产生振荡
R_pchip = PchipInterpolator(
    time_R,
    R_data,
    extrapolate=False
)


def radius(t):
    """
    返回当前半径，单位m
    """
    if t <= t_R_max:
        return float(R_pchip(t))

    # 附件2结束后保持最终半径
    return float(R_data[-1])


print("=" * 60)
print("附件2检查")
print("=" * 60)

print(
    f"初始半径 = {radius(0) * 100:.6f} cm"
)

print(
    f"72h半径 = {radius(72 * 3600) * 100:.6f} cm"
)


# ============================================================
# 3. 基本参数
# ============================================================

h = 25.0
hm = 8e-7

T0 = 28.0
C0 = 2.55


# ============================================================
# 4. 最终收敛数值参数
# ============================================================

dxi = 0.00078125
dt = 7.5

N = (
    int(
        round(1.0 / dxi)
    )
    + 1
)

xi = np.linspace(
    0.0,
    1.0,
    N
)


print()
print("=" * 60)
print("Q4 最终计算参数")
print("=" * 60)

print(
    f"dxi = {dxi:.8f}"
)

print(
    f"dt = {dt:.1f} s"
)

print(
    f"N = {N}"
)


# ============================================================
# 5. 附录4材料物性
# ============================================================

def rho_func(C):
    return (
        760.0
        +
        90.0 * C
    )


def cp_func(C):
    return (
        1850.0
        +
        2150.0
        * C
        / (C + 1.0)
    )


def k_func(C):
    return (
        0.12
        +
        0.20
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
        4.2e-4
        *
        np.exp(
            -0.30 / C_safe
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
# 6. xi方向有限体积几何
# ============================================================

xi_left = np.maximum(
    xi - dxi / 2.0,
    0.0
)

xi_right = np.minimum(
    xi + dxi / 2.0,
    1.0
)

Vxi = 0.5 * (
    xi_right**2
    -
    xi_left**2
)

xif = 0.5 * (
    xi[:-1]
    +
    xi[1:]
)


# ============================================================
# 7. 三对角求解器
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
# 8. 单时间步推进
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

        Rt = radius(t_new)
        Rt2 = Rt**2

        Ta = T_air(t_new)
        Ca = C_air(t_new)


        # ====================================================
        # 温度
        # ====================================================

        rho = rho_func(C_guess)
        cp = cp_func(C_guess)

        k_node = k_func(C_guess)

        k_face = harmonic_mean(
            k_node[:-1],
            k_node[1:]
        )

        GT = (
            xif
            *
            k_face
            /
            dxi
        )

        mass_T = (
            rho
            *
            cp
            *
            Rt2
            *
            Vxi
        )


        diag_T = np.ones(N)
        lower_T = np.zeros(N - 1)
        upper_T = np.zeros(N - 1)

        rhs_T = T_old.copy()


        # 中心
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


        # 内部
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


        # 表面Robin边界
        HT = Rt * h

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
        # 水分
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
            xif
            *
            D_face
            /
            dxi
        )

        mass_C = (
            Rt2
            *
            Vxi
        )


        diag_C = np.ones(N)
        lower_C = np.zeros(N - 1)
        upper_C = np.zeros(N - 1)

        rhs_C = C_old.copy()


        # 中心
        diag_C[0] += (
            dt_now
            *
            GC[0]
            /
            mass_C[0]
        )

        upper_C[0] = (
            -dt_now
            *
            GC[0]
            /
            mass_C[0]
        )


        # 内部
        G_left = GC[:-1]
        G_right = GC[1:]

        M = mass_C[1:-1]

        lower_C[:-1] = (
            -dt_now
            *
            G_left
            /
            M
        )

        upper_C[1:] = (
            -dt_now
            *
            G_right
            /
            M
        )

        diag_C[1:-1] += (
            dt_now
            *
            (G_left + G_right)
            /
            M
        )


        # 表面Robin
        HC = Rt * hm

        lower_C[-1] = (
            -dt_now
            *
            GC[-1]
            /
            mass_C[-1]
        )

        diag_C[-1] += (
            dt_now
            *
            (GC[-1] + HC)
            /
            mass_C[-1]
        )

        rhs_C[-1] += (
            dt_now
            *
            HC
            *
            Ca
            /
            mass_C[-1]
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
# 9. 给定实际半径位置，提取水分浓度
# ============================================================

def get_C_at_physical_r(
    C_field,
    t,
    r_cm
):
    """
    r_cm：实际空间位置，单位cm

    如果该位置已经在收缩后的药材外：
    返回None

    否则由 xi = r/R(t) 插值得到C
    """

    Rt_cm = (
        radius(t)
        *
        100.0
    )

    # 超出当前药材半径
    if r_cm > Rt_cm + 1e-10:
        return None

    # 当前表面
    if abs(
        r_cm - Rt_cm
    ) < 1e-10:
        return float(
            C_field[-1]
        )

    xi_target = (
        r_cm
        /
        Rt_cm
    )

    return float(
        np.interp(
            xi_target,
            xi,
            C_field
        )
    )


# ============================================================
# 10. result4完整输出的空间位置
#
# 0, 0.1, ..., 2.0 cm
# 外部位置留空
# 最后另加“药材表面”
# ============================================================

output_r_cm = np.arange(
    0.0,
    2.0001,
    0.1
)


# ============================================================
# 11. 开始计算
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


# 每60s存一次
save_times = []
save_R_cm = []
save_C = []
save_surface_C = []


# 用于画图
plot_times = [0.0]
plot_Cmax = [C0]
plot_R = [radius(0) * 100]


# 与Q3一致，从60s开始保存
next_save_time = 60.0


tcrit = None
t_dry = None


print()
print("开始求解 Q4 最终模型...")


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

    Cmax_old = np.max(C_old)


    t_new = t + dt


    T, C = advance(
        T_old,
        C_old,
        t_new,
        dt
    )


    t = t_new

    Cmax_new = np.max(C)


    # ========================================================
    # 基本检查
    # ========================================================

    if np.any(~np.isfinite(C)):
        raise RuntimeError(
            "水分浓度出现NaN或Inf"
        )

    if np.any(~np.isfinite(T)):
        raise RuntimeError(
            "温度出现NaN或Inf"
        )

    if np.min(C) < -1e-8:
        raise RuntimeError(
            "出现负水分浓度"
        )


    # ========================================================
    # 每60s保存
    # ========================================================

    if abs(
        t - next_save_time
    ) < 1e-8:

        values = []

        for rr in output_r_cm:

            value = get_C_at_physical_r(
                C,
                t,
                rr
            )

            values.append(value)


        save_times.append(
            int(round(t))
        )

        save_R_cm.append(
            radius(t) * 100
        )

        save_C.append(
            values
        )

        # 当前真实表面
        save_surface_C.append(
            float(C[-1])
        )

        next_save_time += 60.0


    # ========================================================
    # 每5min保存画图数据
    # ========================================================

    if (
        int(round(t))
        %
        300
        ==
        0
    ):

        plot_times.append(t)
        plot_Cmax.append(Cmax_new)
        plot_R.append(
            radius(t) * 100
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
            Cmax_old - 0.15
        ) / (
            Cmax_old - Cmax_new
        )

        tcrit = (
            t_old
            +
            fraction * dt
        )


        # 向上取到下一整分钟
        t_dry = (
            np.ceil(
                tcrit / 60.0
            )
            *
            60.0
        )


        # 若恰好等于整分钟，
        # 为严格满足 C<0.15 再加1分钟
        if abs(
            t_dry - tcrit
        ) < 1e-8:

            t_dry += 60.0


        print()
        print("=" * 60)
        print("检测到理论临界时间")
        print("=" * 60)

        print(
            f"{tcrit:.3f} s"
        )

        print(
            f"{tcrit / 3600:.8f} h"
        )

        print(
            f"临界半径 = "
            f"{radius(tcrit) * 100:.6f} cm"
        )

        print()
        print("实际结束时间：")

        print(
            f"{t_dry:.0f} s"
        )

        print(
            f"{t_dry / 3600:.4f} h"
        )


    # ========================================================
    # 到实际结束时间后停止
    # ========================================================

    if (
        t_dry is not None
        and
        t >= t_dry - 1e-8
    ):

        break


if tcrit is None:
    raise RuntimeError(
        "10天内仍未达到烘干标准"
    )


# ============================================================
# 12. 最终检查
# ============================================================

imax = np.argmax(C)

Cmax_final = np.max(C)

R_final_cm = (
    radius(t_dry)
    *
    100
)

rmax_final_cm = (
    xi[imax]
    *
    R_final_cm
)


print()
print("=" * 60)
print("Q4 最终结果")
print("=" * 60)

print(
    f"理论临界时间 = "
    f"{tcrit / 3600:.8f} h"
)

print(
    f"理论临界时间 = "
    f"{tcrit:.3f} s"
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
    f"最终半径 = "
    f"{R_final_cm:.6f} cm"
)

print(
    f"最终Cmax = "
    f"{Cmax_final:.10f}"
)

print(
    f"最大值实际位置 = "
    f"{rmax_final_cm:.8f} cm"
)


if not (
    Cmax_final < 0.15
):
    raise RuntimeError(
        "最终时刻没有严格满足 Cmax < 0.15"
    )


# ============================================================
# 13. 转数组
# ============================================================

save_times = np.array(
    save_times,
    dtype=int
)

save_R_cm = np.array(
    save_R_cm,
    dtype=float
)

save_C = np.array(
    save_C,
    dtype=object
)

save_surface_C = np.array(
    save_surface_C,
    dtype=float
)


mask = (
    save_times
    <=
    int(t_dry)
)

save_times = save_times[mask]
save_R_cm = save_R_cm[mask]
save_C = save_C[mask]
save_surface_C = save_surface_C[mask]


# ============================================================
# 14. 写入官方result4.xlsx
# ============================================================

wb = load_workbook(
    "result4.xlsx"
)

print(
    "result4工作表：",
    wb.sheetnames
)


# 结果模板包含单个工作表，读取首个工作表
ws = wb[
    wb.sheetnames[0]
]

# ------------------------------------------------------------
# 表头
# ------------------------------------------------------------

ws.cell(
    row=1,
    column=1,
    value="时间\\到药材中心的距离"
)


for j, rr in enumerate(
    output_r_cm
):

    ws.cell(
        row=1,
        column=j + 2,
        value=float(
            round(rr, 1)
        )
    )


surface_col = (
    len(output_r_cm)
    +
    2
)

ws.cell(
    row=1,
    column=surface_col,
    value="药材表面"
)


# ------------------------------------------------------------
# 写数据
# ------------------------------------------------------------

for i, time_s in enumerate(
    save_times
):

    excel_row = i + 2


    ws.cell(
        row=excel_row,
        column=1,
        value=int(time_s)
    )


    # 固定实际位置0~2cm
    for j in range(
        len(output_r_cm)
    ):

        value = save_C[i, j]

        cell = ws.cell(
            row=excel_row,
            column=j + 2
        )


        # 收缩后位于药材外部：
        # Excel保持空白
        if value is None:

            cell.value = None

        else:

            cell.value = round(
                float(value),
                4
            )


    # 单独写真实表面
    ws.cell(
        row=excel_row,
        column=surface_col,
        value=round(
            float(
                save_surface_C[i]
            ),
            4
        )
    )


wb.save(
    "result4.xlsx"
)


print()
print(
    "result4.xlsx 已生成"
)


# ============================================================
# 15. 论文表6
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


report_r_cm = [
    0.0,
    0.5,
    1.0,
    1.5,
    2.0
]


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

    row = {}

    for rr in report_r_cm:

        col_index = int(
            round(
                rr / 0.1
            )
        )

        value = save_C[
            idx,
            col_index
        ]

        if value is None:

            row[
                f"{rr:.1f} cm"
            ] = np.nan

        else:

            row[
                f"{rr:.1f} cm"
            ] = float(value)


    row[
        "药材表面"
    ] = save_surface_C[idx]

    row[
        "当前半径/cm"
    ] = save_R_cm[idx]

    table_rows.append(
        row
    )


table6 = pd.DataFrame(
    table_rows,
    index=report_labels
)

table6.index.name = "时间/h"


print()
print("=" * 100)
print("Q4 每隔6h水分浓度")
print("=" * 100)

print(
    table6.round(4)
)


table6.to_excel(
    "Q4_table_6h.xlsx"
)


# ============================================================
# 16. 最终时刻关键位置
# ============================================================

print()
print("=" * 100)
print("Q4 烘干结束时关键位置水分浓度")
print("=" * 100)


for rr in report_r_cm:

    value = get_C_at_physical_r(
        C,
        t_dry,
        rr
    )


    if value is None:

        print(
            f"r={rr:.1f} cm："
            f"已位于药材外部"
        )

    else:

        print(
            f"r={rr:.1f} cm："
            f"C={value:.6f}"
        )


print(
    f"药材表面："
    f"C={C[-1]:.6f}"
)

print(
    f"最终半径："
    f"{R_final_cm:.6f} cm"
)


# ============================================================
# 17. Cmax阈值图
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
    label="max C"
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
    "Q4_Cmax_threshold_final.png",
    dpi=300
)

plt.close()


# ============================================================
# 18. 半径随时间变化图
# ============================================================

plot_R = np.array(
    plot_R
)


plt.figure(
    figsize=(8, 5)
)

plt.plot(
    plot_times / 3600,
    plot_R
)

plt.xlabel(
    "Time / h"
)

plt.ylabel(
    "Radius / cm"
)

plt.grid()

plt.tight_layout()

plt.savefig(
    "Q4_radius_final.png",
    dpi=300
)

plt.close()


print()
print(
    "Q4_Cmax_threshold_final.png 已生成"
)

print(
    "Q4_radius_final.png 已生成"
)

print(
    "Q4_table_6h.xlsx 已生成"
)