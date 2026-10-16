import time
import numpy as np
import pandas as pd

from scipy.linalg import solve_banded
from scipy.interpolate import PchipInterpolator


# ============================================================
# 0. 总设置
# ============================================================

N_SENS = 1281
DT_SENS = 7.5
PICARD_ITERS = 3

# 是否额外检查 Picard 3次 -> 5次
# 是否执行 Picard 迭代次数稳定性检查
RUN_PICARD_CHECK = True

# Picard 检查采用最终时间步长
DT_PICARD_CHECK = 7.5


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

    # 附件1结束后
    return 50.0


def C_air(t):
    if t <= t_env_max:
        return np.interp(t, time_env, C_env)

    # 附件1结束后
    return 0.05


# ============================================================
# 2. 读取附件2：问题四半径变化
# ============================================================

radius_df = pd.read_excel("附件2.xlsx")

time_R = radius_df["时间"].to_numpy(dtype=float)
R_cm_data = radius_df["半径"].to_numpy(dtype=float)

R_data = R_cm_data / 100.0

order = np.argsort(time_R)

time_R = time_R[order]
R_data = R_data[order]

t_R_max = np.max(time_R)

R_pchip = PchipInterpolator(
    time_R,
    R_data,
    extrapolate=False
)


def radius(t):
    if t <= t_R_max:
        return float(R_pchip(t))

    return float(R_data[-1])


# ============================================================
# 3. 公共参数
# ============================================================

H_BASE = 25.0
HM_BASE = 8e-7

T0 = 28.0
C0 = 2.55

THRESHOLD = 0.15

T_MAX = 10 * 24 * 3600.0


# ============================================================
# 4. 三对角求解器
# ============================================================

def solve_tridiagonal(lower, diag, upper, rhs):

    ab = np.zeros(
        (3, len(diag)),
        dtype=float
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


def harmonic_mean(a, b):

    return (
        2.0
        * a
        * b
        /
        (a + b + 1e-30)
    )


# ============================================================
# 5. 问题三材料物性
# ============================================================

def q3_rho(C):

    return (
        650.0
        +
        128.0 * C
    )


def q3_cp(C):

    return (
        1450.0
        +
        2736.0
        * C
        / (C + 1.0)
    )


def q3_k(C):

    return (
        0.21
        +
        0.38
        * C
        / (C + 1.0)
    )


def q3_D(C, T_celsius, D_scale=1.0):

    C_safe = np.maximum(
        C,
        1e-12
    )

    T_kelvin = (
        T_celsius
        +
        273.15
    )

    D0 = (
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

    return D_scale * D0


# ============================================================
# 6. 问题三求解器
# ============================================================

def solve_q3(
    h_scale=1.0,
    hm_scale=1.0,
    D_scale=1.0,
    N=1281,
    dt=15.0,
    picard_iters=3,
    verbose=False
):

    R = 0.02

    h = H_BASE * h_scale
    hm = HM_BASE * hm_scale

    dr = R / (N - 1)

    r = np.linspace(
        0.0,
        R,
        N
    )

    # --------------------------------------------------------
    # FVM 几何
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # 单步
    # --------------------------------------------------------

    def advance(
        T_old,
        C_old,
        t_new,
        dt_now
    ):

        T_guess = T_old.copy()
        C_guess = C_old.copy()

        for _ in range(picard_iters):

            Ta = T_air(t_new)
            Ca = C_air(t_new)

            # ==================================================
            # 温度方程
            # ==================================================

            rho = q3_rho(C_guess)
            cp = q3_cp(C_guess)

            k_node = q3_k(C_guess)

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

            # 表面 Robin
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

            # ==================================================
            # 水分方程
            # ==================================================

            D_node = q3_D(
                C_guess,
                T_new,
                D_scale=D_scale
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
            lower_C = np.zeros(N - 1)
            upper_C = np.zeros(N - 1)

            rhs_C = C_old.copy()

            # 中心
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

            # 内部
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

            # 表面 Robin
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

    # --------------------------------------------------------
    # 时间推进
    # --------------------------------------------------------

    T = np.full(
        N,
        T0,
        dtype=float
    )

    C = np.full(
        N,
        C0,
        dtype=float
    )

    t = 0.0

    while t < T_MAX:

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

        if not np.all(np.isfinite(T)):
            raise RuntimeError(
                "Q3 温度出现 NaN/Inf"
            )

        if not np.all(np.isfinite(C)):
            raise RuntimeError(
                "Q3 水分出现 NaN/Inf"
            )

        if np.min(C) < -1e-8:
            raise RuntimeError(
                "Q3 出现负水分浓度"
            )

        Cmax_new = np.max(C)

        if (
            Cmax_old > THRESHOLD
            and
            Cmax_new <= THRESHOLD
        ):

            fraction = (
                Cmax_old
                -
                THRESHOLD
            ) / (
                Cmax_old
                -
                Cmax_new
            )

            tcrit = (
                t_old
                +
                fraction * dt
            )

            if verbose:
                print(
                    "Q3 tcrit = "
                    f"{tcrit / 3600:.8f} h"
                )

            return tcrit

    raise RuntimeError(
        "Q3 10天内未达到阈值"
    )


# ============================================================
# 7. 问题四材料物性
# ============================================================

def q4_rho(C):

    return (
        760.0
        +
        90.0 * C
    )


def q4_cp(C):

    return (
        1850.0
        +
        2150.0
        * C
        / (C + 1.0)
    )


def q4_k(C):

    return (
        0.12
        +
        0.20
        * C
        / (C + 1.0)
    )


def q4_D(C, T_celsius, D_scale=1.0):

    C_safe = np.maximum(
        C,
        1e-12
    )

    T_kelvin = (
        T_celsius
        +
        273.15
    )

    D0 = (
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

    return D_scale * D0


# ============================================================
# 8. 问题四求解器
# ============================================================

def solve_q4(
    h_scale=1.0,
    hm_scale=1.0,
    D_scale=1.0,
    N=1281,
    dt=15.0,
    picard_iters=3,
    verbose=False
):

    h = H_BASE * h_scale
    hm = HM_BASE * hm_scale

    dxi = 1.0 / (N - 1)

    xi = np.linspace(
        0.0,
        1.0,
        N
    )

    # --------------------------------------------------------
    # xi方向 FVM 几何
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # 单步
    # --------------------------------------------------------

    def advance(
        T_old,
        C_old,
        t_new,
        dt_now
    ):

        T_guess = T_old.copy()
        C_guess = C_old.copy()

        for _ in range(picard_iters):

            Rt = radius(t_new)
            Rt2 = Rt**2

            Ta = T_air(t_new)
            Ca = C_air(t_new)

            # ==================================================
            # 温度方程
            # ==================================================

            rho = q4_rho(C_guess)
            cp = q4_cp(C_guess)

            k_node = q4_k(C_guess)

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

            # 表面 Robin
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

            # ==================================================
            # 水分方程
            # ==================================================

            D_node = q4_D(
                C_guess,
                T_new,
                D_scale=D_scale
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

            # 表面 Robin
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

    # --------------------------------------------------------
    # 时间推进
    # --------------------------------------------------------

    T = np.full(
        N,
        T0,
        dtype=float
    )

    C = np.full(
        N,
        C0,
        dtype=float
    )

    t = 0.0

    while t < T_MAX:

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

        if not np.all(np.isfinite(T)):
            raise RuntimeError(
                "Q4 温度出现 NaN/Inf"
            )

        if not np.all(np.isfinite(C)):
            raise RuntimeError(
                "Q4 水分出现 NaN/Inf"
            )

        if np.min(C) < -1e-8:
            raise RuntimeError(
                "Q4 出现负水分浓度"
            )

        Cmax_new = np.max(C)

        if (
            Cmax_old > THRESHOLD
            and
            Cmax_new <= THRESHOLD
        ):

            fraction = (
                Cmax_old
                -
                THRESHOLD
            ) / (
                Cmax_old
                -
                Cmax_new
            )

            tcrit = (
                t_old
                +
                fraction * dt
            )

            if verbose:
                print(
                    "Q4 tcrit = "
                    f"{tcrit / 3600:.8f} h"
                )

            return tcrit

    raise RuntimeError(
        "Q4 10天内未达到阈值"
    )


# ============================================================
# 9. 单个模型的敏感性分析
# ============================================================

def run_sensitivity(
    model_name,
    solver
):

    print()
    print("=" * 72)
    print(
        f"{model_name} 参数敏感性分析"
    )
    print("=" * 72)

    print(
        f"N = {N_SENS}, "
        f"dt = {DT_SENS:.1f} s, "
        f"Picard = {PICARD_ITERS}"
    )

    results = {}

    # --------------------------------------------------------
    # 基线
    # --------------------------------------------------------

    print("\n[1/7] baseline")

    t0 = solver(
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["baseline"] = t0

    print(
        f"tcrit = {t0 / 3600:.8f} h"
    )

    # --------------------------------------------------------
    # h -10%
    # --------------------------------------------------------

    print("\n[2/7] h -10%")

    t_h_minus = solver(
        h_scale=0.9,
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["h_-10%"] = t_h_minus

    print(
        f"tcrit = "
        f"{t_h_minus / 3600:.8f} h"
    )

    # --------------------------------------------------------
    # h +10%
    # --------------------------------------------------------

    print("\n[3/7] h +10%")

    t_h_plus = solver(
        h_scale=1.1,
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["h_+10%"] = t_h_plus

    print(
        f"tcrit = "
        f"{t_h_plus / 3600:.8f} h"
    )

    # --------------------------------------------------------
    # hm -10%
    # --------------------------------------------------------

    print("\n[4/7] hm -10%")

    t_hm_minus = solver(
        hm_scale=0.9,
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["hm_-10%"] = t_hm_minus

    print(
        f"tcrit = "
        f"{t_hm_minus / 3600:.8f} h"
    )

    # --------------------------------------------------------
    # hm +10%
    # --------------------------------------------------------

    print("\n[5/7] hm +10%")

    t_hm_plus = solver(
        hm_scale=1.1,
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["hm_+10%"] = t_hm_plus

    print(
        f"tcrit = "
        f"{t_hm_plus / 3600:.8f} h"
    )

    # --------------------------------------------------------
    # D -10%
    # --------------------------------------------------------

    print("\n[6/7] D -10%")

    t_D_minus = solver(
        D_scale=0.9,
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["D_-10%"] = t_D_minus

    print(
        f"tcrit = "
        f"{t_D_minus / 3600:.8f} h"
    )

    # --------------------------------------------------------
    # D +10%
    # --------------------------------------------------------

    print("\n[7/7] D +10%")

    t_D_plus = solver(
        D_scale=1.1,
        N=N_SENS,
        dt=DT_SENS,
        picard_iters=PICARD_ITERS
    )

    results["D_+10%"] = t_D_plus

    print(
        f"tcrit = "
        f"{t_D_plus / 3600:.8f} h"
    )

    # ========================================================
    # 归一化灵敏度
    # ========================================================

    S_h = (
        t_h_plus
        -
        t_h_minus
    ) / (
        0.2 * t0
    )

    S_hm = (
        t_hm_plus
        -
        t_hm_minus
    ) / (
        0.2 * t0
    )

    S_D = (
        t_D_plus
        -
        t_D_minus
    ) / (
        0.2 * t0
    )

    results["S_h"] = S_h
    results["S_hm"] = S_hm
    results["S_D"] = S_D

    # ========================================================
    # 相对变化
    # ========================================================

    results["h_-10%_change"] = (
        (t_h_minus - t0)
        /
        t0
        *
        100.0
    )

    results["h_+10%_change"] = (
        (t_h_plus - t0)
        /
        t0
        *
        100.0
    )

    results["hm_-10%_change"] = (
        (t_hm_minus - t0)
        /
        t0
        *
        100.0
    )

    results["hm_+10%_change"] = (
        (t_hm_plus - t0)
        /
        t0
        *
        100.0
    )

    results["D_-10%_change"] = (
        (t_D_minus - t0)
        /
        t0
        *
        100.0
    )

    results["D_+10%_change"] = (
        (t_D_plus - t0)
        /
        t0
        *
        100.0
    )

    print()
    print("-" * 72)
    print("归一化灵敏度")
    print("-" * 72)

    print(
        f"S_D  = {S_D:.8f}"
    )

    print(
        f"S_hm = {S_hm:.8f}"
    )

    print(
        f"S_h  = {S_h:.8e}"
    )

    print()
    print("绝对值排序：")

    ranking = sorted(
        [
            ("D", abs(S_D)),
            ("hm", abs(S_hm)),
            ("h", abs(S_h))
        ],
        key=lambda x: x[1],
        reverse=True
    )

    print(
        " > ".join(
            [
                x[0]
                for x in ranking
            ]
        )
    )

    return results


# ============================================================
# 10. 输出结果表
# ============================================================

def make_raw_table(
    model_name,
    result
):

    baseline = result["baseline"]

    rows = []

    cases = [
        ("baseline", 1.0, baseline),
        ("h -10%", 0.9, result["h_-10%"]),
        ("h +10%", 1.1, result["h_+10%"]),
        ("hm -10%", 0.9, result["hm_-10%"]),
        ("hm +10%", 1.1, result["hm_+10%"]),
        ("D -10%", 0.9, result["D_-10%"]),
        ("D +10%", 1.1, result["D_+10%"]),
    ]

    for case, scale, tcrit in cases:

        rows.append({
            "模型": model_name,
            "工况": case,
            "尺度系数": scale,
            "tcrit/s": tcrit,
            "tcrit/h": tcrit / 3600.0,
            "相对基线变化/%":
                (
                    (tcrit - baseline)
                    /
                    baseline
                    *
                    100.0
                )
        })

    return pd.DataFrame(rows)


# ============================================================
# 11. Picard 次数检查
# ============================================================

def run_picard_check():

    print()
    print("=" * 72)
    print("Picard 迭代次数检查")
    print("=" * 72)

    rows = []

    for name, solver in [
        ("问题三", solve_q3),
        ("问题四", solve_q4)
    ]:

        print()
        print(
            f"{name}：3 次 Picard"
        )

        t3 = solver(
            N=1281,
            dt=DT_PICARD_CHECK,
            picard_iters=3
        )

        print(
            f"tcrit = {t3 / 3600:.10f} h"
        )

        print(
            f"{name}：5 次 Picard"
        )

        t5 = solver(
            N=1281,
            dt=DT_PICARD_CHECK,
            picard_iters=5
        )

        print(
            f"tcrit = {t5 / 3600:.10f} h"
        )

        diff_s = abs(
            t5 - t3
        )

        relative_pct = (
            diff_s
            /
            t3
            *
            100.0
        )

        print(
            "3次 -> 5次差值 = "
            f"{diff_s:.10e} s"
        )

        rows.append({
            "模型": name,
            "3次Picard tcrit/h":
                t3 / 3600.0,
            "5次Picard tcrit/h":
                t5 / 3600.0,
            "绝对差/s":
                diff_s,
            "相对差/%":
                relative_pct
        })

    return pd.DataFrame(rows)


# ============================================================
# 12. 主程序
# ============================================================

if __name__ == "__main__":

    start = time.time()

    print("=" * 72)
    print("Q3/Q4 参数敏感性分析")
    print("=" * 72)

    print(
        "附件1时间范围："
        f"{time_env[0]:.0f} ~ "
        f"{time_env[-1]:.0f} s"
    )

    print(
        "附件2初始半径："
        f"{radius(0) * 100:.6f} cm"
    )

    print(
        "附件2最终半径："
        f"{radius(t_R_max) * 100:.6f} cm"
    )

    # ========================================================
    # Q3
    # ========================================================

    q3_result = run_sensitivity(
        "问题三",
        solve_q3
    )

    # ========================================================
    # Q4
    # ========================================================

    q4_result = run_sensitivity(
        "问题四",
        solve_q4
    )

    # ========================================================
    # 汇总表
    # ========================================================

    summary = pd.DataFrame(
        [
            {
                "模型": "问题三",
                "S_D": q3_result["S_D"],
                "S_hm": q3_result["S_hm"],
                "S_h": q3_result["S_h"]
            },
            {
                "模型": "问题四",
                "S_D": q4_result["S_D"],
                "S_hm": q4_result["S_hm"],
                "S_h": q4_result["S_h"]
            }
        ]
    )

    raw_q3 = make_raw_table(
        "问题三",
        q3_result
    )

    raw_q4 = make_raw_table(
        "问题四",
        q4_result
    )

    raw_all = pd.concat(
        [raw_q3, raw_q4],
        ignore_index=True
    )

    # ========================================================
    # Picard 检查
    # ========================================================

    picard_df = None

    if RUN_PICARD_CHECK:

        picard_df = run_picard_check()

    # ========================================================
    # 写 Excel
    # ========================================================

    output_file = (
        "sensitivity_results.xlsx"
    )

    with pd.ExcelWriter(
        output_file,
        engine="openpyxl"
    ) as writer:

        summary.to_excel(
            writer,
            sheet_name="灵敏度汇总",
            index=False
        )

        raw_all.to_excel(
            writer,
            sheet_name="原始临界时间",
            index=False
        )

        if picard_df is not None:

            picard_df.to_excel(
                writer,
                sheet_name="Picard检查",
                index=False
            )

    # ========================================================
    # 最终打印
    # ========================================================

    print()
    print("=" * 72)
    print("论文表：归一化灵敏度")
    print("=" * 72)

    print(
        summary.to_string(
            index=False
        )
    )

    print()
    print(
        f"结果已保存：{output_file}"
    )

    elapsed = (
        time.time()
        -
        start
    )

    print(
        "总运行时间："
        f"{elapsed:.1f} s "
        f"({elapsed / 60:.2f} min)"
    )