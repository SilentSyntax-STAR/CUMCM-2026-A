# ================================================================
# Q2 MWORKS.Syslab / Julia 独立交叉验证程序
#
# Python主程序：径向有限体积法(FVM) + SciPy BDF
# 本程序：二阶保守有限差分法(FDM/MOL) + Syslab ode15s
#
# Python主程序先生成：
#   result2_validation_full.xlsx
# 该文件保留未进行四位小数舍入的高精度结果，
# 仅用于与Syslab/Julia独立求解结果进行误差统计。
#
# Python结果不参与Syslab/Julia求解过程。
#
# 输出：
#   Q2_Syslab_validation_summary.xlsx
#   Q2_Syslab_validation_paper.txt
#   Q2_Syslab_vs_Python.png
# ================================================================

using TyBase
using TyMath
using TyPlot
using Statistics
using Printf

# ------------------------------------------------
# 0. 文件定位
# ------------------------------------------------

const SCRIPT_DIR = dirname(abspath(@__FILE__))
const ROOT_DIR = dirname(SCRIPT_DIR)

function first_existing(candidates::Vector{String})
    roots = [
        pwd(),
        SCRIPT_DIR,
        ROOT_DIR,
        joinpath(ROOT_DIR, "data")
    ]

    for root in roots
        for name in candidates
            f = isabspath(name) ? name : joinpath(root, name)
            if isfile(f)
                return abspath(f)
            end
        end
    end

    error(
        "未找到所需文件。请确认以下候选文件存在：\n" *
        join(candidates, "\n")
    )
end

ENV_FILE = first_existing([
    "附件1.xlsx",
    "附件1(1).xlsx",
    "attachment1.xlsx"
])

Q2_PY_FILE = first_existing([
    "result2_validation_full.xlsx"
])

println("环境数据文件：", ENV_FILE)
println("Q2 Python高精度验证结果：", Q2_PY_FILE)

# ------------------------------------------------
# 1. 读取附件1并做分段线性插值
# ------------------------------------------------
env = Float64.(xlsread(ENV_FILE, 1, "A2:C242"))
time_env = vec(env[:, 1])
T_air_data = vec(env[:, 2])
C_air_data = vec(env[:, 3])

function interp_linear_scalar(t::Float64, x::Vector{Float64}, v::Vector{Float64})
    if t <= x[1]
        return v[1]
    elseif t >= x[end]
        return v[end]
    end
    i = searchsortedlast(x, t)
    alpha = (t - x[i]) / (x[i + 1] - x[i])
    return v[i] + alpha * (v[i + 1] - v[i])
end

Tair(t) = interp_linear_scalar(Float64(t), time_env, T_air_data)
Cair(t) = interp_linear_scalar(Float64(t), time_env, C_air_data)

# ------------------------------------------------
# 2. 问题2参数与空间网格
# ------------------------------------------------
const R  = 0.02       # m
const dr = 0.00025    # m = 0.025 cm
const N  = Int(round(R / dr)) + 1
const r  = collect(LinRange(0.0, R, N))

const h   = 25.0      # W/(m^2*K)
const hm  = 8e-7      # m/s
const T0  = 28.0      # degC
const C0  = 2.55      # kg/kg

harmonic_mean(a, b) = 2.0 * a * b / (a + b + 1e-30)

rho(C) = 650.0 + 128.0 * C
cp(C)  = 1450.0 + 2736.0 * max(C,1e-10) / (max(C,1e-10) + 1.0)
k(C)   = 0.21 + 0.38 * max(C,1e-10) / (max(C,1e-10) + 1.0)

function D(C, T_celsius)
    Cs = max(C, 1e-10)
    Tk = T_celsius + 273.15
    if Tk < 250.0
        error("温度单位异常：D公式必须使用 Kelvin")
    end
    return 2.4e-3 * exp(-0.45 / Cs) * exp(-3850.0 / Tk)
end

const out_idx = collect(1:4:N)
const report_r_cols = [1, 6, 11, 16, 21]  # 0,0.5,1.0,1.5,2.0 cm

# ------------------------------------------------
# 3. 半离散ODE：热湿耦合FDM + Robin边界
# ------------------------------------------------
function rhs_q2(t, y)
    T = @view y[1:N]
    C = @view y[N+1:2N]

    dydt = zeros(Float64, 2N)
    dT = @view dydt[1:N]
    dC = @view dydt[N+1:2N]

    Ta = Tair(t)
    Ca = Cair(t)

    rhonode = rho.(C)
    cpnode  = cp.(C)
    knode   = k.(C)
    Dnode   = [D(C[i], T[i]) for i in 1:N]

    kface = [harmonic_mean(knode[i], knode[i+1]) for i in 1:N-1]
    Dface = [harmonic_mean(Dnode[i], Dnode[i+1]) for i in 1:N-1]

    # r=0：圆柱坐标中心极限
    dT[1] = 4.0 * kface[1] * (T[2] - T[1]) / (rhonode[1] * cpnode[1] * dr^2)
    dC[1] = 4.0 * Dface[1] * (C[2] - C[1]) / dr^2

    # 内部节点
    for i in 2:N-1
        rp = 0.5 * (r[i] + r[i+1])
        rm = 0.5 * (r[i-1] + r[i])

        heat_div = (
            rp * kface[i] * (T[i+1] - T[i]) / dr -
            rm * kface[i-1] * (T[i] - T[i-1]) / dr
        ) / (r[i] * dr)
        dT[i] = heat_div / (rhonode[i] * cpnode[i])

        mass_div = (
            rp * Dface[i] * (C[i+1] - C[i]) / dr -
            rm * Dface[i-1] * (C[i] - C[i-1]) / dr
        ) / (r[i] * dr)
        dC[i] = mass_div
    end

    # r=R：第三类边界，利用虚节点消元
    gT = h / knode[N] * (Ta - T[N])
    Tghost = T[N-1] + 2.0 * dr * gT
    dT[N] = knode[N] / (rhonode[N] * cpnode[N]) * (
        (Tghost - 2.0*T[N] + T[N-1]) / dr^2 + gT / R
    )

    gC = hm / Dnode[N] * (Ca - C[N])
    Cghost = C[N-1] + 2.0 * dr * gC
    dC[N] = Dnode[N] * (
        (Cghost - 2.0*C[N] + C[N-1]) / dr^2 + gC / R
    )

    return dydt
end

# ------------------------------------------------
# 4. Syslab ode15s独立复算
# ------------------------------------------------
println("\n================ Q2 Syslab独立复算 ================")
y0 = vcat(fill(T0, N), fill(C0, N))
tspan = collect(0.0:1.0:10800.0)
opts = odeset(reltol=1e-7, abstol=1e-9, maxstep=60.0, stats=:on)
@time tq2, yq2 = ode15s(rhs_q2, tspan, y0, opts)

Tq2 = yq2[:, 1:N]
Cq2 = yq2[:, N+1:2N]
Tq2_out = Tq2[:, out_idx]
Cq2_out = Cq2[:, out_idx]

# ------------------------------------------------
# 5. 读取Python高精度结果，仅用于交叉验证
# ------------------------------------------------
pyT_raw = Float64.(xlsread(Q2_PY_FILE, 1, "A2:V10802"))
pyC_raw = Float64.(xlsread(Q2_PY_FILE, 2, "A2:V10802"))
pyT = pyT_raw[:, 2:end]
pyC = pyC_raw[:, 2:end]

function check_shape_and_time(name, tspan, pyraw, sim)
    if size(pyraw, 1) != length(tspan)
        error("$name 行数不一致：Python=$(size(pyraw,1)), Syslab=$(length(tspan))")
    end
    if size(pyraw, 2) != 22
        error("$name 列数异常，应为时间列+21个半径位置，共22列")
    end
    if maximum(abs.(pyraw[:,1] .- tspan)) > 1e-8
        error("$name 时间列与Syslab输出时间不一致")
    end
    if size(sim, 2) != 21
        error("$name Syslab空间输出列数异常")
    end
end

check_shape_and_time("Q2温度", tspan, pyT_raw, Tq2_out)
check_shape_and_time("Q2水分", tspan, pyC_raw, Cq2_out)

# ------------------------------------------------
# 6. 误差统计
# ------------------------------------------------
function calc_metrics(sim, ref)
    ae = abs.(sim .- ref)
    re = ae ./ max.(abs.(ref), 1e-12)
    return (
        maxabs = maximum(ae),
        mae = mean(ae),
        rmse = sqrt(mean(ae.^2)),
        maxrelpct = 100.0 * maximum(re)
    )
end

# t=0.5,1.0,1.5,2.0,2.5,3.0 h
report_rows = [1801, 3601, 5401, 7201, 9001, 10801]

mT_all = calc_metrics(Tq2_out, pyT)
mC_all = calc_metrics(Cq2_out, pyC)
mT_rep = calc_metrics(Tq2_out[report_rows, report_r_cols], pyT[report_rows, report_r_cols])
mC_rep = calc_metrics(Cq2_out[report_rows, report_r_cols], pyC[report_rows, report_r_cols])

function print_metric(title, m)
    @printf("%-24s  MaxAbs=%12.6e  MAE=%12.6e  RMSE=%12.6e  MaxRel=%9.5f%%\n",
            title, m.maxabs, m.mae, m.rmse, m.maxrelpct)
end

println("\n================ Q2 全输出网格误差 ================")
print_metric("Q2 温度", mT_all)
print_metric("Q2 水分浓度", mC_all)

println("\n================ Q2 论文规定表格点误差 ================")
print_metric("Q2 温度", mT_rep)
print_metric("Q2 水分浓度", mC_rep)

# ------------------------------------------------
# 7. 输出误差汇总Excel与论文文本
# ------------------------------------------------
summary = Any[
    "Problem" "Variable" "Scope" "MaxAbs" "MAE" "RMSE" "MaxRelPercent";
    "Q2" "Temperature" "AllGrid"   mT_all.maxabs mT_all.mae mT_all.rmse mT_all.maxrelpct;
    "Q2" "Moisture"    "AllGrid"   mC_all.maxabs mC_all.mae mC_all.rmse mC_all.maxrelpct;
    "Q2" "Temperature" "ReportPts" mT_rep.maxabs mT_rep.mae mT_rep.rmse mT_rep.maxrelpct;
    "Q2" "Moisture"    "ReportPts" mC_rep.maxabs mC_rep.mae mC_rep.rmse mC_rep.maxrelpct
]
SUMMARY_XLSX = joinpath(
    ROOT_DIR,
    "Q2_Syslab_validation_summary.xlsx"
)

SUMMARY_TXT = joinpath(
    ROOT_DIR,
    "Q2_Syslab_validation_paper.txt"
)

VALIDATION_PNG = joinpath(
    ROOT_DIR,
    "Q2_Syslab_vs_Python.png"
)

xlswrite(SUMMARY_XLSX, summary)
println("已生成：", SUMMARY_XLSX)

open(SUMMARY_TXT, "w") do io
    println(io, "MWORKS.Syslab（Julia）问题2独立交叉验证结果（论文规定表格点）：")
    @printf(io, "温度：最大绝对误差 %.6g °C，RMSE %.6g °C，最大相对误差 %.6g%%。\n",
            mT_rep.maxabs, mT_rep.rmse, mT_rep.maxrelpct)
    @printf(io, "水分浓度：最大绝对误差 %.6g kg/kg，RMSE %.6g kg/kg，最大相对误差 %.6g%%。\n",
            mC_rep.maxabs, mC_rep.rmse, mC_rep.maxrelpct)
end
println("已生成：", SUMMARY_TXT)

# ------------------------------------------------
# 8. Python vs Syslab叠加验证图
#    r=0,1,2 cm；Python实线，Syslab虚线
# ------------------------------------------------
sel = collect(1:60:length(tspan))
cols = [1, 11, 21]
x = tspan[sel] ./ 3600.0

fig = figure()
subplot(2,1,1)
plot(x, pyT[sel, cols[1]], "-b"; linewidth=1.8)
hold("on")
plot(x, Tq2_out[sel, cols[1]], "--b"; linewidth=1.3)
plot(x, pyT[sel, cols[2]], "-g"; linewidth=1.8)
plot(x, Tq2_out[sel, cols[2]], "--g"; linewidth=1.3)
plot(x, pyT[sel, cols[3]], "-r"; linewidth=1.8)
plot(x, Tq2_out[sel, cols[3]], "--r"; linewidth=1.3)
hold("off")
grid("on")
xlabel("Time / h")
ylabel("Temperature / °C")
title("Q2 cross-validation — Temperature")
legend(["Python r=0", "Syslab r=0", "Python r=1", "Syslab r=1", "Python r=2", "Syslab r=2"];
       loc="best", ncol=2)

subplot(2,1,2)
plot(x, pyC[sel, cols[1]], "-b"; linewidth=1.8)
hold("on")
plot(x, Cq2_out[sel, cols[1]], "--b"; linewidth=1.3)
plot(x, pyC[sel, cols[2]], "-g"; linewidth=1.8)
plot(x, Cq2_out[sel, cols[2]], "--g"; linewidth=1.3)
plot(x, pyC[sel, cols[3]], "-r"; linewidth=1.8)
plot(x, Cq2_out[sel, cols[3]], "--r"; linewidth=1.3)
hold("off")
grid("on")
xlabel("Time / h")
ylabel("Moisture concentration / kg/kg")
title("Q2 cross-validation — Moisture")
legend(["Python r=0", "Syslab r=0", "Python r=1", "Syslab r=1", "Python r=2", "Syslab r=2"];
       loc="best", ncol=2)

tightlayout()
saveas(fig, VALIDATION_PNG)
println("已生成：", VALIDATION_PNG)

println("\nQ2验证完成。")
