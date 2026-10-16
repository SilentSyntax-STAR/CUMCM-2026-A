# ================================================================
# 问题3
# MWORKS.Syslab / Julia 独立交叉验证程序
# Python主程序：有限体积法(FVM) + 全隐式时间推进 + Picard + solve_banded
# 本程序：二阶有限差分(FDM) + Method of Lines(MOL) + Syslab ode15s
# 两套程序使用同一连续热湿耦合模型、同一物性和同一初边界条件，
# 但采用不同的空间离散与时间积分路线。
# Python结果只在Julia独立求解完成后用于误差比较，不参与Julia求解。
# 所有输出统一写入脚本目录下：Q3_Julia_Output/
# ================================================================

using TyBase
using TyMath
using TyPlot
using Statistics
using Printf

const SCRIPT_DIR = dirname(abspath(@__FILE__))
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "Q3_Julia_Output")
mkpath(OUTPUT_DIR)

function first_existing(candidates::Vector{String}; required::Bool=true)
    roots = [SCRIPT_DIR, pwd()]
    for root in roots
        for name in candidates
            f = isabspath(name) ? name : joinpath(root, name)
            if isfile(f)
                return abspath(f)
            end
        end
    end
    if required
        error("未找到输入文件。候选文件：\n" * join(candidates, "\n"))
    end
    return ""
end

function interp_linear_scalar(xq::Float64, x::Vector{Float64}, y::Vector{Float64})
    if xq <= x[1]
        return y[1]
    elseif xq >= x[end]
        return y[end]
    end
    i = searchsortedlast(x, xq)
    i = min(i, length(x)-1)
    alpha = (xq - x[i]) / (x[i+1] - x[i])
    return y[i] + alpha * (y[i+1] - y[i])
end

function write_q3_full_csv(path::String, times::Vector{Int}, output_r_cm, output_idx,
                           Cmat, output_dt::Float64)
    open(path, "w") do io
        print(io, "time_s")
        for rr in output_r_cm
            @printf(io, ",%.1f", rr)
        end
        println(io)

        for tt in times
            row = Int(round(tt / output_dt)) + 1
            @printf(io, "%d", tt)
            for j in output_idx
                @printf(io, ",%.12g", Cmat[row, j])
            end
            println(io)
        end
    end
end

function write_q3_table5_csv(path::String, report_hours, report_r_cm, report_idx,
                             Cmat, output_dt::Float64)
    open(path, "w") do io
        print(io, "time_h")
        for rr in report_r_cm
            @printf(io, ",%.1fcm", rr)
        end
        println(io)

        for hh in report_hours
            tsec = Int(round(hh * 3600.0))
            row = Int(round(tsec / output_dt)) + 1
            @printf(io, "%.6f", hh)
            for j in report_idx
                @printf(io, ",%.10f", Cmat[row, j])
            end
            println(io)
        end
    end
end

function main()
    env_file = first_existing(["附件1.xlsx", "附件1(1).xlsx", "attachment1.xlsx"])
    py_result_file = first_existing([
    "result3.xlsx"
]; required=false)

    println("环境数据：", env_file)
    if py_result_file == ""
        println("Python结果：未找到，将只执行Julia独立求解。")
    else
        println("Python结果：", py_result_file)
    end
    println("输出目录：", OUTPUT_DIR)

    env = Float64.(xlsread(env_file, 1, "A2:C242"))
    time_env = vec(env[:, 1])
    T_env = vec(env[:, 2])
    C_env = vec(env[:, 3])
    t_env_max = time_env[end]

    T_air(t) = Float64(t) <= t_env_max ?
        interp_linear_scalar(Float64(t), time_env, T_env) : 50.0
    C_air(t) = Float64(t) <= t_env_max ?
        interp_linear_scalar(Float64(t), time_env, C_env) : 0.05

    # ---------------- 参数 ----------------
    R = 0.02
    H = 25.0
    HM = 8.0e-7
    T0 = 28.0
    C0 = 2.55
    THRESHOLD = 0.15

    # 独立复算采用N=321，与Python主程序N=1281不同
    N = 321
    r = collect(LinRange(0.0, R, N))
    dr = R / (N - 1)

    OUTPUT_DT = 30.0
    T_MAX = 72.0 * 3600.0

    @printf("独立验证网格：N=%d, dr=%.7f cm\n", N, dr*100.0)
    @printf("ode15s输出间隔：%.1f s\n", OUTPUT_DT)

    output_r_cm = collect(0.0:0.1:2.0)
    output_idx = [Int(round((x/100.0)/dr)) + 1 for x in output_r_cm]

    report_r_cm = [0.0, 0.5, 1.0, 1.5, 2.0]
    report_idx = [Int(round((x/100.0)/dr)) + 1 for x in report_r_cm]

    rho(C) = 650.0 + 128.0*C
    cp(C) = 1450.0 + 2736.0*max(C,1.0e-12)/(max(C,1.0e-12)+1.0)
    kfun(C) = 0.21 + 0.38*max(C,1.0e-12)/(max(C,1.0e-12)+1.0)
    Dfun(C,Tc) = 2.4e-3 * exp(-0.45/max(C,1.0e-12)) * exp(-3850.0/(Tc+273.15))

    # ============================================================
    # 独立复算采用FDM/MOL离散
    # PDE写成：
    # (1/r)d/dr(r k T_r) = k T_rr + (k_r + k/r) T_r
    # (1/r)d/dr(r D C_r) = D C_rr + (D_r + D/r) C_r
    # 内部使用二阶中心差分；表面Robin边界用虚节点消元。
    # ============================================================
    function rhs_q3(t, y)
        T = @view y[1:N]
        C = @view y[N+1:2N]
        dydt = zeros(Float64, 2N)
        dT = @view dydt[1:N]
        dC = @view dydt[N+1:2N]

        Ta = T_air(t)
        Ca = C_air(t)

        # r=0：圆柱算子极限；对称条件使一阶导项消失
        k0 = kfun(C[1])
        D0 = Dfun(C[1], T[1])
        dT[1] = 4.0*k0*(T[2]-T[1])/(rho(C[1])*cp(C[1])*dr^2)
        dC[1] = 4.0*D0*(C[2]-C[1])/dr^2

        # 内部节点：经典二阶中心有限差分
        @inbounds for i in 2:N-1
            Ci = C[i]
            ki = kfun(Ci)
            kim = kfun(C[i-1])
            kip = kfun(C[i+1])

            Di = Dfun(Ci, T[i])
            Dim = Dfun(C[i-1], T[i-1])
            Dip = Dfun(C[i+1], T[i+1])

            Tr = (T[i+1]-T[i-1])/(2.0*dr)
            Trr = (T[i+1]-2.0*T[i]+T[i-1])/dr^2
            kr = (kip-kim)/(2.0*dr)
            heat_op = ki*Trr + (kr + ki/r[i])*Tr
            dT[i] = heat_op/(rho(Ci)*cp(Ci))

            Cr = (C[i+1]-C[i-1])/(2.0*dr)
            Crr = (C[i+1]-2.0*C[i]+C[i-1])/dr^2
            Dr = (Dip-Dim)/(2.0*dr)
            mass_op = Di*Crr + (Dr + Di/r[i])*Cr
            dC[i] = mass_op
        end

        # r=R：Robin边界 + 二阶虚节点
        kN = kfun(C[N])
        kNm1 = kfun(C[N-1])
        kNm2 = kfun(C[N-2])
        DN = Dfun(C[N], T[N])
        DNm1 = Dfun(C[N-1], T[N-1])
        DNm2 = Dfun(C[N-2], T[N-2])

        gT = H/kN*(Ta-T[N])                    # T_r(R)
        Trr = 2.0*(T[N-1]-T[N]+dr*gT)/dr^2
        kr = (3.0*kN - 4.0*kNm1 + kNm2)/(2.0*dr)
        dT[N] = (kN*Trr + (kr + kN/R)*gT)/(rho(C[N])*cp(C[N]))

        gC = HM/DN*(Ca-C[N])                   # C_r(R)
        Crr = 2.0*(C[N-1]-C[N]+dr*gC)/dr^2
        Dr = (3.0*DN - 4.0*DNm1 + DNm2)/(2.0*dr)
        dC[N] = DN*Crr + (Dr + DN/R)*gC

        return dydt
    end

    println("\n================ Q3 Syslab独立复算 ================")
    y0 = vcat(fill(T0,N), fill(C0,N))
    tspan = collect(0.0:OUTPUT_DT:T_MAX)
    opts = odeset(reltol=1e-7, abstol=1e-9, maxstep=60.0, stats=:on)
    @time tq3, yq3 = ode15s(rhs_q3, tspan, y0, opts)

    Tq3 = yq3[:, 1:N]
    Cq3 = yq3[:, N+1:2N]
    Cmax = vec(maximum(Cq3, dims=2))

    @printf("T_MAX末端 Cmax = %.10f kg/kg\n", Cmax[end])

    icross = 0
    for i in 2:length(Cmax)
        if Cmax[i-1] > THRESHOLD && Cmax[i] <= THRESHOLD
            icross = i
            break
        end
    end
    icross == 0 && error("在T_MAX内未检测到 Cmax<=0.15，请检查模型或增大T_MAX。")

    t1, t2 = tq3[icross-1], tq3[icross]
    c1, c2 = Cmax[icross-1], Cmax[icross]
    frac = (c1-THRESHOLD)/(c1-c2)
    tcrit = t1 + frac*(t2-t1)

    tdry = ceil(tcrit/60.0)*60.0
    if abs(tdry-tcrit) < 1.0e-10
        tdry += 60.0
    end

    @printf("理论临界时间 = %.6f s = %.8f h\n", tcrit, tcrit/3600.0)
    @printf("实际结束时间 = %.0f s = %.4f h\n", tdry, tdry/3600.0)

    # Python当前正式基准，仅用于求解后核验
    PY_TCRIT_H = 57.49294733
    PY_TDRY_H = 57.5000
    err_tcrit_h = abs(tcrit/3600.0 - PY_TCRIT_H)
    @printf("Python基准临界时间 = %.8f h\n", PY_TCRIT_H)
    @printf("临界时间差 = %.6f h = %.3f s\n", err_tcrit_h, err_tcrit_h*3600.0)

    # 60s正式输出
    save_times = collect(60:60:Int(round(tdry)))
    full_csv = joinpath(OUTPUT_DIR, "Q3_moisture_60s_0p1cm.csv")
    write_q3_full_csv(full_csv, save_times, output_r_cm, output_idx, Cq3, OUTPUT_DT)

    # 表5
    report_hours = Float64[]
    hh = 6.0
    while hh < tdry/3600.0 - 1.0e-12
        push!(report_hours, hh)
        hh += 6.0
    end
    push!(report_hours, tdry/3600.0)
    table5_csv = joinpath(OUTPUT_DIR, "Q3_table5_julia.csv")
    write_q3_table5_csv(table5_csv, report_hours, report_r_cm, report_idx, Cq3, OUTPUT_DT)

    # 与Python正式四位小数结果比较
    compare_ok = false
    maxabs = NaN
    mae = NaN
    rmse = NaN
    timeerr = NaN
    if py_result_file != ""
        nrows = length(save_times)
        pyraw_any = xlsread(py_result_file, 1, "A2:V$(nrows+1)")
        if size(pyraw_any,1) == nrows && size(pyraw_any,2) == 22
            pyraw = Float64.(pyraw_any)
            pytime = pyraw[:,1]
            pyC = pyraw[:,2:end]
            jout = zeros(Float64, nrows, length(output_idx))
            for i in 1:nrows
                row = Int(round(save_times[i]/OUTPUT_DT)) + 1
                jout[i,:] .= Cq3[row, output_idx]
            end
            j4 = round.(jout; digits=4)
            diffv = j4 .- pyC
            av = abs.(diffv)
            maxabs = maximum(av)
            mae = mean(av)
            rmse = sqrt(mean(diffv.^2))
            timeerr = maximum(abs.(pytime .- Float64.(save_times)))
            compare_ok = true
            @printf("Python四位小数全时空比较：MaxAE=%.6g, MAE=%.6g, RMSE=%.6g kg/kg\n",
                    maxabs, mae, rmse)
        else
            println("找到Python结果，但尺寸与预期不一致，跳过全时空比较。")
        end
    end

    # 汇总文本
    summary_path = joinpath(OUTPUT_DIR, "Q3_validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "Q3 MWORKS.Syslab / Julia 独立交叉验证")
        println(io, "数值路线：二阶FDM + MOL + ode15s")
        @printf(io, "N = %d\n", N)
        @printf(io, "dr = %.8f cm\n", dr*100.0)
        @printf(io, "理论临界时间 = %.8f h\n", tcrit/3600.0)
        @printf(io, "实际结束时间 = %.4f h\n", tdry/3600.0)
        @printf(io, "Python基准临界时间 = %.8f h\n", PY_TCRIT_H)
        @printf(io, "临界时间绝对差 = %.8f h = %.3f s\n", err_tcrit_h, err_tcrit_h*3600.0)
        @printf(io, "T_MAX末端Cmax = %.10f kg/kg\n", Cmax[end])
        if compare_ok
            @printf(io, "Python正式结果时间列最大误差 = %.6e s\n", timeerr)
            @printf(io, "四位小数全时空MaxAE = %.6g kg/kg\n", maxabs)
            @printf(io, "四位小数全时空MAE = %.6g kg/kg\n", mae)
            @printf(io, "四位小数全时空RMSE = %.6g kg/kg\n", rmse)
        end
    end

    # 阈值图（绘图失败不影响数值结果输出）
    try
        fig = figure()
        plot(tq3./3600.0, Cmax; linewidth=1.8)
        hold("on")
        yline(THRESHOLD, "--")
        xline(tcrit/3600.0, "--")
        hold("off")
        xlabel("时间 / h")
        ylabel("最大水分浓度 / (kg/kg)")
        grid("on")
        saveas(fig, joinpath(OUTPUT_DIR, "Q3_Cmax_threshold_julia.png"))
    catch err
        println("提示：阈值图绘制失败，但数值结果已正常输出。")
        println(err)
    end

    println("\nQ3验证完成，所有文件均已写入：", OUTPUT_DIR)
    println("  ", full_csv)
    println("  ", table5_csv)
    println("  ", summary_path)
end

main()
