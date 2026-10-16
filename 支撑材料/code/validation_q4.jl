# ================================================================
# 问题4
# MWORKS.Syslab / Julia 独立交叉验证程序
# Python主程序：移动坐标FVM + 全隐式时间推进 + Picard + solve_banded
# 本程序：移动坐标二阶FDM + Method of Lines(MOL) + Syslab ode15s
# 两套程序使用同一连续移动边界模型、同一R(t)、物性和初边界条件，
# 但采用不同的空间离散与时间积分路线。
# Python结果只在Julia独立求解完成后用于误差比较，不参与Julia求解。
# 所有输出统一写入脚本目录下：Q4_Julia_Output/
# ================================================================

using TyBase
using TyMath
using TyPlot
using Statistics
using Printf

const SCRIPT_DIR = dirname(abspath(@__FILE__))
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "Q4_Julia_Output")
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
    alpha = (xq-x[i])/(x[i+1]-x[i])
    return y[i] + alpha*(y[i+1]-y[i])
end

function pchip_slopes(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n >= 2 || error("PCHIP至少需要两个数据点。")
    h = diff(x)
    delta = diff(y) ./ h
    d = zeros(Float64,n)
    if n == 2
        d[1] = delta[1]; d[2] = delta[1]
        return d
    end
    for i in 2:n-1
        dm, dp = delta[i-1], delta[i]
        if dm == 0.0 || dp == 0.0 || dm*dp <= 0.0
            d[i] = 0.0
        else
            w1 = 2.0*h[i] + h[i-1]
            w2 = h[i] + 2.0*h[i-1]
            d[i] = (w1+w2)/(w1/dm+w2/dp)
        end
    end
    d1 = ((2.0*h[1]+h[2])*delta[1]-h[1]*delta[2])/(h[1]+h[2])
    if d1*delta[1] <= 0.0
        d1 = 0.0
    elseif delta[1]*delta[2] < 0.0 && abs(d1) > 3.0*abs(delta[1])
        d1 = 3.0*delta[1]
    end
    d[1] = d1

    dn = ((2.0*h[end]+h[end-1])*delta[end]-h[end]*delta[end-1])/(h[end]+h[end-1])
    if dn*delta[end] <= 0.0
        dn = 0.0
    elseif delta[end]*delta[end-1] < 0.0 && abs(dn) > 3.0*abs(delta[end])
        dn = 3.0*delta[end]
    end
    d[end] = dn
    return d
end

function pchip_eval(x::Vector{Float64}, y::Vector{Float64}, d::Vector{Float64}, xq::Float64)
    if xq <= x[1]
        return y[1]
    elseif xq >= x[end]
        return y[end]
    end
    i = min(searchsortedlast(x,xq), length(x)-1)
    h = x[i+1]-x[i]
    s = (xq-x[i])/h
    h00 = 2.0*s^3 - 3.0*s^2 + 1.0
    h10 = s^3 - 2.0*s^2 + s
    h01 = -2.0*s^3 + 3.0*s^2
    h11 = s^3 - s^2
    return h00*y[i] + h10*h*d[i] + h01*y[i+1] + h11*h*d[i+1]
end

function main()
    env_file = first_existing(["附件1.xlsx", "附件1(1).xlsx", "attachment1.xlsx"])
    radius_file = first_existing(["附件2.xlsx", "附件2(1).xlsx", "attachment2.xlsx"])
    py_result_file = first_existing([
    "result4.xlsx"
]; required=false)

    println("环境数据：", env_file)
    println("半径数据：", radius_file)
    if py_result_file == ""
        println("Python结果：未找到，将只执行Julia独立求解。")
    else
        println("Python结果：", py_result_file)
    end
    println("输出目录：", OUTPUT_DIR)

    env = Float64.(xlsread(env_file,1,"A2:C242"))
    time_env = vec(env[:,1]); T_env = vec(env[:,2]); C_env = vec(env[:,3])
    t_env_max = time_env[end]

    T_air(t) = Float64(t) <= t_env_max ?
        interp_linear_scalar(Float64(t),time_env,T_env) : 50.0
    C_air(t) = Float64(t) <= t_env_max ?
        interp_linear_scalar(Float64(t),time_env,C_env) : 0.05

    rad = Float64.(xlsread(radius_file,1,"A2:B146"))
    time_R_raw = vec(rad[:,1]); R_cm_raw = vec(rad[:,2])
    p = sortperm(time_R_raw)
    time_R = time_R_raw[p]
    R_data = R_cm_raw[p] ./ 100.0
    t_R_max = time_R[end]
    R_slope = pchip_slopes(time_R,R_data)

    radius(t) = Float64(t) <= t_R_max ?
        pchip_eval(time_R,R_data,R_slope,Float64(t)) : R_data[end]

    # 半径插值自检
    rcheck_t = collect(0.0:60.0:t_R_max)
    rcheck = [radius(tt) for tt in rcheck_t]
    minimum(rcheck) > 0.0 || error("半径插值出现非正值。")
    any(diff(rcheck) .> 1.0e-12) && error("半径插值出现非物理反向增长。")
    println("半径PCHIP自检：PASS")

    # ---------------- 参数 ----------------
    H = 25.0
    HM = 8.0e-7
    T0 = 28.0
    C0 = 2.55
    THRESHOLD = 0.15

    # 与Python主网格N=1281不同的独立验证网格
    N = 321
    xi = collect(LinRange(0.0,1.0,N))
    dxi = 1.0/(N-1)

    OUTPUT_DT = 30.0
    T_MAX = 72.0*3600.0

    @printf("独立验证网格：N=%d, dxi=%.8f\n", N, dxi)
    @printf("ode15s输出间隔：%.1f s\n", OUTPUT_DT)

    rho(C) = 760.0 + 90.0*C
    cp(C) = 1850.0 + 2150.0*max(C,1.0e-12)/(max(C,1.0e-12)+1.0)
    kfun(C) = 0.12 + 0.20*max(C,1.0e-12)/(max(C,1.0e-12)+1.0)
    Dfun(C,Tc) = 4.2e-4 * exp(-0.30/max(C,1.0e-12)) * exp(-3850.0/(Tc+273.15))

    # ============================================================
    # 移动坐标xi=r/R(t)上的独立FDM/MOL离散
    # 连续模型与Python主程序保持一致，不另加入Rdot项。
    # 内部用二阶中心差分，表面Robin条件用虚节点消元。
    # ============================================================
    function rhs_q4(t,y)
        T = @view y[1:N]
        C = @view y[N+1:2N]
        dydt = zeros(Float64,2N)
        dT = @view dydt[1:N]
        dC = @view dydt[N+1:2N]

        Rt = radius(t)
        Rt2 = Rt^2
        Ta = T_air(t)
        Ca = C_air(t)

        # xi=0
        k0 = kfun(C[1])
        D0 = Dfun(C[1],T[1])
        dT[1] = 4.0*k0*(T[2]-T[1])/(rho(C[1])*cp(C[1])*Rt2*dxi^2)
        dC[1] = 4.0*D0*(C[2]-C[1])/(Rt2*dxi^2)

        # 内部节点
        @inbounds for i in 2:N-1
            Ci = C[i]
            ki = kfun(Ci)
            kim = kfun(C[i-1])
            kip = kfun(C[i+1])

            Di = Dfun(Ci,T[i])
            Dim = Dfun(C[i-1],T[i-1])
            Dip = Dfun(C[i+1],T[i+1])

            Txi = (T[i+1]-T[i-1])/(2.0*dxi)
            Txixi = (T[i+1]-2.0*T[i]+T[i-1])/dxi^2
            kxi = (kip-kim)/(2.0*dxi)
            heat_op = ki*Txixi + (kxi + ki/xi[i])*Txi
            dT[i] = heat_op/(rho(Ci)*cp(Ci)*Rt2)

            Cxi = (C[i+1]-C[i-1])/(2.0*dxi)
            Cxixi = (C[i+1]-2.0*C[i]+C[i-1])/dxi^2
            Dxi = (Dip-Dim)/(2.0*dxi)
            mass_op = Di*Cxixi + (Dxi + Di/xi[i])*Cxi
            dC[i] = mass_op/Rt2
        end

        # xi=1 表面Robin
        kN = kfun(C[N]); kNm1 = kfun(C[N-1]); kNm2 = kfun(C[N-2])
        DN = Dfun(C[N],T[N]); DNm1 = Dfun(C[N-1],T[N-1]); DNm2 = Dfun(C[N-2],T[N-2])

        gT = Rt*H/kN*(Ta-T[N])                 # T_xi(1)
        Txixi = 2.0*(T[N-1]-T[N]+dxi*gT)/dxi^2
        kxi = (3.0*kN - 4.0*kNm1 + kNm2)/(2.0*dxi)
        dT[N] = (kN*Txixi + (kxi+kN)*gT)/(rho(C[N])*cp(C[N])*Rt2)

        gC = Rt*HM/DN*(Ca-C[N])                # C_xi(1)
        Cxixi = 2.0*(C[N-1]-C[N]+dxi*gC)/dxi^2
        Dxi = (3.0*DN - 4.0*DNm1 + DNm2)/(2.0*dxi)
        dC[N] = (DN*Cxixi + (Dxi+DN)*gC)/Rt2

        return dydt
    end

    println("\n================ Q4 Syslab独立复算 ================")
    y0 = vcat(fill(T0,N),fill(C0,N))
    tspan = collect(0.0:OUTPUT_DT:T_MAX)
    opts = odeset(reltol=1e-7, abstol=1e-9, maxstep=60.0, stats=:on)
    @time tq4, yq4 = ode15s(rhs_q4,tspan,y0,opts)

    Tq4 = yq4[:,1:N]
    Cq4 = yq4[:,N+1:2N]
    Cmax = vec(maximum(Cq4,dims=2))
    @printf("T_MAX末端 Cmax = %.10f kg/kg\n", Cmax[end])

    icross = 0
    for i in 2:length(Cmax)
        if Cmax[i-1] > THRESHOLD && Cmax[i] <= THRESHOLD
            icross = i
            break
        end
    end
    icross == 0 && error("在T_MAX内未检测到 Cmax<=0.15，请检查模型或增大T_MAX。")

    t1,t2 = tq4[icross-1],tq4[icross]
    c1,c2 = Cmax[icross-1],Cmax[icross]
    frac = (c1-THRESHOLD)/(c1-c2)
    tcrit = t1 + frac*(t2-t1)

    tdry = ceil(tcrit/60.0)*60.0
    if abs(tdry-tcrit) < 1.0e-10
        tdry += 60.0
    end

    tdry_row = Int(round(tdry/OUTPUT_DT))+1
    Cmax_dry = maximum(@view Cq4[tdry_row,:])
    Rfinal_cm = radius(tdry)*100.0

    @printf("理论临界时间 = %.6f s = %.8f h\n", tcrit,tcrit/3600.0)
    @printf("实际结束时间 = %.0f s = %.4f h\n", tdry,tdry/3600.0)
    @printf("结束时半径 = %.6f cm\n", Rfinal_cm)
    @printf("结束时Cmax = %.10f kg/kg\n", Cmax_dry)

    # Python当前正式基准，仅用于求解后核验
    PY_TCRIT_H = 51.09391594
    PY_TDRY_H = 51.1000
    PY_RFINAL_CM = 1.200000
    PY_CMAX_FINAL = 0.1499887558

    err_tcrit_h = abs(tcrit/3600.0-PY_TCRIT_H)
    err_r = abs(Rfinal_cm-PY_RFINAL_CM)
    err_cmax = abs(Cmax_dry-PY_CMAX_FINAL)
    @printf("Python基准临界时间 = %.8f h\n", PY_TCRIT_H)
    @printf("临界时间差 = %.6f h = %.3f s\n", err_tcrit_h,err_tcrit_h*3600.0)
    @printf("最终半径差 = %.6e cm\n", err_r)
    @printf("最终Cmax差 = %.6e kg/kg\n", err_cmax)

    # 固定实际半径位置插值
    function C_at_physical_r(Crow,tsec::Float64,r_cm::Float64)
        Rt_cm = radius(tsec)*100.0
        if r_cm > Rt_cm + 1.0e-10
            return NaN
        elseif abs(r_cm-Rt_cm) < 1.0e-10
            return Crow[end]
        end
        xq = r_cm/Rt_cm
        # xi等距，线性插值
        if xq <= 0.0
            return Crow[1]
        elseif xq >= 1.0
            return Crow[end]
        end
        pos = xq/dxi + 1.0
        i = floor(Int,pos)
        i = max(1,min(i,N-1))
        alpha = pos-i
        return (1.0-alpha)*Crow[i] + alpha*Crow[i+1]
    end

    output_r_cm = collect(0.0:0.1:2.0)
    save_times = collect(60:60:Int(round(tdry)))

    # 完整60s输出CSV：time + 0~2cm + surface + radius
    full_csv = joinpath(OUTPUT_DIR,"Q4_moisture_60s_0p1cm_surface_radius.csv")
    open(full_csv,"w") do io
        print(io,"time_s")
        for rr in output_r_cm
            @printf(io,",%.1f",rr)
        end
        print(io,",surface,radius_cm\n")

        for tt in save_times
            row = Int(round(tt/OUTPUT_DT))+1
            Crow = @view Cq4[row,:]
            @printf(io,"%d",tt)
            for rr in output_r_cm
                val = C_at_physical_r(Crow,Float64(tt),rr)
                if isnan(val)
                    print(io,",")
                else
                    @printf(io,",%.12g",val)
                end
            end
            @printf(io,",%.12g,%.8f\n",Crow[end],radius(Float64(tt))*100.0)
        end
    end

    # 表6
    report_hours = Float64[]
    hh = 6.0
    while hh < tdry/3600.0 - 1.0e-12
        push!(report_hours,hh)
        hh += 6.0
    end
    push!(report_hours,tdry/3600.0)
    report_r_cm = [0.0,0.5,1.0,1.5,2.0]

    table6_csv = joinpath(OUTPUT_DIR,"Q4_table6_julia.csv")
    open(table6_csv,"w") do io
        println(io,"time_h,0.0cm,0.5cm,1.0cm,1.5cm,2.0cm,surface,radius_cm")
        for hh2 in report_hours
            tt = Int(round(hh2*3600.0))
            row = Int(round(tt/OUTPUT_DT))+1
            Crow = @view Cq4[row,:]
            @printf(io,"%.6f",hh2)
            for rr in report_r_cm
                val = C_at_physical_r(Crow,Float64(tt),rr)
                if isnan(val)
                    print(io,",")
                else
                    @printf(io,",%.10f",val)
                end
            end
            @printf(io,",%.10f,%.8f\n",Crow[end],radius(Float64(tt))*100.0)
        end
    end

    # 与Python result4正式四位小数输出比较：只比较双方都有数值的公共单元格
    compare_ok = false
    maxabs = NaN; mae = NaN; rmse = NaN; ncommon = 0
    if py_result_file != ""
        nrows = length(save_times)
        pyraw = xlsread(py_result_file,1,"A2:W$(nrows+1)")
        if size(pyraw,1) == nrows && size(pyraw,2) >= 23
            abs_sum = 0.0
            sq_sum = 0.0
            maxv = 0.0
            countv = 0
            for i in 1:nrows
                tt = save_times[i]
                rowj = Int(round(tt/OUTPUT_DT))+1
                Crow = @view Cq4[rowj,:]
                # 21个固定实际位置
                for j in 1:length(output_r_cm)
                    valj = C_at_physical_r(Crow,Float64(tt),output_r_cm[j])
                    cell = pyraw[i,j+1]
                    if !isnan(valj) && cell isa Number
                        vj4 = round(valj;digits=4)
                        e = vj4 - Float64(cell)
                        ae = abs(e)
                        abs_sum += ae
                        sq_sum += e^2
                        maxv = max(maxv,ae)
                        countv += 1
                    end
                end
                # 表面列W（第23列）
                cell_s = pyraw[i,23]
                if cell_s isa Number
                    e = round(Crow[end];digits=4)-Float64(cell_s)
                    ae = abs(e)
                    abs_sum += ae
                    sq_sum += e^2
                    maxv = max(maxv,ae)
                    countv += 1
                end
            end
            if countv > 0
                maxabs = maxv
                mae = abs_sum/countv
                rmse = sqrt(sq_sum/countv)
                ncommon = countv
                compare_ok = true
                @printf("Python四位小数公共单元格比较：n=%d, MaxAE=%.6g, MAE=%.6g, RMSE=%.6g kg/kg\n",
                        ncommon,maxabs,mae,rmse)
            end
        else
            println("找到Python结果，但尺寸与预期不一致，跳过全时空比较。")
        end
    end

    # 汇总文本
    summary_path = joinpath(OUTPUT_DIR,"Q4_validation_summary.txt")
    open(summary_path,"w") do io
        println(io,"Q4 MWORKS.Syslab / Julia 独立交叉验证")
        println(io,"数值路线：移动坐标二阶FDM + MOL + ode15s")
        @printf(io,"N = %d\n",N)
        @printf(io,"dxi = %.10f\n",dxi)
        @printf(io,"理论临界时间 = %.8f h\n",tcrit/3600.0)
        @printf(io,"实际结束时间 = %.4f h\n",tdry/3600.0)
        @printf(io,"结束半径 = %.8f cm\n",Rfinal_cm)
        @printf(io,"结束Cmax = %.10f kg/kg\n",Cmax_dry)
        @printf(io,"Python基准临界时间 = %.8f h\n",PY_TCRIT_H)
        @printf(io,"临界时间绝对差 = %.8f h = %.3f s\n",err_tcrit_h,err_tcrit_h*3600.0)
        @printf(io,"最终半径绝对差 = %.8e cm\n",err_r)
        @printf(io,"最终Cmax绝对差 = %.8e kg/kg\n",err_cmax)
        if compare_ok
            @printf(io,"公共比较单元格数 = %d\n",ncommon)
            @printf(io,"四位小数公共单元格MaxAE = %.6g kg/kg\n",maxabs)
            @printf(io,"四位小数公共单元格MAE = %.6g kg/kg\n",mae)
            @printf(io,"四位小数公共单元格RMSE = %.6g kg/kg\n",rmse)
        end
    end

    # 图1：Cmax阈值
    try
        f1 = figure()
        plot(tq4./3600.0,Cmax;linewidth=1.8)
        hold("on")
        yline(THRESHOLD,"--")
        xline(tcrit/3600.0,"--")
        hold("off")
        xlabel("时间 / h")
        ylabel("最大水分浓度 / (kg/kg)")
        grid("on")
        saveas(f1,joinpath(OUTPUT_DIR,"Q4_Cmax_threshold_julia.png"))
    catch err
        println("提示：Q4阈值图绘制失败，但数值结果已正常输出。")
        println(err)
    end

    # 图2：半径变化
    try
        rplot = [radius(tt)*100.0 for tt in tq4]
        f2 = figure()
        plot(tq4./3600.0,rplot;linewidth=1.8)
        xlabel("时间 / h")
        ylabel("药材半径 / cm")
        grid("on")
        saveas(f2,joinpath(OUTPUT_DIR,"Q4_radius_julia.png"))
    catch err
        println("提示：Q4半径图绘制失败，但数值结果已正常输出。")
        println(err)
    end

    println("\nQ4验证完成，所有文件均已写入：",OUTPUT_DIR)
    println("  ",full_csv)
    println("  ",table6_csv)
    println("  ",summary_path)
end

main()
