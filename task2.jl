using Statistics
using Plots
using LaTeXStrings
using PlotThemes
using LinearAlgebra
using SparseArrays

gr();
theme(:default;)

include("integrator.jl")

L = 1
C = 900
K = 210
ρ = 2700
Δt = 0.1
Nₜ = 10000

Δx = 0.01
# Nₓ = floor(Int,L/Δx+1)
xs = 0:Δx:L

Δts = 0.001:0.0005:0.7
ts =  0:Nₜ
t=100


λ = K/(C*ρ)
a = λ*Δt/Δx^2
T0 = sin.((π*xs) / L)
Nₓ = length(T0)

save_folder = "saves_t2"

function fullFCTS()
    Teval = FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ);
    hmap_ftcs = heatmap(xs, ts.*Δt, Teval, ylabel="Time",xlabel="Rod domain", title="Rod evolution with FTCS")
    savefig(save_folder*"/rod_FTCS.pdf")
    display(hmap_ftcs)
    return Teval
end;

function fullEB()

    A =Tridiagonal([-a for _ in 1:Nₓ-1], [(1+2*a) for _ in 1:Nₓ], [-a for _ in 1:Nₓ-1])
    A[1,2] = 0;
    A[1,1] = 1;
    A[end,end-1] = 0;
    A[end,end] = 1;

    Teval = euler_backward_all(T0, A ,Nₜ, Nₓ);
    hmap_ftcs = heatmap(xs, ts.*Δt, Teval, ylabel="Time",xlabel="Rod domain", title="Rod evolution with Euler Backward")
    savefig(save_folder*"/rod_euler_back.pdf")
    display(hmap_ftcs)
    return Teval
end

function fullCN()
    a = λ*Δt/(2*Δx^2)
    A =Tridiagonal([-a for _ in 1:Nₓ-1], [(1+2*a) for _ in 1:Nₓ], [-a for _ in 1:Nₓ-1])
    A[1,2] = 0;
    A[1,1] = 1;
    A[end,end] = 1
    A[end,end-1] = 0
    B =Tridiagonal([a for _ in 1:Nₓ-1], [(1-2*a) for _ in 1:Nₓ], [a for _ in 1:Nₓ-1])
    B[1,2] = 0;
    B[1,1] = 1;
    B[end,end] = 1
    B[end,end-1] = 0
    Teval = crank_nicolson_all(T0, A, B, Nₜ, Nₓ);
    hmap_ftcs = heatmap(xs, ts.*Δt, Teval, ylabel="Time",xlabel="Rod domain", title="Rod evolution with Crank-Nicolson")
    savefig(save_folder*"/rod_crank_nic.pdf")
    display(hmap_ftcs)
    return Teval
end

function fullDF()
    a = 2*λ*Δt/(Δx^2)
    Teval = dufort_frankel(T0, a, Nₜ, Nₓ);
    hmap_ftcs = heatmap(xs, ts.*Δt, Teval, ylabel="Time",xlabel="Rod domain", title="Rod evolution with Dufort-Frankel")
    savefig(save_folder*"/rod_dufort_frankel.pdf")
    display(hmap_ftcs)
    return Teval
end

function errorDevFCTS(disp=false)

    ϵ_Ts = Float64[]
    res = 0
    Threads.@threads for Δt in Δts
        Nₜ = floor(Int,t/Δt+2)
        T_all = FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ);
        T = T_all[end,:]
        res = betterϵ(T,Nₜ*Δt, xs, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    error_plot = plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for FCTS scheme", yaxis=:log)
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_fcts.pdf"))

    if disp
        display(error_plot)
    end

    return ϵ_Ts
end;

function errorDevEB(disp=false)
    ϵ_Ts = Float64[]
    res = 0

    Threads.@threads for Δt in Δts
        Nₜ = floor(Int,t/Δt+2)
        a = λ*Δt/Δx^2
        A =Tridiagonal([-a for _ in 1:Nₓ-1], [(1+2*a) for _ in 1:Nₓ], [-a for _ in 1:Nₓ-1])
        A[1,2] = 0;
        A[1,1] = 1;
        A[end,end-1] = 0;
        A[end,end] = 1;

        T = euler_backward(T0,A, Nₜ);
        res = betterϵ(T,Nₜ*Δt, xs, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    error_plot = plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for Euler Backward scheme")
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_eb.pdf"))

    if disp
        display(error_plot)
    end

    return  ϵ_Ts
end;

function errorDevCN(disp=false)
    ϵ_Ts = Float64[]
    res = 0

    Threads.@threads for Δt in Δts
        # Nₜ = floor(Int,t/Δt)
        # aray creation
        # ts = 0:Δt:t
        Nₜ = floor(Int,t/Δt+2)
        # Nₜ = length(ts)

        a = λ*Δt/(2*Δx^2)
        A =Tridiagonal([-a for _ in 1:Nₓ-1], [(1+2*a) for _ in 1:Nₓ], [-a for _ in 1:Nₓ-1])
        A[1,2] = 0;
        A[1,1] = 1;
        A[end,end] = 1
        A[end,end-1] = 0
        B =Tridiagonal([a for _ in 1:Nₓ-1], [(1-2*a) for _ in 1:Nₓ], [a for _ in 1:Nₓ-1])
        B[1,2] = 0;
        B[1,1] = 1;
        B[end,end] = 1
        B[end,end-1] = 0

        # T_all = crank_nicolson_all(T0,A,B, Nₜ, Nₓ)
        # T = T_all[end,:]
        T = crank_nicolson(T0,A,B, Nₜ);
        res = betterϵ(T,Nₜ*Δt, xs, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    error_plot = plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for Crank-Nicolson scheme")
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_cn.pdf"))

    if disp
        display(error_plot)
    end

    return  ϵ_Ts
end;

function errorDevDF(disp=false)
    ϵ_Ts = Any[]
    res = 0

    Threads.@threads for Δt in Δts
        Nₜ = floor(Int,t/Δt)
        a = 2*λ*Δt/(Δx^2)
        T_all = dufort_frankel(T0,a,Nₜ,Nₓ)
        T = T_all[end,:]
        res = betterϵ(T,Nₜ*Δt, xs, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    error_plot = plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for Dufort-Frankel scheme")
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_df.pdf"))

    if disp
        display(error_plot)
    end

    return  ϵ_Ts

end
# @time task2a()

function composeErrors(fctsErr, ebErr, cnErr, dfErr, disp=false)
    data = [fctsErr, ebErr, cnErr, dfErr]
    labels = ["FCTS" "Euler Backward" "Crank-Nicolson" "Dufort-Frankel"]
    p =plot(Δts, data, label=labels)
    title!("Error comparison")
    
    xlabel!(L"Time resolution $\Delta t$")
    ylabel!(L"Error value $\epsilon(t=100)$")
    ylims!((0,0.00002))
    savefig(save_folder*"/error_comp_diffusion.pdf")

    if disp
        display(p)
    end
    
end



d=true
# @time err_fcts = errorDevFCTS(d);
# @time err_eb = errorDevEB(d);
# @time err_cn = errorDevCN(d);
# @time err_df = errorDevDF(d);

minimum(err_fcts)

# fullFCTS()
# fullEB()
# fullCN()
# fullDF()

# display(plot(Δts, err_cn))
# composeErrors(err_fcts, err_eb, err_cn, err_df, true);


