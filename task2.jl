using Statistics
using Plots
using LaTeXStrings
using PlotThemes
using LinearAlgebra
plotlyjs();
theme(:default;)

include("integrator.jl")

L = 1
C = 900
K = 210
ρ = 2700
Δt = 0.1
Δx = 0.01
Nₜ = 10000
Nₓ = floor(Int,L/Δx+1)

Δts = 0.0001:0.0005:0.7
ts =  0:Nₜ
t=100

xs = 0:Δx:L

λ = K/(C*ρ)
a = λ*Δt/Δx^2
T0 = sin.((π*x) / L)

save_folder = "saves_t2"

function task2a()


    Teval = FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ);
    hmap_ftcs = heatmap(xs, ts.*Δt, Teval, ylabel="Time",xlabel="Rod domain", title="Rod evolution with FTCS")
    savefig(save_folder*"/rod_FTCS.pdf")


    # Δts = 0.6:0.0005:1
    ϵ_Ts = Any[]
    res = 0

    Threads.@threads for Δt in Δts
        Nₜ = floor(Int,t/Δt+2)
        Teval = FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ);
        res = ϵ_T(Teval,t, Δx,Δt, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    error_plot = plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for FCTS scheme")
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_fcts.pdf"))

    display(error_plot)
end;

function errorDevEB()
    ϵ_Ts = Any[]
    res = 0

    Threads.@threads for Δt in Δts
        Nₜ = floor(Int,t/Δt+2)
        a = λ*Δt/Δx^2
        A =Tridiagonal([-a for _ in 1:Nₓ-1], [(1+2*a) for _ in 1:Nₓ], [-a for _ in 1:Nₓ-1])
        A[1,2] = 0;
        T = euler_backward(T0,A, Nₜ, Nₓ);
        res = betterϵ(T,t, xs, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    error_plot = plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for Euler Backward scheme")
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_eb.pdf"))

    display(error_plot)

    return  ϵ_Ts
end;

# @time task2a()
@time err_eb = errorDevEB();

