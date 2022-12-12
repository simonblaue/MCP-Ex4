using Statistics
using Plots
using LaTeXStrings
using PlotThemes
using LinearAlgebra
plotlyjs();
theme(:vibrant;)

include("integrator.jl")

L = 1
C = 900
K = 210
ρ = 2700
Δt = 0.1
Δx = 0.01
Nₜ = 10000
Nₓ = floor(Int,L/Δx+1)

x = 0:Δx:L

λ = K/(C*ρ)
a = λ*Δt/Δx^2
T0 = sin.((π*x) / L)

save_folder = "saves_t2"

function task2a()

    t =  0:Nₜ
    Teval = FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ);
    hmap_ftcs = heatmap(x, t.*Δt, Teval, ylabel="Time", title="Rod evolution with FTCS")
    savefig(string(save_folder,"/rod_FTCS.pdf"))

    t = 100
    Δts = 0.0001:0.0005:0.7
    ϵ_Ts = Any[]
    res = 0

    Threads.@threads for Δt in Δts
        Nₜ = floor(Int,t/Δt+2)
        Teval = FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ);
        res = ϵ_T(Teval,t, Δx,Δt, L,K,C,ρ)
        push!(ϵ_Ts, res)
    end

    plot(Δts, ϵ_Ts, legend=false)
    title!("Error development for FCTS scheme")
    xlabel!("Time resolution Δt")
    ylabel!("Error ϵ(t=100)")
    savefig(string(save_folder,"/error_development_fcts.pdf"))

end;

task2a()