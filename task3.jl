
using Plots
using LaTeXStrings
using PlotThemes

using LinearAlgebra
using SparseArrays

plotlyjs();
theme(:vibrant;)

###############

function iterative()
    #### Var declaration

    ϵ = 0.2
    μ = 0.1
    Δx = 0.4
    Δt = 0.1
    L = 52
    xs = -Δx:Δx:L+Δx

    Nₓ = length(xs) 
    Nₜ = 2000

    u₁ = [1/2*(1-tanh((x-25)/5)) for x in xs]
    println(length(u₁))

    u₁[1] = 1
    u₁[2] = 1
    u₁[end] = 0
    u₁[end-1] = 0



    A = Tridiagonal(fill(1,Nₓ-1), fill(1,Nₓ), fill(1,Nₓ-1))
    A[1,:] .= 0
    A[2,:] .= 0
    A[end,:] .= 0
    A[end-1,:] .= 0

    B = Tridiagonal(fill(1,Nₓ-1), fill(0,Nₓ), fill(-1,Nₓ-1))
    B[1,:] .= 0
    B[2,:] .= 0
    B[end,:] .= 0
    B[end-1,:] .= 0

    C = zeros(Int, Nₓ,Nₓ)
    for i in 2:Nₓ-1
        C[i,i+1] = -2
        C[i,i-1] = 2
        if i>2
            C[i,i-2] = -1
        end
        if i<Nₓ-1
            C[i,i+2] = 1
        end
    end
    C[1,:] .= 0
    C[2,:] .= 0
    C[end,:] .= 0
    C[end-1,:] .= 0

    #### Done init

    ### First step

    u₂ = u₁ - (ϵ/6 * Δt/Δx) * (A * u₁) .* (B * u₁) - (μ/2 * Δt/Δx^3) * (C * u₁) 

    ## 
    u = zeros(Nₜ, Nₓ)
    u[1,:] = u₁
    u[2,:] = u₂

    ### evolution of u
    for i in 2:Nₜ-1
        u[i+1,:] = u[i-1,:] - (ϵ/3*Δt/Δx) * (A * u[i,:]) .* (B * u[i,:]) - (μ*Δt/Δx^3) * (C * u[i,:])
        u[1] = 1
        u[2] = 1
        u[end] = 0
        u[end-1] = 0
    end;

    #####
    return u[:,3:end-2]
        
end


function analytical(L, t_end)
    xs = 0:0.4:L
    ts = 0:0.1:t_end

    c = 0.5
    ξ₀ = 2

    u = Float64[-c/2 * sech(1/2*sqrt(c)*((x-c*t)-ξ₀))^2 for x in xs, t in ts ]

    return u
end

function anim(u, name)
    n = size(u)[1]
    println(n)
    anim = @animate for i in 1:n
        plot(u[i,:], label="t=$i")
        # ylims!(0,2)
        title!(name)
        xlabel!("x")
        ylabel!("u")
    end

    display(gif(anim, "anim/$name.gif", fps=30))

end

L = 52
t_end = 200

# @time u= analytical(L ,t_end)
u = iterative()
anim(u[1:10:end,:], "iterative")

;


display(u[end,:])