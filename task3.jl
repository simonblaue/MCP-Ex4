
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
    xs = 0:Δx:L+Δx

    Nₓ = length(xs) 
    Nₜ = 2000
    n = 250

    u₁ = [1/2*(1-tanh((x-25)/5)) for x in xs]
    println(length(u₁))
    # plot(xs, u₀, legend=false)

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
    display(C)
    #### Done init

    ### First step

    u₂ = u₁ - (ϵ/6 * Δt/Δx) * (A * u₁) .* (B * u₁) - (μ/2 * Δt/Δx^3) * (C * u₁) 
    # u₂[2] =  - (ϵ/6 * Δt/Δx) * (A[2] * u₁) .* (B[2] * u₁) - (μ/2 * Δt/Δx^3) * (u₁[4] + 2*u₁[1] - 2*u₁[3]- u₁[1])
    # u₂[end-1] =  - (ϵ/6 * Δt/Δx) * (A[end-1] * u₁) .* (B[end-1] * u₁) - (μ/2 * Δt/Δx^3) * (u₁[end] + 2*u₁[end-2] - 2*u₁[3]- u₁[1])


    plot(xs, u₂, legend=false)

    ## 

    u = zeros(Nₜ, Nₓ)
    u[1,:] = u₁
    u[2,:] = u₂

    ### evolution of u


    for i in 2:Nₜ-1
        u[i+1,:] = u[i-1,:] - (ϵ/3*Δt/Δx) * (A * u[i,:]) .* (B * u[i,:]) - (μ*Δt/Δx^3) * (C * u[i,:])
        # u[i+1,1] = 1.
        # u[i+1,end] = 0.
        # u₂[2] = u₂[1]
        # u₂[end-1] = u₂[end] 
    end;
    #####

    u_save = [u[i,:] for i in 1:250:Nₜ]

    return u
        
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
    n = length(u[1,:])
    Nₓ = length(u[:,1])
    println(n)
    anim = @animate for i in 1:n
        plot(u[i,:], label="t=$i", dpi=800)
        xaxis!()
        title!(name)
        xlabel!("x")
        ylabel!("u")
    end

    display(gif(anim, "anim/$name.gif", fps=10))

end

L = 52
t_end = 200

# @time u= analytical(L ,t_end)
u = iterative()
display(plot(u[2,:]))
# anim(u, "iterative")
println("done")

