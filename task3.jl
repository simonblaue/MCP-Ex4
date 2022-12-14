
using Plots
using LaTeXStrings
using PlotThemes

using LinearAlgebra
using SparseArrays

# plotlyjs();
theme(:vibrant;)



#############

const ϵ = 0.2
const μ = 0.1
const Δx = 0.4
const Δt = 0.1
const L = 52

const xs = 0:Δx:L

const Nₓ = length(xs) 
const Nₜ = 2000

const t_end = Δt*Nₜ

###############

####### Iterative simulation ####

function getGroundtate(xs)

    u₁ = [1/2*(1-tanh((x-25)/5)) for x in xs]
    u₁[1] = 1
    u₁[2] = 1
    u₁[end] = 0
    u₁[end-1] = 0

    return u₁
end

function getMatricies()
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

    return A, B, C
    
end

function firstStepMatricies(u₁, ϵ, μ, Δt, Δx, A, B, C)

    u₂ = u₁ - (ϵ/6 * Δt/Δx) * (A * u₁) .* (B * u₁) - (μ/2 * Δt/Δx^3) * (C * u₁) 

    u₂[1] = 1
    u₂[end] = 0
    u₂[2] = u₁[2] - pref1 * (u₁[3]+u₁[2]+u₁[1]) * (u₁[3]-u₁[1] ) - pref2 * (u₁[4]+2*u₁[1]-2*u₁[3]-u₁[1])
    u₂[end-1] = u₁[end-1] - pref1 * (u₁[end]+u₁[end-1]+u₁[end-2]) * (u₁[end]-u₁[end-2] ) - pref2 * (u₁[end]+2*u₁[end-2]-2*u₁[end]-u₁[end-3])

    return u₂
end

function firstStep(u₁, ϵ, μ, Δt, Δx,Nₓ)
    pref1 = ϵ/6 * Δt/Δx
    pref2 = μ/2 * Δt/Δx^3

    u₂ = zeros(Nₓ)
    for j in 3:Nₓ-2
        u₂[j] = u₁[j] - pref1 * (u₁[j+1]+u₁[j]+u₁[j-1]) * (u₁[j+1]-u₁[j-1] ) - pref2 * (u₁[j+2]+2*u₁[j-1]-2*u₁[j+1]-u₁[j-2])
    end
    
    u₂[1] = 1
    u₂[end] = 0
    u₂[2] = u₁[2] - pref1 * (u₁[3]+u₁[2]+u₁[1]) * (u₁[3]-u₁[1] ) - pref2 * (u₁[4]+2*u₁[1]-2*u₁[3]-u₁[1])
    u₂[end-1] = u₁[end-1] - pref1 * (u₁[end]+u₁[end-1]+u₁[end-2]) * (u₁[end]-u₁[end-2] ) - pref2 * (u₁[end]+2*u₁[end-2]-2*u₁[end]-u₁[end-3])

    return  u₂
end

function matrixstep(uᵢ₋₁,uᵢ, ϵ, μ, Δt, Δx, A, B, C)

    uᵢ₊₁ = uᵢ₋₁ - ϵ/3 * Δt/Δx * (A*uᵢ) .* (B*uᵢ) - μ * Δt/Δx^3 * (C*uᵢ)

    uᵢ₊₁[1] = 1
    uᵢ₊₁[end] = 0
    uᵢ₊₁[2] = uᵢ₋₁[2] - pref1 * (uᵢ[3]+uᵢ[2]+uᵢ[1]) * (uᵢ[3]-uᵢ[1] ) - pref2 * (uᵢ[4]+2*uᵢ[1]-2*uᵢ[3]-uᵢ[1])
    uᵢ₊₁[end-1] = uᵢ₋₁[end-1] - pref1 * (uᵢ[end]+uᵢ[end-1]+uᵢ[end-2]) * (uᵢ[end]-uᵢ[end-2] ) - pref2 * (uᵢ[end]+2*uᵢ[end-2]-2*uᵢ[end]-uᵢ[end-3])

    return uᵢ₊₁
end;

function loopstep(uᵢ₋₁,uᵢ, ϵ, μ, Δt, Δx, Nₓ)

    pref1 = ϵ/3 * Δt/Δx
    pref2 = μ * Δt/Δx^3

    uᵢ₊₁ = zeros(Nₓ)
    for j in 3:Nₓ-2
        uᵢ₊₁[j] = uᵢ₋₁[j] - pref1 * (uᵢ[j+1]+uᵢ[j]+uᵢ[j-1]) * (uᵢ[j+1]-uᵢ[j-1] ) - pref2 * (uᵢ[j+2]+2*uᵢ[j-1]-2*uᵢ[j+1]-uᵢ[j-2])
    end
    
    uᵢ₊₁[1] = 1
    uᵢ₊₁[end] = 0
    uᵢ₊₁[2] = uᵢ₋₁[2] - pref1 * (uᵢ[3]+uᵢ[2]+uᵢ[1]) * (uᵢ[3]-uᵢ[1] ) - pref2 * (uᵢ[4]+2*uᵢ[1]-2*uᵢ[3]-uᵢ[1])
    uᵢ₊₁[end-1] = uᵢ₋₁[end-1] - pref1 * (uᵢ[end]+uᵢ[end-1]+uᵢ[end-2]) * (uᵢ[end]-uᵢ[end-2] ) - pref2 * (uᵢ[end]+2*uᵢ[end-2]-2*uᵢ[end]-uᵢ[end-3])

    return  uᵢ₊₁
end

function iterative()
    
    u = zeros(Nₜ, Nₓ)

    u[1,:] = getGroundtate(xs)
    u[2,:] = firstStep(u[1,:], ϵ, μ, Δt, Δx,Nₓ)

    ### evolution of u
    for i in 2:Nₜ-1
       u[i+1,:] = loopstep(u[i-1,:],u[i,:], ϵ, μ, Δt, Δx, Nₓ)
    end;

    #####
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
    n = size(u)[1]
    println(n)
    anim = @animate for i in 1:n
        t = round(i*Δt)
        plot(u[i,:], legend=false)
        ylims!(0,2)
        title!(L"t="*"$t")
        xlabel!("x")
        ylabel!("u")
    end

    display(gif(anim, "anim/$name.gif", fps=40))

end




# @time u= analytical(L ,t_end)
u = iterative()

anim(u[1:10:end,:], "iterative")

;


display(u[end,:])