
using Plots
using LaTeXStrings
using PlotThemes

using LinearAlgebra
using SparseArrays
gr();
# plotlyjs();
theme(:vibrant;)


include("task3Matricies.jl")

#############

const ϵ = 0.2
const μ = 0.1
const Δx = 0.4
const Δt = 0.05
const L = 52

const xs = 0:Δx:L
const Nₓ = length(xs) 

const Nₜ = 2000
const ts = 0:Δt:Nₜ

const t_end = Δt*Nₜ

const save_folder = "saves_t3/"

###############

####### Iterative simulation ####

function getInitSoliton(xs)

    u₁ = [1/2*(1-tanh((x-25)/5)) for x in xs]
    return u₁
end

function getCollidingSoliton(xs)
    u₁ = [ 0.8*(1-tanh(3*x/12 -3)^2)+ 0.3(1-tanh(4.5*x/26 - 4.5)^2)  for x in xs]
    return u₁
end

function firstStep(u₁, ϵ, μ, Δt, Δx,Nₓ,boundarys)
    pref1 = ϵ/6 * Δt/Δx
    pref2 = μ/2 * Δt/Δx^3

    u₂ = zeros(Nₓ)
    for j in 3:Nₓ-2
        u₂[j] = u₁[j] - pref1 * (u₁[j+1]+u₁[j]+u₁[j-1]) * (u₁[j+1]-u₁[j-1] ) - pref2 * (u₁[j+2]+2*u₁[j-1]-2*u₁[j+1]-u₁[j-2])
    end
    
    u₂[1] = boundarys[1]
    u₂[end] = boundarys[end]
    u₂[2] = u₁[2] - pref1 * (u₁[3]+u₁[2]+u₁[1]) * (u₁[3]-u₁[1] ) - pref2 * (u₁[4]+2*u₁[1]-2*u₁[3]-u₁[1])
    u₂[end-1] = u₁[end-1] - pref1 * (u₁[end]+u₁[end-1]+u₁[end-2]) * (u₁[end]-u₁[end-2] ) - pref2 * (u₁[end]+2*u₁[end-2]-2*u₁[end]-u₁[end-3])

    return  u₂
end


function loopstep(uᵢ₋₁,uᵢ, ϵ, μ, Δt, Δx, Nₓ, boundarys)

    pref1 = ϵ/3 * Δt/Δx
    pref2 = μ * Δt/Δx^3

    uᵢ₊₁ = zeros(Nₓ)
    for j in 3:Nₓ-2
        uᵢ₊₁[j] = uᵢ₋₁[j] - pref1 * (uᵢ[j+1]+uᵢ[j]+uᵢ[j-1]) * (uᵢ[j+1]-uᵢ[j-1] ) - pref2 * (uᵢ[j+2]+2*uᵢ[j-1]-2*uᵢ[j+1]-uᵢ[j-2])
    end
    
    uᵢ₊₁[1] = boundarys[1]
    uᵢ₊₁[end] = boundarys[end]
    uᵢ₊₁[2] = uᵢ₋₁[2] - pref1 * (uᵢ[3]+uᵢ[2]+uᵢ[1]) * (uᵢ[3]-uᵢ[1] ) - pref2 * (uᵢ[4]+2*uᵢ[1]-2*uᵢ[3]-uᵢ[1])
    uᵢ₊₁[end-1] = uᵢ₋₁[end-1] - pref1 * (uᵢ[end]+uᵢ[end-1]+uᵢ[end-2]) * (uᵢ[end]-uᵢ[end-2] ) - pref2 * (uᵢ[end]+2*uᵢ[end-2]-2*uᵢ[end]-uᵢ[end-3])

    return  uᵢ₊₁
end

function iterative(u₁, boundarys=[1.,0.])
    
    u = zeros(Nₜ, Nₓ)

    u[1,:] = u₁
    u[2,:] = firstStep(u[1,:], ϵ, μ, Δt, Δx,Nₓ,boundarys)

    ### evolution of u
    for i in 2:Nₜ-1
       u[i+1,:] = loopstep(u[i-1,:],u[i,:], ϵ, μ, Δt, Δx, Nₓ,boundarys)
    end;

    #####
    return u
        
end

function analytical(L, t_end)
    xs = 0:0.4:L
    ts = 0:0.1:t_end

    c = 0.5
    ξ₀ = 2

    u = [-c/2 * sech(1/2*sqrt(c)*((x-c*t)-ξ₀))^2 for x in xs, t in ts ]

    return u
end


function anim(u,every, name)
    u_anim = u[1:every:end, :]
    n = size(u_anim)[1]
    anim = @animate for i in 1:n
        t = i*every
        plot(xs, u_anim[i,:], legend=false)
        ylims!(0,2)
        title!(L"t="*"$t")
        xlabel!("x")
        ylabel!("u")
    end

    display(gif(anim, "anim/$name.gif", fps=40))

end

function d3plot(u)
    plotlyjs()
    plt = plot(u,st=:surface,camera=(-30,30))

    display(plt)
end;

function timeEvplot(u, savename,every=250, disp=false)
    list_plots = Any[]
    plt = plot(xs, u[1,:], title="Time evolution", label="Timestep 1",xformatter=_->"")
    push!(list_plots,plt)
    for i in every+1:every:Nₜ-every
        plt = plot(xs, u[i,:], label="Timestep $i", xformatter=_->"" )
        # ylabel!("Disturbance")
        # title!("Timestep $i")
        push!(list_plots,plt)
    end
    plt = plot(xs, u[Nₜ,:], label="Timestep $Nₜ", xlabel=L"Domain $x$")
    push!(list_plots,plt)


    final_plot = plot(
        list_plots...,
        layout = (length(list_plots), 1),
        size=(600,1200), 
        ylims=(-0.2,1.35),
        link=:both,
        ylabel=L"Disturbance $u$")
    savefig(save_folder*savename*".pdf")
    if disp
        display(final_plot)
    end
    # display(p)

end

# @time u= analytical(L ,t_end)
# u₁ = getCollidingSoliton(xs)
# @time u = iterative(u₁,[0.,0.])

u₁ = getInitSoliton(xs)
@time u = iterative(u₁,[1.,0.])

# @time u = iterativeMatricies(Nₜ, Nₓ, xs, ϵ, μ, Δt, Δx);

println("calculated u waiting for animation...")
anim(u,10, "singleSoliton");
# d3plot(u)
# timeEvplot(u, "colliding_solitons",1500, true)
println("Done ⚡️")

