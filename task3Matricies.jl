
function getMatricies()
    A = Tridiagonal(fill(1,Nₓ-1), fill(1,Nₓ), fill(1,Nₓ-1))

    B = Tridiagonal(fill(1,Nₓ-1), fill(0,Nₓ), fill(-1,Nₓ-1))


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

    A[1,1] = 0
    A[end,end] = 0
    B[1,1] = 0
    B[end,end] = 0
    C[1,1] = 0
    C[end,end] = 0

    return A, B, C
    
end

function firstStepMatricies(u₁, ϵ, μ, Δt, Δx, A, B, C)

    pref1 = (ϵ/6 * Δt/Δx)
    pref2 = (μ/2 * Δt/Δx^3)

    u₂ = u₁ - pref1 * (A * u₁) .* (B * u₁) - pref2 * (C * u₁) 

    u₂[1] = 1
    u₂[end] = 0
    u₂[2] = u₁[2] - pref1 * (u₁[3]+u₁[2]+u₁[1]) * (u₁[3]-u₁[1] ) - pref2 * (u₁[4]+2*u₁[1]-2*u₁[3]-u₁[1])
    u₂[end-1] = u₁[end-1] - pref1 * (u₁[end]+u₁[end-1]+u₁[end-2]) * (u₁[end]-u₁[end-2] ) - pref2 * (u₁[end]+2*u₁[end-2]-2*u₁[end]-u₁[end-3])

    return u₂
end

function matrixstep(uᵢ₋₁,uᵢ, ϵ, μ, Δt, Δx, A, B, C)

    pref1 = ϵ/3 * Δt/Δx
    pref2 = μ * Δt/Δx^3

    uᵢ₊₁ = uᵢ₋₁ - pref1 * (A*uᵢ) .* (B*uᵢ) - pref2 * (C*uᵢ)

    uᵢ₊₁[1] = 1
    uᵢ₊₁[end] = 0
    uᵢ₊₁[2] = uᵢ₋₁[2] - pref1 * (uᵢ[3]+uᵢ[2]+uᵢ[1]) * (uᵢ[3]-uᵢ[1] ) - pref2 * (uᵢ[4]+2*uᵢ[1]-2*uᵢ[3]-uᵢ[1])
    uᵢ₊₁[end-1] = uᵢ₋₁[end-1] - pref1 * (uᵢ[end]+uᵢ[end-1]+uᵢ[end-2]) * (uᵢ[end]-uᵢ[end-2] ) - pref2 * (uᵢ[end]+2*uᵢ[end-2]-2*uᵢ[end]-uᵢ[end-3])

    return uᵢ₊₁
end;

function iterativeMatricies(Nₜ, Nₓ, xs, ϵ, μ, Δt, Δx)
    u = zeros(Nₜ, Nₓ)

    u[1,:] = getGroundtate(xs)
    A,B,C = getMatricies()
    u[2,:] = firstStepMatricies(u[1,:], ϵ, μ, Δt, Δx,A,B,C)

    ### evolution of u
    for i in 2:Nₜ-1
       u[i+1,:] = matrixstep(u[i-1,:],u[i,:], ϵ, μ, Δt, Δx, A,B,C)
    end;

    #####
    return u
    
end