
############### Task 1

function ϵ(ϕ)
  return abs.([ϕ[i,j] - (ϕ[i+1,j] + ϕ[i-1,j] + ϕ[i,j+1] + ϕ[i,j-1])/4  for i in 2:length(ϕ[:,1])-1, j in 2:length(ϕ[1,:])-1])
end;


function laplace_inf_sum(x,y,l,terms) 
  summ = 0
  for n in 1:2:2*terms
    summ += 1/n * sin(n*π*y/l) * exp(-n*π*x)
  end
  return  400/π *summ
end;


function jacobi(ϕ0)
  ϕ = copy(ϕ0)
  k = 0
  N = length(ϕ[:,1])
  ϵs_max = Any[]
  ϵs_av = Any[]
  while( maximum(ϵ(ϕ)) > 1e-3 && k<1e5)
    ϕ_old = (ϕ)
    for i in 2:N-1
      for j in 2:N-1
        ϕ[i,j] = (ϕ_old[i+1,j] + ϕ_old[i-1,j] + ϕ_old[i,j+1] + ϕ_old[i,j-1])/4
      end
    end
    k += 1
    push!(ϵs_max, maximum(ϵ(ϕ)))
    push!(ϵs_av, mean(ϵ(ϕ)))
  end

  end_err = maximum(ϵ(ϕ))
  println("Stopped after $k steps with an error of $end_err ")
  return ϕ,k, ϵs_max, ϵs_av
end;

function gauß_seidel(ϕ0)
  ϕ = copy(ϕ0)
  k = 0
  N = length(ϕ[:,1])
  ϵs_max = Any[]
  ϵs_av = Any[]
  while( maximum(ϵ(ϕ)) > 1e-3 && k<1e5)
    for i in 2:N-1
      for j in 2:N-1
        ϕ[i,j] = (ϕ[i+1,j] + ϕ[i-1,j] + ϕ[i,j+1] + ϕ[i,j-1])/4
      end
    end
    k += 1
    push!(ϵs_max, maximum(ϵ(ϕ)))
    push!(ϵs_av, mean(ϵ(ϕ)))
  end

  end_err = maximum(ϵ(ϕ))
  println("Stopped after $k steps with an error of $end_err ")
  return ϕ,k, ϵs_max, ϵs_av
end;

function SOR(ϕ0, ω)
  ϕ = copy(ϕ0)
  k = 0
  N = length(ϕ[:,1])
  ϵs_max = Any[]
  ϵs_av = Any[]
  while( maximum(ϵ(ϕ)) > 1e-3 && k<1e5)
    for i in 2:N-1
      for j in 2:N-1
        ϕ[i,j] = (1-ω)*ϕ[i,j] + ω*(ϕ[i+1,j] + ϕ[i-1,j] + ϕ[i,j+1] + ϕ[i,j-1])/4
      end
    end
    k += 1
    push!(ϵs_max, maximum(ϵ(ϕ)))
    push!(ϵs_av, mean(ϵ(ϕ)))
  end

  end_err = maximum(ϵ(ϕ))
  println("Stopped after $k steps with an error of $end_err ")
  return ϕ,k, ϵs_max, ϵs_av
end;

function get_A_laplace()
  A = zeros(N*N,N*N)
  for r in 1:N*N
    A[r,r] = -4
    if r>1
      A[r,r-1] = 1
    end
    if r>N
      A[r,r-N] = 1
    end
    if r<N*N-N
      A[r,r+N] = 1
    end
    if r < N*N
      A[r,r+1] = 1
    end
  end;
  return A
end;

################ Task 2

## Integrators 

function FCTS(T0, Δt, Δx, λ, Nₜ, Nₓ)
  a = λ*Δt/(Δx^2)
  T = copy(T0)
  Teval = zeros(Nₓ, Nₜ+1)
  Teval[:,1] = copy(T)

  for i in 1:Nₜ
    T_old = copy(T) 
    for j in 2:length(T)-1
      T[j] = (1-2*a)*T_old[j] + a*(T_old[j-1]+T[j+1])  
    end
    Teval[:,i+1] = copy(T)
  end
  return Teval'
end;

function euler_backward(T0, A, Nₜ, Nₓ)
  A⁻¹ = inv(A)
  T = zeros(Nₓ, Nₜ+1)
  T[:,1] = copy(T0)
  for i in 1:Nₜ
      T[:,i+1] = A⁻¹ * vec(T[:,i])
  end
  return T'
end;

function crank_nicolson(T0, A, B, Nₜ, Nₓ)
  A⁻¹ = inv(A)
  T = zeros(Nₓ, Nₜ+1)
  T[:,1] = copy(T0)
  for i in 1:Nₜ
      T[:,i+1] = A⁻¹ * B * vec(T[:,i])
  end
  return T'
end;

function dufort_frankel(T0,a, Nₜ, Nₓ)
  T = zeros(Nₓ, Nₜ+1)
  T[:,1] = copy(T0)
  T[:,2] = copy(T0)
  A = Diagonal([(1-a)/(1+a) for _ in 1:Nₓ])
  B = Tridiagonal([a/(1+a) for _ in 1:Nₓ-1], [0. for _ in 1:Nₓ], [a/(1+a) for _ in 1:Nₓ-1])
  for i in 2:Nₜ
    T[:,i+1] = A*T[:,i-1] + B * T[:,i]
  end
  return T'
end;

## Error calcultaion

function T_exact(x,t, L,K,C,ρ)
  return sin((π*x)/L)*exp(-(π^2*K*t)/(L^2*C*ρ))
end;


function ϵ_T(Teval,t, Δx, Δt, L,K,C,ρ)
  Nₓ = floor(Int,L/Δx)
  T = (x,t) -> Teval[floor(Int,t/Δt),floor(Int,x/Δx)]
  summe = sum([abs(T(j*Δx,t) - T_exact(j*Δx,t,L,K,C,ρ)) for j in 1:Nₓ-1])
  return 1/Nₓ * summe
end;