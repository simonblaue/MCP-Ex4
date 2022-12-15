using Statistics
using Plots;
using LaTeXStrings
using PlotThemes

plotlyjs();
theme(:wong;)

include("integrator.jl")

## Initialize
Δl = 0.005
l = 1

xs = 0:Δl:l
ys = 0:Δl:l

N = length(xs)

ϕ₀ = zeros(length(xs),length(ys))
ϕ₀[1,:] .= 100;

save_folder="saves_t1"

### Functions

function task1a(xs,ys,l)
   ns = [1,10,100,1000]

   # setup plots
   plots = Any[]
   res = Any[]

   clims = (0,100)

   Threads.@threads for n in ns
      ϕ = [laplace_inf_sum(x,y,l,n) for x in xs, y in ys]
      push!(res,ϕ)
      push!(plots, heatmap(xs,ys,ϕ, title="n=$n", colorbar=false, clims=clims))
   end



   blank = plot(foreground_color_subplot=:white)
   l = @layout [grid(2, 2) a{0.075w}]
   plot(plots[1],plots[2],plots[3],plots[4], blank, layout=l, link=:all)
   p_all = scatter!([0], [0], zcolor=[NaN], clims=clims, label="",background_color_subplot=:transparent, markerstrokecolor=:transparent, framestyle=:none, inset=bbox(0.1, 0, 0.6, 0.9, :center, :right), subplot=6)
   savefig(save_folder*"/comp_lapplace_series_heatmap.pdf")
end

function task1b(ϕ₀)
   x_j, k_j, ϵs_max_j,ϵs_av_j = jacobi(ϕ₀);
   x_gs, k_gs, ϵs_max_gs ,ϵs_av_gs = gauß_seidel(ϕ₀);
   αs = [0.5 1.0 1.25 1.5 1.75 1.99]
   res = Any[]

   for α in αs
      push!(res, SOR(ϕ₀, α))
   end

   # Plot maximal error
   data = [ϵs_max_j, ϵs_max_gs,res[1][3],res[3][3],res[4][3],res[5][3], res[6][3]]
   labels = ["Jacobi" "Gauss-Seidel" "SOR 0.5" "SOR 1.25" "SOR 1.5" "SOR 1.75" "SOR 1.99"]
   plot(data, xlim=(0,100), label=labels, title="Maximal error vs step", xaxis=:log, yaxis=:log, xticks=([1, 10, 100, 1000, 10000], [1, 10, 100, 1000, 10000]))
   xlabel!("Step")
   ylabel!("Maximal error")
   savefig(string(save_folder, "/max_errors_comp.pdf"))

   # Plot av error
   data = [ϵs_av_j, ϵs_av_gs, res[1][4],res[2][4],res[3][4],res[4][4],res[5][4],res[6][4]]
   labels = ["Jacobi" "Gauß Seidel" "SOR 0.5" "SOR 1.0" "SOR 1.25" "SOR 1.5" "SOR 1.75" "SOR 1.99"]
   plot(data, xlim=(0,400), label=labels, title="Average error vs step", xaxis=:log, yaxis=:log, xticks=([1, 10, 100, 1000, 10000], [1, 10, 100, 1000, 10000]))
   xlabel!("Step")
   ylabel!("Average error")
   savefig(string(save_folder, "/av_errors_comp.pdf"))

   #Plot SOR 1.99 heatmap
   heatmap(xs,ys, res[6,1], title="Domain after SOR 1.99 method")
   savefig(string(save_folder, "/sor199_heatmap.pdf"))

   # Plot convegence steps
   data = [k_j, k_gs, res[1][2],res[2][2],res[3][2],res[4][2],res[5][2],res[6][2]]
   labels = ["Jacobi", "Gauß Seidel", "SOR 0.5", "SOR 1.0", "SOR 1.25", "SOR 1.5", "SOR 1.75", "SOR 1.99"]
   scatter(labels,data, title="Number of steps until convergence", legend=false)
   xlabel!("Method")
   ylabel!("Steps until convergence")
   savefig(string(save_folder, "/number_of_convergence_steps.pdf"))

   # Difference plot
   iterative_heatmap = heatmap(xs,ys, x_gs, title="After Gauß Seidel Iteration")
   series_heatmap_1000 = heatmap(xs, ys, [laplace_inf_sum(x,y,l,1000) for x in xs, y in ys], title="Infinite sum with n=1000")
   difference = x_gs-[laplace_inf_sum(x,y,l,1000) for x in xs, y in ys]
   difference_heatmap = heatmap(xs,ys, difference, title="Difference", clims=(minimum(difference),maximum(difference[2:end,2:end])))
   savefig(string(save_folder, "/difference_laplace_heatmap.pdf"))

   # heatmap comp
   plot(iterative_heatmap, series_heatmap_1000, layout=grid(1,2, widths=(4/8,4/8)), size=(800,300))
   savefig(string(save_folder, "/comp_laplace_heatmap.pdf"))
end

function task1c(ϕ₀)

   broke_SOR = SOR(ϕ₀, 2.0)

   heatmap(xs, ys, broke_SOR[1], title="SOR for α=2.0")
   savefig(save_folder * "/broken_SOR_heatmap.pdf")

   plot(broke_SOR[3],title="Maximal error vs step", label="SOR for α = 2.0", xaxis=:log, yaxis=:log, xticks=([1, 10, 100, 1000, 10000], [1, 10, 100, 1000, 10000]))
   xlabel!("Timesteps")
   ylabel!("Max error ϵ")
   savefig(save_folder * "/broken_SOR_error.pdf")

end

# task1a(xs,ys,l)
task1b(ϕ₀)
# task1c(ϕ₀)s

