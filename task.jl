using Statistics
using Plots
using LaTeXStrings
using PlotThemes
plotlyjs();
include("integrator.jl")
theme(:vibrant;)


## Initialize
Δl = 0.005
l = 1

xs = 0:Δl:l
ys = 0:Δl:l

N = length(xs)

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
savefig("saves/comp_lapplace_series_heatmap.pdf")
display(p_all)
