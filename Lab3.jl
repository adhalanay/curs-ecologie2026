### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ f4b392b0-b56c-11f0-3efc-5d06c6d0c331
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using DifferentialEquations
	using PlutoTeachingTools
	# using LoopVectorization
	TableOfContents()
end

# ╔═╡ 2535d5ed-642f-491f-9d39-8d904d21e9c9
begin
	total=82037000.0
	μ =770744.0
	δ=864330.0
	f(x)=(μ/total-δ/total)*x
	x0=1.0
    xs=[x0] 
for t ∈ range(1,50,length=50)
    x = f(last(xs))
    push!(xs,x)
end
fig1=Figure()
ax1 = Axis(fig1[1,1])
scatter!(total*xs,color=:red)
fig1
end

# ╔═╡ 807a6e1a-90d2-4d64-a1f4-b542f2c0811a
begin
	k=0.5
	b=0.4
	x1 = 1.0
	f1(x)=(1.0-k)*x+b
	    xs1=[x1] 
for t ∈ range(1,50,length=50)
    x = f1(last(xs1))
    push!(xs1,x)
end
fig2=Figure()
ax2 = Axis(fig2[1,1])
scatter!(xs1,color=:blue)
fig2
end

# ╔═╡ 673c91f2-6420-4a4f-9262-a5880d271571
begin
	M = 5000.0
	m = 100.0
	k1 = 0.0001
	f2(x)=x+k1*(M-x)*(x-m)
	x2 = 15001.0
	xs2=[x2]
	for t ∈ range(1,50,length=50)
    x = f2(last(xs2))
    push!(xs2,x)
end
fig3=Figure()
ax3 = Axis(fig3[1,1],limits=(nothing,(-x2,x2)))
scatter!(xs2,color=:green)
fig3
end

# ╔═╡ d8c09454-ba34-4cc6-b2f2-32c5057f9f05
begin
	r=3
f3(x)= r*x*(1-x)
f4(x)=x-f3(f3(x))
ts=LinRange(0.0,4.0,500)
ys=[f4(t) for t in ts]
fig4=Figure()
	ax4=Axis(fig4[1,1])
	scatter!(ys,markersize=5)
	fig4
end

# ╔═╡ 8e5e55d6-3345-4e5d-b8cb-96e2e220ac19
begin
	μ0=2
	function g(x) 
		if x<= 1/2 
			μ0*x
		else
			μ0*(1-x)
		end
	end
	ys1=[g(t) for t in ts]
  fig5=Figure()
	ax5=Axis(fig5[1,1])
	scatter!(ys1,markersize=5)
	fig5
end

# ╔═╡ Cell order:
# ╠═f4b392b0-b56c-11f0-3efc-5d06c6d0c331
# ╠═2535d5ed-642f-491f-9d39-8d904d21e9c9
# ╠═807a6e1a-90d2-4d64-a1f4-b542f2c0811a
# ╠═673c91f2-6420-4a4f-9262-a5880d271571
# ╠═d8c09454-ba34-4cc6-b2f2-32c5057f9f05
# ╠═8e5e55d6-3345-4e5d-b8cb-96e2e220ac19
