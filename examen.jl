### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 353516e2-f5f5-11f0-24a3-69addd20ae34
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using DifferentialEquations
	using PlutoTeachingTools
	using DynamicalSystems
	# using LoopVectorization
	TableOfContents()
end

# ╔═╡ 04d5bc60-11cf-481e-ab7b-40aa784458a9
@bind μ PlutoUI.Slider(0.5:-0.01:-0.5,show_value=true)

# ╔═╡ 29f56919-83d5-4105-904d-8f9a4a4cb2de
begin
	xs = LinRange(-2.0, 2.0, 500)
    ys= LinRange(-2.0, 2.0, 500)
	s(x)=Point2f(μ-x[1]^2,-x[2])
	t1(x)=Point2f(μ*x[1]-x[1]^2,-x[2])
	t2(x)=Point2f(μ*x[1]-x[1]^3,-x[2])
	h(x) = Point2f(-x[2]+μ*x[1]-x[1]*x[2]^2,x[1]+μ*x[2]-x[1]^2)
	fig1 = Figure()
		ax = Axis(fig1[1,1])
		streamplot!(ax,h,xs,ys,stepsize = 0.01,gridsize=(60,60))
		# scatter!(ax,[-sqrt(μ),sqrt(μ)],[0,0],markersize=10,color=:green)
		Makie.scatter!(ax,[0],[0],markersize=10,color=:green)
	fig1
end

# ╔═╡ 12244e57-be2b-497c-9d97-c34e5884a393
@bind c PlutoUI.Slider(-2.0:0.1:2.0,show_value=true)

# ╔═╡ 21bb679e-0e20-45a8-a9fc-ca490c6c2887
begin
	function hopf!(du,u,p,t)
		p=[c]
		du[1]=p[1]*u[1]+u[2]
		du[2]=-u[1]+p[1]*u[2]-u[1]^2*u[2]
		return nothing
	end
	us=[[0.0,0.0],[-1.0,1.0],[-1.0,-1.0],[1.0,1.0]]
	t_final=20.0
	s_t = 0.01
	sols=[]
	probs=[]
	for u ∈ us 
		prob=CoupledODEs(hopf!,u)
		push!(probs,prob)
	end
	fig2=Figure()
	ax2 = Axis(fig2[1,1])
	for prob ∈ probs
		sol,t = trajectory(prob,t_final;Ttr=2.2, Δt=s_t)
		scatter!(ax2,sol,colormap=:thermal)
	end
	fig2
end

# ╔═╡ Cell order:
# ╠═353516e2-f5f5-11f0-24a3-69addd20ae34
# ╠═04d5bc60-11cf-481e-ab7b-40aa784458a9
# ╠═29f56919-83d5-4105-904d-8f9a4a4cb2de
# ╠═12244e57-be2b-497c-9d97-c34e5884a393
# ╠═21bb679e-0e20-45a8-a9fc-ca490c6c2887
