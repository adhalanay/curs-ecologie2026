### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 3288dd8a-c3d1-11f0-2899-a34f81731429
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

# ╔═╡ f1c3504d-7536-48b2-9569-eeeb497135d7
@bind a PlutoUI.Slider(0.1:0.01:1.8,show_value=true)

# ╔═╡ 58308c55-0d54-4725-bc3c-03aaaa9163ea
@bind b PlutoUI.Slider(0.2:0.01:0.7,show_value=true)

# ╔═╡ fc719755-95c0-4f62-97d1-a6e0bf7caf38
begin
	function henon_rule(u,p,n)
		x, y = u # system state
        a, b = p # system parameters
        xn = 1.0 - a*x^2 + y
        yn = b*x
        return SVector(xn, yn)
	end
	p0 = [a, b]
	init_vals = [[0.2,0.3],[0.5,0.3],[0.1,0.2],[0.3,0.1]]
	u0 = [0.2, 0.3]
	prbs= []
	total_time = 1000
	for u in init_vals
		henon = DeterministicIteratedMap(henon_rule, u, p0)
    	X, t = trajectory(henon, total_time)
		push!(prbs,(X,t))
	end
	fig1 = Figure()
		ax1=Axis(fig1[1,1])
		ax2=Axis(fig1[1,2])
		for (X,t) in prbs
			scatter!(ax1,X,markersize=10)
			for var in columns(X)
				scatter!(ax2,t,var)
			end
		end
	fig1
end

# ╔═╡ Cell order:
# ╠═3288dd8a-c3d1-11f0-2899-a34f81731429
# ╠═f1c3504d-7536-48b2-9569-eeeb497135d7
# ╠═58308c55-0d54-4725-bc3c-03aaaa9163ea
# ╠═fc719755-95c0-4f62-97d1-a6e0bf7caf38
