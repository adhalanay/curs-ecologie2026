### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 0e204a80-baf0-11f0-22da-a308a465b653
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

# ╔═╡ 4fcf1efa-ed42-4531-8dc0-8480c1dc3d57
begin
	f1(x) = Point2f(24*x[2]-2*x[2]^2-3*x[1]*x[2],15*x[1]-x[1]^2-3*x[1]*x[2])
	fig1 = Figure()
	ax1 = Axis(fig1[1,1],title=L"x'=24x-2x^2-3xy,y=15x-y^2-3xy ",limits = (-0.1,15,-0.1,15))
	xs1 = LinRange(0, 15, 50)
    ys1 = LinRange(0, 15, 50)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax1,f1,xs1,ys1,colormap = :magma,linewidth = 1)
	scatter!(ax1,(0,12),markersize=20)
	scatter!(ax1, (15,0),markersize=20)
	scatter!(ax1,(0,0),markersize=20)
	scatter!(ax1,(6,3),markersize=20)
	# fig7
	n1(x) = 12.0-3.0/2.0*x
	n2(x) = 5.0-1.0/3.0*x
	ds = [n1(x) for x in xs1]
	ms = [n2(x) for x in xs1]
	lines!(ax1,[(0,0),(0,10)],linewidth=4)
	lines!(ax1,[(0,0),(10,0)],linewidth=4)
	lines!(ax1,xs1,ds,linewidth=4)
	lines!(ax1,xs1,ms,linewidth=4)
	fig1
end

# ╔═╡ 294c670f-3402-488d-b20e-2d81b76f9af0
begin
	f2(x) = Point2f(0.3*x[1]-0.02*x[1]^2-0.05*x[1]*x[2],0.2*x[2]-0.04*x[2]^2-0.02*x[1]*x[2])
	fig2 = Figure()
	ax2 = Axis(fig2[1,1],title=L"x'=24x-2x^2-3xy,y=15x-y^2-3xy ",limits = (-0.1,20,-0.1,20))
	xs2 = LinRange(0, 20, 50)
    ys2 = LinRange(0, 20, 50)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax2,f2,xs2,ys2,colormap = :magma,linewidth = 1)
	scatter!(ax2,(0,5),markersize=20)
	scatter!(ax2, (6,0),markersize=20)
	scatter!(ax2,(0,0),markersize=20)
	scatter!(ax2,(30/7,20/7),markersize=20)
	n3(x) = 10.0-5.0/3.0*x
	n4(x) = 5.0-1.0/2.0*x
	ds3 = [n3(x) for x in xs2]
	ms4 = [n4(x) for x in xs2]
	lines!(ax2,[(0,0),(0,20)],linewidth=4)
	lines!(ax2,[(0,0),(20,0)],linewidth=4)
	lines!(ax2,xs2,ds3,linewidth=4)
	lines!(ax2,xs2,ms4,linewidth=4)
	fig2

end

# ╔═╡ d288c832-41ec-49b0-94c5-c23dd5b5fd0b
begin
	f3(x) = Point2f(5/(1+x[2]^4)-x[1],5/(1+x[1]^4)-x[2])
	fr(c) = 5/(1+c^4)
	fc(r) = 5/(1+r^4)
	rs = LinRange(0,5,50)
	frs = [fc(r) for r in rs]
	cs = LinRange(0,5,50)
	fcs = [fr(c) for c in cs]
	fig3=Figure()
	ax3=Axis(fig3[1,1])
	lines!(cs,fcs)
	lines!(frs,rs)
	xs3 = LinRange(0, 10, 50)
    ys3 = LinRange(0, 10, 50)
	streamplot!(ax3,f3,xs3,ys3)
	fig3
end

# ╔═╡ Cell order:
# ╠═0e204a80-baf0-11f0-22da-a308a465b653
# ╠═4fcf1efa-ed42-4531-8dc0-8480c1dc3d57
# ╠═294c670f-3402-488d-b20e-2d81b76f9af0
# ╠═d288c832-41ec-49b0-94c5-c23dd5b5fd0b
