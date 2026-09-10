### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 48136d93-2bd9-4a02-bd6b-a6f1d605c3bf
begin
using Pkg
Pkg.activate(".")
end;

# ╔═╡ 36e85730-9603-11f0-1552-d1341786a12b
using PlutoUI, GLMakie

# ╔═╡ 871673ac-b1ef-42bf-96ca-a956348c3693
html"<button onclick='present()'>present</button>"

# ╔═╡ ef0cd7a1-13c0-4e3a-baa3-62ecaafa037d
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ a227f331-c1a1-49b1-a7e8-0fbd9d47bd3f
md"""
$$\int_0^\infty e^{-x^2}dx=\sqrt{\pi}$$
"""

# ╔═╡ 5a52aea3-2418-437a-9115-1a5ffbc9862b
@bind r PlutoUI.Slider(0.01:0.01:4.0, show_value=true)

# ╔═╡ cec2384d-50ee-4963-91c3-3df516d5fc5d
f(x)=r*x*(1-x);

# ╔═╡ ee9106b0-70b9-4971-bed7-970617fc5505
begin
	n=400
	x0=0.7
	xs=[]
	for i∈1:n
		x1=f(x0)
		push!(xs,x1)
		x0=x1
	end
end

# ╔═╡ a3505a12-fe09-4dc1-8bc6-fc3b6fabf28a
begin
	ys = collect(1:n)
	fig = Figure(size=(1200,600))
		ax = Axis(fig[1,1])
		pts = Point2f.(ys,xs)
	    scatter!(pts,color = 1:n, markersize=20)
	fig
end

# ╔═╡ Cell order:
# ╟─871673ac-b1ef-42bf-96ca-a956348c3693
# ╟─ef0cd7a1-13c0-4e3a-baa3-62ecaafa037d
# ╟─48136d93-2bd9-4a02-bd6b-a6f1d605c3bf
# ╟─36e85730-9603-11f0-1552-d1341786a12b
# ╟─a227f331-c1a1-49b1-a7e8-0fbd9d47bd3f
# ╠═cec2384d-50ee-4963-91c3-3df516d5fc5d
# ╠═ee9106b0-70b9-4971-bed7-970617fc5505
# ╟─a3505a12-fe09-4dc1-8bc6-fc3b6fabf28a
# ╠═5a52aea3-2418-437a-9115-1a5ffbc9862b
