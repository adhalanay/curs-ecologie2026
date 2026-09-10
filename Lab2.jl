### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 4ec76850-aa67-11f0-3196-d9e3917b8038
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using DifferentialEquations
	TableOfContents()
end

# ╔═╡ 684aed4a-fe23-4694-a810-310173123f10
md"""
The population of an endangered species is declining at a rate of 2.5%ayear. Ifthere
 are currently 4000 individuals of this species, how many will there be in 20 years?
"""

# ╔═╡ d42ecacc-739a-4ecd-a933-aef6859b7df6
X(t)=4000*exp(-0.025*t)

# ╔═╡ da79387b-2357-49b6-ad1b-3649fcb63922
md"""
If money in a bank account earns an interest rate of 1.5%, compounded continuously,
 and the initial balance is $1000, how much money will be in the account in ten years?
"""

# ╔═╡ b4148ca9-be1a-4ae2-a76b-4680b376a482
begin
	t0=20
	print("Dupa $t0 ani vom avea $(X(t0)) indivizi")
end

# ╔═╡ fbb1b6ac-a864-4159-892d-8289e228eefe
begin
	t1 = 10
	X1(t) = 1000*exp(0.015*t)
	print("Dupa $t1 ani vom avea $(X1(t1)) dolari")
end

# ╔═╡ 9edc6339-9454-4227-b49d-d951cfc163af
md"""
Radioactive iodine, used to treat thyroid cancer, has a half-life of eight days. Find its decay constant, r."""

# ╔═╡ a5b5f252-529d-4ce4-9fa7-faa0522aef6d
begin
	r = -log(0.5)/8
	print("Rata de descompunere este $r")
end

# ╔═╡ 04c1b218-42b4-4c57-8de1-8f51aba13687
begin
	k=log(10000)/20
	C=exp(10*k)
	print("Mary trebuie sa cumpere $(round(C)) nuferi")
end

# ╔═╡ 95090872-88f2-435e-819d-117ef9929b95
md"""
You have \$10,000 and can put it either in an account bearing 3.9% interest compounded
 monthly or one bearing 4% interest compounded annually. If the money will be in the
 account for five years, which one should you choose?
"""

# ╔═╡ 71d62b4a-ee10-48ed-a22e-cdadff15f6f5
begin
	Xl(t)=10000*exp(0.039/12*t)
	Xa(t)=10000*exp(0.04*t)
	if Xl(60)-Xa(5) > 0.0
		print("Lunar e mai bine")
	else
		print("Anual e mai bine")
	end
	print(Xl(60)-Xa(5))
end

# ╔═╡ f94e23c5-411e-49dc-b0ea-638061b93056
begin
	t2 = range(-1.5,1.5,length=101)
	f(x)=2*x^2-x
	Xs = [f(x) for x in t2]  
	fig1=Figure()
    ax1=Axis(fig1[1,1],title="Graficul lui f")
    lines!(t2,Xs)
    fig1
end

# ╔═╡ 7ee2a7e1-e440-464f-ae62-df3db2de1f8b
begin

f1(u,p,t)=u*(1-1.5*u)*(1-u)
tspan = (0.0,80.0)
u0=1.1
prob1 = ODEProblem(f1,u0,tspan)
sol1 = solve(prob1, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
fig6=Figure()
ax6=Axis(fig6[1,1])
series!(sol1.t,hcat(sol1.u...),linewidth=5)
fig6
end

# ╔═╡ bc4cf710-61a1-455f-8551-bca427ce447a
md"""
Ecuatia lui von Bertalanffi
"""

# ╔═╡ 3db172ff-3714-4c5e-847d-6b4894eb36ee
begin

f2(u,p,t)=0.2*(5-u)
u1=4.6
prob2 = ODEProblem(f2,u0,tspan)
sol2 = solve(prob2, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
fig5=Figure()
ax5=Axis(fig5[1,1])
series!(sol2.t,hcat(sol2.u...),linewidth=5)
fig5
end

# ╔═╡ Cell order:
# ╠═4ec76850-aa67-11f0-3196-d9e3917b8038
# ╟─684aed4a-fe23-4694-a810-310173123f10
# ╠═d42ecacc-739a-4ecd-a933-aef6859b7df6
# ╟─da79387b-2357-49b6-ad1b-3649fcb63922
# ╠═b4148ca9-be1a-4ae2-a76b-4680b376a482
# ╠═fbb1b6ac-a864-4159-892d-8289e228eefe
# ╟─9edc6339-9454-4227-b49d-d951cfc163af
# ╟─a5b5f252-529d-4ce4-9fa7-faa0522aef6d
# ╠═04c1b218-42b4-4c57-8de1-8f51aba13687
# ╟─95090872-88f2-435e-819d-117ef9929b95
# ╠═71d62b4a-ee10-48ed-a22e-cdadff15f6f5
# ╠═f94e23c5-411e-49dc-b0ea-638061b93056
# ╠═7ee2a7e1-e440-464f-ae62-df3db2de1f8b
# ╠═bc4cf710-61a1-455f-8551-bca427ce447a
# ╠═3db172ff-3714-4c5e-847d-6b4894eb36ee
