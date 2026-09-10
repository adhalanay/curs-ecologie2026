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

# ╔═╡ e5b1174c-d529-11f0-2a30-6f5ef7b91715
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

# ╔═╡ 522e922b-287c-4c34-87cc-b29bcca0cd40
begin
	using GLMakie
	h(x,y) = log(x)-x-y+log(y)
	f(x) = Point2f(x[1]*x[2]-x[1],-x[1]*x[2]+x[2])
	S = LinRange(0.5, 1.5, 100)
    T = LinRange(0.5, 1.5, 100)
	H = [h(s,t) for s ∈ S, t ∈ T]
	fig1=Figure(size=(1000,500))
	ax1 = Axis3(fig1[1,1])
	surface!(ax1,S,T,H)
	scatter!(ax1,[1],[1],[h(1,1)],markersize = 20.0,color=:red)
	ax2 = Axis(fig1[1,2])
	streamplot!(ax2,f,S,T)
	fig1
end

# ╔═╡ 11d49c9f-f3d6-4d1b-bfd4-08dd1491cc97
WidthOverDocs()

# ╔═╡ 7f2464da-120e-4774-b364-5ed61eb97272
md"""
# Oscilații în sisteme dinamice

## Sisteme conservative

O clasă de sisteme care posedă centre este ce a **sistemelor conservative**.Acestea se caracterizează prin prezența unei funcții ``H`` constante de-a lungul fiecărei traiectorii, adică ``\frac{dH}{dt}.`` Astfel de sisteme nu pot avea echilibre stabile și nici cicluri limită. Punctele de echilibru pot fi doar centre sau uncte șa. Chiar dacă sistemele conservative nu sînt structural stabile sînt totuși utile în unele situații. Pentru ele putem demonstra existența centrelor, chiar dacă nu se poate aplica 
teorema Hartmann-Grobmann. Avem următoare:

**Teoremă** Fie ``V(x,y)`` un cîmp de vectori ``2-``dimensional și ``(x^*,y^*)``un echilibru izolat. Dacă sistemul este conservativ, adică există o funcție ``H(x,y)`` constantă de-a lungul traiectoriilot și ``(x^*,y^*)`` este un punct de extrem local pentru ``H``, atunci ``(x^*,y^*)`` este un centru și toate orbitele din jurul lui sînt închise.

Situația este reprezentată în imaginea de mai jos

$(PlutoUI.LocalResource("./COURSE_NOTE_6_Conservative_and_gradient_systems_8.png",:width=>500))

(descărcată de [aici](https://venturi.soe.ucsc.edu/sites/default/files/COURSE_NOTE_6_Conservative_and_gradient_systems.pdf)).

Revenim la sistemul Lotka-Volterra:
```math
\begin{aligned}
S' & = & ST & - & S \\
T' & = & -ST & + & T 
\end{aligned}
```
Acest sistem are două puncte fixe ``(0,0)`` și respectiv ``(1,1)``. Primul este un punct șa. Pentru cel de-al doilea valorile proprii ale Jacobienei sînt ``\pm i,`` deci nu se aplică teorema Hartmann-Grobmann. Dacă scriem forma generală a sistemului 
```math
\begin{aligned}
S' & = & aST & - & dS \\
T' & = & -dST & + & cT 
\end{aligned}
```
atunci cantitatea conservată este ``H(t)=c \text{ln } S(t) - dS(t) - a T(t) + b\text{ln} T(t).`` ``H`` are un punct de maxim pentru ``(S,T)=\left(\frac{c}{d},\frac{b}{a}\right)``. Avem graficul lui ``H`` și alături simularea sistemului.
"""

# ╔═╡ 79367e8b-7051-4752-9b30-bb32b92783ea
md"""

Un alt exemplu este dat de ecuația pendulului. Sistemul este
```math
\begin{aligned}
\theta' & = & \omega \\
\omega' & = & -\sin \theta,
\end{aligned}
```
unde ``\theta`` și ``\omega``sînt unghiurile făcute de firul pendului și de normala din centrul de greutate al greutății, ambele cu verticala. Punctele de echilibruvor fi de forma ``(0,k\pi)`` cu ``k \in \mathbb{Z}.`` La fel ca în cazul precedent valorile proprii ale Jacobienei pentru ``(0,0)`` vor fi ``\pm i.`` Cantitatea conservată este energia totală
```math
E = \frac{1}{2}\omega^2-\cos\theta,
```
deci acesta este un centru. Trecînd la ``(0,\pi)`` valorile proprii vor fi ``\pm 1``, deci acesta va fi un punct șa. Simularea este:
"""

# ╔═╡ af271326-a8e9-42e1-9e7e-04ff41b186d6
begin
	e(x,y) = 1/2*x^2-cos(y)
	g(x) = Point2f(x[2],-sin(x[1]))
	xs = LinRange(-2.0, 2.0, 500)
    ys= LinRange(-2.0, 2.0, 500)
	zs = [e(x,y) for x ∈ xs, y ∈ ys]
	fig2=Figure(size=(1000,500))
	ax21 = Axis3(fig2[1,1],elevation=0)
	surface!(ax21,xs,ys,zs)
	scatter!(ax21,[0],[0],[e(0,0)],markersize = 20.0,color=:red)
	ax22 = Axis(fig2[1,2])
	streamplot!(ax22,g,5*xs,5*ys,stepsize = 0.1)
	scatter!(ax22,[-2π,-π,0,π,2π],[0,0,0,0,0],markersize = 10, color=:red)
	fig2
end

# ╔═╡ 63f160e5-b075-431c-9648-ff7caf3a25e7
md"""
Se observă că avem orbite care pornesc și se termină în punctele șa. Aceste orbite se numesc **homoclinice**.Ele nu sînt periodice ci se întorc de unde au pleca, dar într-un timp infinit.

## Bifurcația Poincaré-Andronov-Hopf 

Bifurcația Hopf reprezintă una din căile principale către comportament periodic. Să luăm ca exemplu clarinetul lui Raileigh cu un parametru:
```math
\begin{aligned}
x' & = & v \\
v' & = & -x - c(v^3-v)
\end{aligned}
```
Singurul punct de echilibru este ``(0,0)``. Matricea Jacobiană este
```math
J(x,v)=\left(\begin{array}{ccc} 0 & 1 \\-1 & -c(3v^2-1) \end{array}\right).
```
Evaluată în punctul de echilibru avem
```math
J(0,0)=\left(\begin{array}{ccc} 0 & 1 \\-1 & -c \end{array}\right)
```
cu ecuația caracteristică ``\lambda^2-c\lambda+1=0.`` Rădăcinile acesteia sînt 
```math
\frac{c \pm \sqrt{c^2-4}}{2}.
```
Studiem ce se întîmplă cînd variem ``c`` în jurul lui ``0``. Expresia de sub radical va fi negativă deci valorile proprii sînt complexe. Pentru ``c<0`` partea reală este negativă, deci echilibru este o spirală stabilă. Pentru ``c=0`` avem un centru, apoi echilibrul devine instabil, dar apare un ciclu limită stabil. 
"""

# ╔═╡ 0a520d03-8843-4a3d-9e34-3105f8177101
@bind c PlutoUI.Slider(-1.0:0.1:1.0, show_value=true)

# ╔═╡ c5159ecb-2a91-4552-9f04-93206866f6af
begin
r(x)=Point2f(x[2],-x[1]-c*(x[2]^3-x[2]))
	fig3=Figure(size=(1000,500))
	ax3=Axis(fig3[1,1])
	streamplot!(ax3,r,xs,ys,stepsize = 0.01,gridsize=(60,60))
	scatter!([0],[0],markersize=10,color=:green)
	fig3
end

# ╔═╡ a2a48931-17cb-46c7-9ad3-79639e5758a0
md"""
# Alte bifurcații ale punctelor de echilibru pentru sisteme continue

Ne reamintim că pentru ecuații ``x'=f(x,c)`` cu ``f`` funcție reală, orbitele fie tind către un echilibru stabil, fie tind la infinit. În acest caz avem două tipuri de bifurcații:
  - **șa-nod** în care două echilibre unul stabil și unul instabil se apropie se ciocnesc și apoi dispar;
  - **transcritică** două puncte critice se apropie, se ciocnesc, apoi schimbă stabilitatea;
  - **trident** dintr-un punct fix apar încă două cu stabilitatea opusă față de cel inițial.

Aceleași bifurcații apar și în dimensiuni mai mari. Pentru bifurcația șa-nod exemplul tipic este
```math
\begin{aligned}
x'& = & \mu & - & x^2 \\
y' & = & -y
\end{aligned}
```
Se vede că avem două puncte fixe ``(-\sqrt{\mu},0)`` și ``(\sqrt{\mu},0)``. Primu este punct șa, iar cel de-al doilea este un nod stabil. Cînd ``\mu`` scade atunci punctele se apropie, se unesc pentru ``\mu=0``, apoi dispar.

Exemplul pentru bifurcația transcritică este
```math
\begin{aligned}
x'& = & \mu x& - & x^2 \\
y' & = & -y,
\end{aligned}
```
iar pentru trident este
```math
\begin{aligned}
x'& = & \mu x & - & x^3 \\
y' & = & -y,
\end{aligned}
```

Bifurcația Hopf este tipică pentru dimensiunile superioare și nu are un corespondent în dimensiune 1. 
"""

# ╔═╡ 64b18bec-91c2-40af-8309-2423a9365963
@bind μ PlutoUI.Slider(0.5:-0.01:-0.5,show_value=true)

# ╔═╡ 72200ad4-4b5d-4e3d-8f8d-950db1fe0888
begin
s(x)=Point2f(μ-x[1]^2,-x[2])
t1(x)=Point2f(μ*x[1]-x[1]^2,-x[2])
	t2(x)=Point2f(μ*x[1]-x[1]^3,-x[2])
	fig4=Figure(size=(2000,500))
	ax41=Axis(fig4[1,1],title="Bifurcația șa-nod")
	ax42=Axis(fig4[1,2],title="Bifurcația transcritică")
	ax43=Axis(fig4[1,3],title="Bifurcația trident")
	streamplot!(ax41,s,xs,ys,stepsize = 0.01,gridsize=(60,60))
	streamplot!(ax42,t1,xs,ys,stepsize = 0.01,gridsize=(60,60))
	streamplot!(ax43,t2,xs,ys,stepsize = 0.01,gridsize=(60,60))
	# scatter!([sqrt],[0],markersize=10,color=:green)
	fig4
end

# ╔═╡ 55df461d-095b-4f16-890f-bd38ff22c78b
md"""
# Bifurcații ale punctelor de echilibru pentru sisteme discrete

Dacă avem sistemul discret ``x_{t+1}=f(x_t,\alpha)`` definit pe ``\mathbb{R}^n``. Bifurcațiile apar într-unul din următoarele cazuri:
  - apare o valoare proprie ``\lambda_1 = 1``: bifurcația **fold** (sau tangentă);
  - apare o valoare proprie ``\lambda_1 = -1``: bifurcația **flop** (sau dublarea perioadei);
  - apare o valoare proprie ``\lambda_1`` complexă de modul 1: bifurcația **Neimark-Sacker** (sau torică);

Aceste bifurcații sînt destul de asemănătoare cu cazul continuu:
  - în cazul bifurcației tangente apar două puncte de echilibru unul stabil și unul instabil;
  - În cazul bifurcației de dublare a perioadei un punct fix stabil se transformă într-o orbită periodică de lungime ``2``;
  - În cazul bifurcației Neimark-Sacker o spirală stabilă dă naștere unui cerc invariant, iar punctul fix devine instabil. 
"""

# ╔═╡ da7805bc-08cb-4c77-a8f9-0917983161be
@bind r0 PlutoUI.Slider(0.3:0.01:0.5,show_value=true)

# ╔═╡ 60467996-85ea-4341-ae20-310ef2707d9d
begin
	function henon_rule(u,p,n)
		x, y = u # system state
        r,b = p # system parameters
        xn = r - x^2 + 0.3*y
        yn = b*x
        return SVector(xn, yn)
	end
	p0 = [r0,1.0]
	init_vals = [[0.2,0.3],[0.5,0.3],[0.1,0.2],[0.3,0.1]]
	u0 = [0.2, 0.3]
	prbs= []
	total_time = 1000
	for u in init_vals
		henon = DeterministicIteratedMap(henon_rule, u, p0)
    	X, t = trajectory(henon, total_time)
		push!(prbs,(X,t))
	end
	fig5 = Figure(size=(1000,700))
	Label(fig5[0,:],text="Dublarea Perioadei",fontsize=40,justification = :center,tellwidth = false)
		ax51=Axis(fig5[1,1])
		ax52=Axis(fig5[1,2])
		for (X,t) in prbs
			scatter!(ax51,X,markersize=10)
			for var in columns(X)
				scatter!(ax52,t,var)
			end
		end
	fig5
end

# ╔═╡ 15adc7a7-8275-45de-a27e-4b998feb51e4
@bind r1 PlutoUI.Slider(0.7:0.01:1.9,show_value=true)

# ╔═╡ 67f69caf-0a29-4b4e-a4da-4b5d2af0d3e1
let
begin
	function henon_rule1(u,p,n)
		x, y = u # system state
        a,b = p # system parameters
        xn = a*x - b*y - x*(x^2+y^2)
        yn = b*x + a*y - y*(x^2+y^2)
        return SVector(xn, yn)
	end
	p0 = [r1,0.5]
	init_vals = [[0.2,0.3],[0.5,0.3],[0.1,0.2],[0.3,0.1]]
	u0 = [0.2, 0.3]
	prbs= []
	total_time = 1000
	for u in init_vals
		henon = DeterministicIteratedMap(henon_rule1, u, p0)
    	X, t = trajectory(henon, total_time)
		push!(prbs,(X,t))
	end
	fig5 = Figure(size=(1000,700))
	Label(fig5[0,:],text="Neimark-Sacker",fontsize=40,justification = :center,tellwidth = false)
		ax51=Axis(fig5[1,1])
		ax52=Axis(fig5[1,2])
		for (X,t) in prbs
			scatter!(ax51,X,markersize=10)
			for var in columns(X)
				scatter!(ax52,t,var)
			end
		end
	fig5
end
end

# ╔═╡ f9c5113c-df24-4b02-ab3f-9bd53c175995
@bind r2 PlutoUI.Slider(0.01:-0.0001:-0.2,show_value=true)

# ╔═╡ 0212bc42-e41d-476e-aba4-c87d61cfd3ad
let
begin
	function henon_rule(u,p,n)
		x, y = u # system state
        r,b = p # system parameters
        xn = r - x^2 + 0.3*y
        yn = b*x
        return SVector(xn, yn)
	end
	p0 = [r2,1.0]
	init_vals = [[0.2,0.3],[0.5,0.3],[0.1,0.2],[0.3,0.1]]
	u0 = [0.2, 0.3]
	prbs= []
	total_time = 100
	for u in init_vals
		henon = DeterministicIteratedMap(henon_rule, u, p0)
    	X, t = trajectory(henon, total_time)
		push!(prbs,(X,t))
	end
	fig5 = Figure(size=(1000,700))
	Label(fig5[0,:],text="Bifurcația Fold",fontsize=40,justification = :center,tellwidth = false)
		ax51=Axis(fig5[1,1])
		ax52=Axis(fig5[1,2])
		for (X,t) in prbs
			scatter!(ax51,X,markersize=10)
			for var in columns(X)
				scatter!(ax52,t,var)
			end
		end
	fig5
end
end

# ╔═╡ Cell order:
# ╠═e5b1174c-d529-11f0-2a30-6f5ef7b91715
# ╠═11d49c9f-f3d6-4d1b-bfd4-08dd1491cc97
# ╟─7f2464da-120e-4774-b364-5ed61eb97272
# ╟─522e922b-287c-4c34-87cc-b29bcca0cd40
# ╟─79367e8b-7051-4752-9b30-bb32b92783ea
# ╠═af271326-a8e9-42e1-9e7e-04ff41b186d6
# ╟─63f160e5-b075-431c-9648-ff7caf3a25e7
# ╟─0a520d03-8843-4a3d-9e34-3105f8177101
# ╠═c5159ecb-2a91-4552-9f04-93206866f6af
# ╟─a2a48931-17cb-46c7-9ad3-79639e5758a0
# ╠═64b18bec-91c2-40af-8309-2423a9365963
# ╠═72200ad4-4b5d-4e3d-8f8d-950db1fe0888
# ╟─55df461d-095b-4f16-890f-bd38ff22c78b
# ╠═da7805bc-08cb-4c77-a8f9-0917983161be
# ╟─60467996-85ea-4341-ae20-310ef2707d9d
# ╠═15adc7a7-8275-45de-a27e-4b998feb51e4
# ╟─67f69caf-0a29-4b4e-a4da-4b5d2af0d3e1
# ╟─f9c5113c-df24-4b02-ab3f-9bd53c175995
# ╟─0212bc42-e41d-476e-aba4-c87d61cfd3ad
