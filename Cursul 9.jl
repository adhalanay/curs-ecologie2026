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

# ╔═╡ e84d2096-cfa6-11f0-2d5f-3fa172654109
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

# ╔═╡ e21e98bb-51a2-4d7c-8844-75c01a79fac5
WidthOverDocs()

# ╔═╡ 8b2a6ae7-cc9d-4207-a7c3-7512530056f9
md"""
# Aplicații de mai multe variabile

Pentru a studia stabilitatea sistemelor neliniare avem nevoie să studiem funcțiile arbitrare de forma ``f:\mathbb{R}^n \to \mathbb{R}^n.`` În cazul continu un astfel de sistem neliniar este de forma
```math
\begin{aligned}
x'_1 & = & f_1(x_1,x_2,\dots,x_n) \\
x'_2 & = & f_2(x_1,x_2,\dots,x_n) \\
\vdots \\
x'_n & = & f_n(x_1,x_2,\dots,x_n).
\end{aligned}
```
Punctele sale de echilibru $(x_1^*,\dots,x_n^*)$ se determină din condiția $x_1'=0,\dots,x_n'=0$. În acest curs vom studia stabilitatea echilibrelor. Ne amintim  că stabilitatea punctelor de echilibru pentru ecuațiile 
diferențiale scalare $x'=f(x)$ este determinată de panta tangentei la graficul funcției $f$, adică de aproximarea liniară a acestuia în punctul de echilibru. Ne propunem să dezvoltăm o teorie similară pentru funcții de mai multe variabile. 

Presupunem că avem acum o funcție de două variabile și cu valori reale $f:\mathbb{R}^2 \to \mathbb{R}$. Graficul acestei funcții 
definește o suprafață. De exemplu dacă $z=f(x,y)=5-x^2-y^2$, graficul este

$(PlutoUI.LocalResource("cuadr.png"))

"""

# ╔═╡ 2a6e5eed-19fd-455c-a35a-a9121dc47ff2
md"""

# Liniarizarea sistemelor neliniare

Vrem să vedem care este aproximarea liniară a unui astfel de grafic într-un punct $(x_0,y_0)$. Dacă funcția $f$ este liniară, atunci imaginea sa este un plan. Deci aproximarea liniară a graficului va trebui să fie tot un plan. Acest plan va fi numit planul tangent la graficul lui $f$. În general un plan în spațiu este determinat de două direcții, adică doi vectori liniar-independeți. Pentru a obține aceste direcții vom tăia graficul cu anumite plane convenabile. Pentru început putem lua planul 
$YZ$ dat de $X=1$ și obținem curba $z=1-y^2$ din planul $YZ$. Funcția $z \mapsto 1-y^2$ este o funcție de o variabilă și deci putem calcula
```math
\left.\frac{dz}{dy}\right|_{y=y_0}=-2y_0.
```
Aproximarea liniară a acestei curbe în punctul $(1,y_0,1-y_0^2)$ este $dz=-2y_0 dy$. Derivata pe care am obținut-o poate fi obținută și derivînd funcția $f$, după ce am presupus $x$ constant. Notația este
```math
\left.\frac{\partial f}{\partial y}\right|_{y=y_0}.
```
Deci aproximarea liniară de-a lungul curbei $x=const$ este $\left.\frac{\partial f}{\partial y}\right|_{y=y_0}.$ La fel dacă 
luăm curba $y=const$, atunci aproximarea în lungul ei este 
```math
\left.\frac{\partial f}{\partial x}\right|_{x=x_0}=-2x_0.
```
Am găsit astfel a doua dreaptă care definește planul tangent este $\left.\frac{\partial f}{\partial x}\right|_{x=x_0}.$ Prin 
analogie cu cazul scalar, obținem că ecuația planului tangent este
```math
z-z_0=\left.\frac{\partial f}{\partial x}\right|_{(x_0,y_0)}(x-x_0)+\left.\frac{\partial f}{\partial y}\right|_{(x_0,y_0)}(y-y_0).
```

Trecem acum la o funcție $\mathbb{R}^2 \to \mathbb{R}^2$, $V(x,y)=\left(f(x,y),g(x,y)\right).$ Presupunem că $g(x,y)=5+x^2-y^2.$
Planul tangent într-un punct $(x_0,y_0)$ este de asemenea
```math
z-z_0=\left.\frac{\partial g}{\partial x}\right|_{(x_0,y_0)}(x-x_0)+\left.\frac{\partial g}{\partial y}\right|_{(x_0,y_0)}(y-y_0).
```

Lipind cele două componente obținem că aproximarea liniară a funcției. Putem scrie că aceasta este
```math
\left(\Delta x,\Delta y\right)\mapsto\left(\Delta z, \Delta w\right)=\left.\left(\frac{\partial f}{\partial x}
\Delta x+\frac{\partial f}{\partial y}\Delta y,\frac{\partial g}{\partial x}\Delta x+\frac{\partial g}{\partial y}\Delta y\right)\right|_{(x_0,y_0)}.
```
Prin urmare aplicația liniară care aproximează funcția $V$ are matricea
```math
\left.\left(\begin{array}{cc}\frac{\partial f}{\partial x} & \frac{\partial f}{\partial y} \\ 
\frac{\partial g}{\partial x} & \frac{\partial g}{\partial y}\end{array}\right)\right|_{(x_0,y_0),}
```
numită **matricea Jacobiană** a lui $V$.

Putem spune că matricea Jacobiană este aproximarea liniară a funcției într-un punct. Ne interesează acum să studiem comportamentul unui sistem de ecuații neliniare în jurul unui punct de echilibru. Acest lucru nu poate fi făcut întotdeauna. Condițiile cînd aproximarea este corectă sînt date de 

**Teorema Hartmann-Grobman.**  Fie sistemul de ecuații diferențiale ``x'=f(x)``, unde ``x=\left(x_1,\dots,x_n\right)`` și ``f:\mathbb{R}^b \to \mathbb{R}^n``. Fie ``x^*`` un punct de echilibru. Dacă nici o valoare proprie a matricei Jacobiene nu are partea reală ``0`` atunci comportamentul local în jurul lui ``x^*`` este același cu cel al sistemuliui liniarizat.

Mai concret, pentru o funcție $V:\mathbb{R}^n \to \mathbb{R}^n$, matricea jacobiană
```math
D_{V,(x_1,\dots,x_n)_0}=\left.\left(\begin{array}{cccc}\frac{\partial f_1}{\partial x_1} & \dots &\frac{\partial f_1}{
\partial x_n} \\
\vdots & \ddots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \dots & \frac{\partial f_n}{\partial x_n}\end{array}\right)\right|_{(x_1,\dots,x_n)_0.}
```
Matricea definește un cîmp de vectori $(x_1,\dots,x_n) \mapsto D_{V,(x_1,\dots,x_n)},$ notată pe scurt cu $D.$ Valorile 
proprii sînt rădăcinile polinomului
```math
\det(D-\lambda I_n).
```
Valorile proprii descompun spațiul în subspații în care $D$  acționează astfel:
 - echilibru stabil $\lambda < 0;$
 - echilibru instabil $\lambda > 0;$
 - spirală stabilă $\lambda=a+ib$ cu $a < 0;$
 - spirală instabilă $\lambda=a+ib$ cu $a > 0.$
 
Rezumînd, procedăm în felul următor:
 1. Calculăm liniarizarea sistemului $V;$
 2. Folosind metoda valorilor proprii determinăm stabilitate sistemului liniarizat;
 3. Dacă nici o valoare proprie nu este $0$ sau nu are partea reală $0$, atunci comportamentul sistemului neliniar este 
    echivalent cu cel al liniarizării. 

## Exemple
Ecuația unui oscilator cu frecare introdusă de Rayleigh este
```math
\begin{aligned}
x' & = & v \\
v' & = & -x -(v^3-v).
\end{aligned}
```
Ea are un singur punct de echilibru $(0,0)$. Pentru a studia stabilitatea studiem matricea jacobiană în acesta. Matricea este
```math
J=\left(\begin{array}{cc} 0 & -1 \\ 1 & 1\end{array}\right).
```
Polinomul caracteristic este $\lambda^2-\lambda+1$, cu rădăcinile $\frac{1\pm\sqrt{3}i}{2}.$ Deci traiectoriile sînt spirale 
instabile.
"""

# ╔═╡ 66dae57c-9218-4da5-a108-682adf36bff1
begin
	f1(x) = Point2f(x[2],-x[1]-(x[2]^3-x[2]))
	fig1 = Figure()
	ax1 = Axis(fig1[1,1],title="x\'= v ,v\'=-x-(v^3-v) ")
	xs1 = LinRange(-2, 2, 10)
    ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax1,f1,xs1,ys1,colormap = :magma,linewidth = 1)
	fig1
end

# ╔═╡ 76cb5516-ce62-4b1f-86c4-a51de9f0ca0a
md"""
### Sistemul căprioare-elani
Am văzut mai demult modelul de competiție pentru resurse între căprioare și elani
```math
\begin{aligned}
D' & = & 3D & - & 2MD & - & D^2 \\
M' & = & 2M & - & DM & - & M^2.
\end{aligned}
```
Am văzut că echilibrul netrivial este $(1,1)$. Matricea jacobiană este
```math
J=\left(\begin{array}{cc} -1 & -1 \\ -2 & -1\end{array}\right).
```
Polinomul caracteristic este $\lambda^2+2\lambda-1$, cu rădăcinile $-1\pm\sqrt{2}.$ Una este pozitivă și una negativă, deci
punctul de echilibru este un punct șa, care este un echilibru instabil, deci speciile nu pot coexista.
"""

# ╔═╡ 989f46f6-ce35-4c24-93c1-b594334046e5
begin
	f2(x) = Point2f(3*x[1]-2*x[1]*x[2]-x[1]^2,2*x[2]-x[1]*x[2]-x[2]^2)
	fig2 = Figure()
	ax2 = Axis(fig2[1,1],title="Căprioare-Elani")
	xs2 = LinRange(0, 3, 10)
    ys2 = LinRange(0, 3, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	scatter!(ax2,1,1,markersize=20,color=:red)
	streamplot!(ax2,f2,xs2,ys2,colormap = :magma,linewidth = 1)
	fig2
end

# ╔═╡ d9f80311-d217-4a5d-b54e-a35ebc22fd7f
md"""
# Cînd nu funcționează Hartmann-Grobman-Cazul valorilor proprii nule

Să vedem ce  se întîmplă cînd avem un sistem în care una dintre valorile proprii este ``0``. Considerăm sistemul
```math
\begin{aligned}
x' & = & x & - & 2y \\
y' & = & 3x & - & 6y.
\end{aligned}
```
Matricea sistemului este 
```math
M=\left(\begin{array}{cc} 1 & -2 \\ 3 & -6\end{array}\right).
```
Valorile proprii sînt ``0`` și ``-5.`` Cu direcțiile proprii ``y = 0.5 x`` respectiv ``y=3x.`` Din simulare 
"""

# ╔═╡ 96f62795-1e06-4125-b6f4-0a645f7af7c7
begin
	f3(x) = Point2f(x[1]-2*x[2],3*x[1]-6*x[2])
	fig3 = Figure()
	ax3 = Axis(fig3[1,1])
	# xs2 = LinRange(0, 3, 10)
    # ys2 = LinRange(0, 3, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	scatter!(ax3,0,0,markersize=20,color=:red)
	streamplot!(ax3,f3,xs1,ys1,colormap = :magma,linewidth = 1)
	lines!(ax3,[-2,2],[-1,1],linewidth=3)
	fig3
end

# ╔═╡ 11617e30-e047-4a2b-a973-dee84b22ee55
md"""
se vede că avem nu un punct fix ci o dreaptă întreagă de puncte fixe. Fiecare valoare inițială va fi atrasă de de un punct de pe dreapta albastră. Problema cu aceste sisteme este că ele nu sînt robuste. Dacă perturbăm puțin sistemul obținem situația de mai jos:
"""

# ╔═╡ a3463893-f65b-4882-8f53-38808de8b92c
@bind ϵ PlutoUI.Slider(-0.6:0.01:0.6,show_value=true)

# ╔═╡ bff4de5d-f28a-4d7e-9d55-7aed912edd3f
begin
	f4(x) = Point2f(x[1]-2*x[2],(3+ϵ)*x[1]-6*x[2])
	fig4 = Figure()
	ax4 = Axis(fig4[1,1])
	# xs2 = LinRange(0, 3, 10)
    # ys2 = LinRange(0, 3, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	# scatter!(ax4,0,0,markersize=20,color=:red)
	streamplot!(ax4,f4,xs1,ys1,colormap = :magma,linewidth = 3)
	# lines!(ax3,[-2,2],[-1,1],linewidth=3)
	fig4
end

# ╔═╡ 73c5aa5f-ff3a-4293-83b4-4ea2f3843015
md"""
Se vede că pe măsură ce ``\epsilon`` crește și devine din negativ pozitiv, originea (care este singurul punct fix) trece dintr-un punct șa într-o spirală stabilă. Sistemele care sînt robuste și nu își schimbă comportamentul la schimbări ale parametrilor se numesc **structural stabile.** Toate sistemele de ecuații diferențiale folosite pentru modelarea unor fenomene naturale trebuie să fie structural stabile.

# Cînd nu funcționează Hartmann-Grobman-Cazul valorilor proprii pur imaginare

Să luăm ``f(x,y)=\frac{x^2+y^2}{1+x^2+y^2}`` și sistemul 
```math
\begin{aligned}
x' & = & -y & + & f(x,y)x \\
y' & = & x & + & f(x,y)y
\end{aligned}
```
Simularea este 
"""

# ╔═╡ 560ab720-b3da-4ac0-abc8-972cf5ed0ab7
begin
	f(x,y) = (x^2+y^2)/(1+x^2+y^2)
	f5(x) = Point2f(-x[2]+f(x[1],x[2])*x[1],x[1]+f(x[1],x[2])*x[2])
	fig5 = Figure()
	ax5 = Axis(fig5[1,1])
	# xs2 = LinRange(0, 3, 10)
    # ys2 = LinRange(0, 3, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax5,f5,xs1,ys1,colormap = :magma,linewidth = 1)
	fig5
end

# ╔═╡ 05dfd111-0a31-4f7c-9f62-0ccd86916202
md"""
arătînd că originea este o spirală instabilă. Pe de altă parte Jacobianul are valorile proprii ``\pm i`` și deci prezice că originea ar fi un centru.

# Sisteme discrete neliniare

``\newcommand{\vecx}{\overrightarrow{x}}``

Pentru sistemul discret ``\vecx_{n+1}=f(\vecx_n)`` stabilitatea unui punct fix ``\vecx^*`` se poate determina cu ajutorul matricei jacobiene a lui ``f`` evaluată în ``\vecx^*``. Mai precis avem:
 - dacă toate valorile proprii au modulul mai mic decît ``1`` atunci echilibrul este stabil;
 - dacă cel puțin una dintre valorile proprii are modulul mai mare decît ``1`` atunci echilibrul este instabil.

Este mai ușor uneori în loc să calculăm valorile proprii să folosim condțiile lui Jury. În cazul dimensiunii 2 polinomul caracteristic este 
```math
\lambda^2 - \text{tr} J \lambda +\det J = 0.
```
Atunci echilbrul este stabil dacă și numai dacă 
```math
|\text{tr} J | < 1 + \det J <2.
```
## Exemplu: un model de interacțiune gazdă-parazitoid

Parazitoizi sînt insecte care depun ouăle pe sau în corpul gazdelor. Ouăle dau la naștere la larve care ucid gazda. Gazdele parazitate dau naștere la următoarea generație de parazitoizi, iar cele neparazitate la următoarea generație de gazde. Fie ``H_n`` și
``P_n`` numărul de gazde, respectiv parazitoizi. Fie ``f\left(H_n,P_n\right)`` proporția de gazde ne parazitate. Modelul va fi
```math
\begin{aligned}
H_{n+1} & = & kf(H_n,P_n) H_n \\
P_{n+1} & = & c(1-f(H_n,P_n)) P_n.
\end{aligned}
```
Funcția ``f(H_n,P_n) = e^{-aP_n}, k > 1`` este rata de reproducere al gazdelor, ``c`` este numărul mediu de ouă depuse de parazitoizi în gazde.

Modelul are două puncte fixe ``(H_1^*,P_1^*)=(0,0)`` și 
```math
(H_2^*,P_2^*) = \left(\frac{k \text{ln}k}{ac(k-1)},\frac{\text{ln}k}{a}\right).
```

Punctul ```(0,0)``` reprezintă cazul în care cele două populații sînt dispărute, așa că echilibrul interesant este cel netrivial. Pentru acesta matricea jacobiană este 
```math
\left(\begin{array}{cc} 1 & -\frac{k \text{ln}k}{c(k-1)} \\
\frac{c(k-1)}{k} & \frac{\text{ln} k }{k-1} 
\end{array}\right).
```
Avem 
```math
\begin{aligned}
\text{tr} J & = & 1 + \frac{\text{ln} k}{k-1} \\
\det J & = & \text{ln } k + \frac{\text{ln} k}{k-1}.
\end{aligned}
```
Cum ``k > 1`` prima inegalitate a lui Jury este satisfăcută, dar ``\det J > 1,`` deci a doua inegalitate nu este satisfăcută nici o dată. Deci echilibrul netrivial nu este stabil. 
"""

# ╔═╡ 61275058-934b-47ed-9448-c7d4b88336c6
begin
	a = 0.005
	k = 1.05
	c = 3.0
	fn(x,y) = exp(-a*y)
	function nicholson(u,p,n)
		x, y = u # system state
        k,c = p # system parameters
        xn = k*fn(x,y)*x
        yn = c*(1-fn(x,y))*x
        return SVector(xn, yn)
	end
	p0 = [k,c]
	init_vals = 10*[[5.0,3.0]] #,[8.0,2.0],[5.0,1.0],[3.0,2.0]]
	u0 = [50.0, 10.0]
	prbs= []
	total_time = 200
	for u in init_vals
		nichols = DeterministicIteratedMap(nicholson, u, p0)
    	X, t = trajectory(nichols, total_time)
		push!(prbs,(X,t))
	end
	fig6 = Figure(size=(1000,500))
		ax6=Axis(fig6[1,1],title="Portretul de fază")
		ax61=Axis(fig6[1,2],title="Seriile temporale")
		for (X,t) in prbs
			lines!(ax6,X,linewidth=5)
			for var in columns(X)
				scatter!(ax61,t,var)
			end
		end
	fig6
end

# ╔═╡ Cell order:
# ╟─e84d2096-cfa6-11f0-2d5f-3fa172654109
# ╟─e21e98bb-51a2-4d7c-8844-75c01a79fac5
# ╟─8b2a6ae7-cc9d-4207-a7c3-7512530056f9
# ╟─2a6e5eed-19fd-455c-a35a-a9121dc47ff2
# ╟─66dae57c-9218-4da5-a108-682adf36bff1
# ╟─76cb5516-ce62-4b1f-86c4-a51de9f0ca0a
# ╟─989f46f6-ce35-4c24-93c1-b594334046e5
# ╟─d9f80311-d217-4a5d-b54e-a35ebc22fd7f
# ╟─96f62795-1e06-4125-b6f4-0a645f7af7c7
# ╟─11617e30-e047-4a2b-a973-dee84b22ee55
# ╠═a3463893-f65b-4882-8f53-38808de8b92c
# ╟─bff4de5d-f28a-4d7e-9d55-7aed912edd3f
# ╟─73c5aa5f-ff3a-4293-83b4-4ea2f3843015
# ╟─560ab720-b3da-4ac0-abc8-972cf5ed0ab7
# ╟─05dfd111-0a31-4f7c-9f62-0ccd86916202
# ╟─61275058-934b-47ed-9448-c7d4b88336c6
