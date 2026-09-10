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

# ╔═╡ e2576524-c9f3-11f0-180f-1b51b8a5b421
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

# ╔═╡ 2a08f8d0-a81c-40f1-86c0-9984f1a27cea
begin
	using LinearAlgebra
	let
		M=[0.65 0.5;0.25 0.9]
println(eigvals(M))
println(eigvecs(M))
init1 = [10,5]
init2 = [10,-5]
init3 = [5,20]
init4 = [10,30]
function orbits(init,iters)
    xs = []
    ys = []
    val = init
    push!(xs,val[1])
    push!(ys,val[2])
    for i in 1:iters
        val = M*val
        push!(xs,val[1])
        push!(ys,val[2])
    end
    xs,ys
end
n = 15
xs1,ys1 = orbits(init1,n)
xs2,ys2 = orbits(init2,n)
xs3,ys3 = orbits(init3,n)
xs4,ys4 = orbits(init4,n)
fig3 = Figure(size=(500,500))
ax1 = Axis(fig3[1,1])
scatter!(ax1,xs1,ys1)
scatter!(ax1,xs2,ys2)
scatter!(ax1,xs3,ys3)
scatter!(ax1,xs4,ys4)
fig3
	end
end

# ╔═╡ 71b64ca1-62e8-431b-bdf3-58d12810d900
WidthOverDocs()

# ╔═╡ d92c50aa-19c6-4486-9612-4a04ed932d85
md""" a = $(@bind a NumberField(-2.0:0.1:1.5))"""

# ╔═╡ fced5079-d314-45d4-a956-f765df74b99c
md""" b = $(@bind b NumberField(-2.0:0.1:1.5))"""

# ╔═╡ d3721a4f-9d81-4220-86dd-82d70b181b64
md""" c = $(@bind c NumberField(-2.0:0.1:1.5))"""

# ╔═╡ aac19b27-6864-444f-81fd-f2317eaab47b
md""" d = $(@bind d NumberField(-2.0:0.1:1.5))"""

# ╔═╡ d00845bc-f88b-4afa-a7ab-fec17b35c033
md""" numărul de pași: $(@bind timp PlutoUI.Slider(1:2000,show_value=true))"""

# ╔═╡ 6391e7e3-c7a7-43b4-b6b7-4a5ca4d6b398
begin
	function gen(u,p,n)
		x,y = u
		a,b,c,d = p
		xn = a*x+b*y
		yn = c*x+d*y
		return SVector(xn,yn)
	end
	u0 = [30.0,80.0]
	p0 = [a,b,c,d]
	sislin = DeterministicIteratedMap(gen, u0, p0)
	X,t = trajectory(sislin,timp)
	fig1 = Figure(size=(1000,500))
	ax1 = Axis(fig1[1,1],title="Sisteme liniare cu parametrii $a, $b, $c, $d")
	scatter!(ax1,u0[1],u0[2],markersize=15,color=:red)
	lines!(ax1,X)
	ax10 = Axis(fig1[1,2], title="Seriile temporale")
	for v in columns(X)
		scatter!(ax10,v)
	end
	fig1
end

# ╔═╡ 252de602-5ce3-48e7-8321-3da20f88a59a
begin
	let
	function gen(u,p,n)
		x,y = u
		a,b,c,d = p
		xn = a*x+b*y
		yn = c*x+d*y
		return SVector(xn,yn)
	end
	u0 = [30.0,80.0]
	p0 = [cos(π/25),sin(π/25),-sin(π/25),cos(π/25)]
	sislin = DeterministicIteratedMap(gen, u0, p0)
	X,t = trajectory(sislin,timp)
	fig2 = Figure(size=(1000,500))
	ax2 = Axis(fig2[1,1],title="Sisteme liniare cu parametrii $(p0[1]),$(p0[2]),$(p0[3]),$(p0[4])")
	scatter!(ax2,u0[1],u0[2],markersize=15,color=:red)
	lines!(ax2,X)
	ax20 = Axis(fig2[1,2], title="Seriile temporale")
	for v in columns(X)
		scatter!(ax20,v)
	end
	fig2
	end
end

# ╔═╡ 5c4a5a48-8cd4-4444-865d-1b8c5517a4fd
md"""
# Sisteme liniare discrete

``\newcommand{\vecv}{\overrightarrow{v}}``
Să presupunem pentru început că avem un sistem dat de o matrice diagonală. Adică un sistem 
```math
\begin{aligned}
x_{n+1} & = & a \cdot x_{n} \\
y_{n+1} & = & d \cdot y_{n}
\end{aligned}
```
cu valorile inițiale  pozitive $(x_0,y_0)$. Dacă facem cîteva iterații vedem că $x_1=a\cdot x_0$, $x_2=ax_1=a^2x_0$, $x_3=a^3x_0$ etc. Deci în general
$x_n=a^nx_0$ și la fel $y_n=d^ny_0$. La nivel de matrice vom avea că 
```math
A^n = \left(\begin{array}{cc} a^n & 0 \\ 0 & d^n\end{array}\right)
```

Studiem comportamentul pe termen lung al lui $x_n$ (cel al lui $y_n$ este identic). 
 - Dacă $|a| > 1$ și $x_0 > 0$, atunci 

```math
      \lim_{n \to \infty} x_n=x_0 \lim_{n\to\infty} a^n = \infty;
```

 - Dacă $|a|=1,$ atunci $x_n=\pm x_0$ și avem un punct fix sau unul periodic cu perioadă $2$;
 - Dacă $|a| < 1$ atunci 

```math
\lim_{n \to \infty} x_n =0.
```

La fel se întîmplă pentru $y_n$. Comportamentul va fi dat de compunerea celor două componente:
 - dacă ambii parametri au modulul mai mare decît $1$, atunci originea este un echilibru instabil;
 - dacă ambii parametri au modulul mai mic decăt $1$, atunci originea este un echilibru stabil;
 - dacă $|a| \leq 1$ și $|d|>1$ sau invers, atunci comportamentul pe termen lung va fi determinat de coficientul dominant (cel cu modulul mai mare), iar originea va fi un punct șa.
 Semnul fiecărui coeficient ne spune dacă sistemul va oscila în jurul unei axe. De exemplu dacă $a<0$, atunci sistemul va oscila
 în jurul axei $Oy$, adică $x_n$ va lua alternativ valori pozitive și negative. 

## Valori proprii reale

Studiem acum cazul unei matrice generale. Vom presupune pentru început că matricea are două valori proprii reale $\lambda_1$
și $\lambda_2$. Fie matricea
```math
A=\left(\begin{array}{cc} a & b \\ c & d\end{array}\right).
```
Atunci sistemul discret este 
```math
\begin{aligned}
x_{n+1} &=& ax_n & + & by_n \\
y_{n+1} &=& cx_n & + & dy_n.
\end{aligned}
```

Pentru că $A$ are două valori proprii va avea deci doi vectori proprii liniar independenți $\vecv_1$ corespunzător lui $\lambda_1$ și respectiv 
$\vecv_2$ corespunzător lui $\lambda_2.$ Dacă valoarea inițială are descompunerea 
```math
\left(\begin{array}{cc} x_0 \\ y_0 \end{array}\right)=z_1\vecv_1+z_2\vecv_2,
```
atunci
```math
A^n\left(\begin{array}{cc} x_0 \\ y_0 \end{array}\right)=\lambda_1^n z_1\vecv_1+\lambda_2^nz_2\vecv_2,
```
deci ne regăsim în situația pe care am studiat-o mai sus. Prin urmare comportamentul sistemului este dat de valorile proprii, mai precis de comparația modulului acestora cu $1$ și de semnul lor. Mai precis dacă numim valoarea proprie cea mai mare în modul **dominantă**. Presupunem că dominantă este $\lambda_2=:\lambda_d$. Avem 
  - dacă $|\lambda_d| > 1$ atunci originea este un echilibru instabil;
  - dacă $|\lambda_d| < 1$, atunci originea este un echilibru stabil;
  - dacă $|\lambda_d| = 1$, sistemul este stabil, dar vectorul propriu corespunzător lui $\lambda_d$ este conservat, deci s-ar 
  putea să conveargă la un echilibru netrivial.

## Valori proprii complexe

Să presupunem acum că matricea $A$ are două valori proprii complexe. Acestea vor fi conjugate $\lambda_1=\alpha+i\beta$ 
și $\lambda_2=\alpha-i\beta.$ Considerăm acum matricea noastră ca acționînd pe $\mathbb{C}^2.$ Fie $\vecv$ un vector propriu
pentru $\lambda_1$. Acest vector este de forma $\vecv=\vecv_1+i\vecv_2$, iar $\vecv_1-i\vecv_2$ este vector propriu pentru $\lambda_2$. Avem
```math
A\vecv=A\vecv_1+iA\vecv_2=\lambda_1\vecv
```
Egalăm părțile reale între ele și pe cele imaginare și obținem
```math
\begin{aligned}
A\vecv_1 & = & \alpha\vecv_1 & -& \beta\vecv_2 \\
A\vecv_2 & = & \beta\vecv_1 & +& \alpha\vecv_2.
\end{aligned}
```
Deci în raport cu baza $\left\{\vecv_1,\vecv_2\right\}$ $A$ devine 
```math
\left(\begin{array}{cc} \alpha & \beta \\
-\beta & \alpha \end{array}\right).
```
În funcție de modulul valorii proprii $\lambda_1$, care este egal cu cel al lui $\lambda_2$ avem
  - dacă $|\lambda_1| > 1$ traiectoria este o spirală care pleacă din origine, deci originea este un echilibru instabil;
  - dacă $|\lambda_1| < 1$ traiectoria este o spirală care merge în origine, deci originea este un echilibru stabil;
  - dacă $|\lambda_1| = 1$ traiectoria este o rotație în jurul originii.

Mai precis în ultimul caz vom avea ``\lambda=\cos(\theta)+i \sin(\theta)``. Dacă ``\theta`` este de forma ``\frac{2\pi}{n}`` atunci se obține o orbită periodică de perioadă ``n``, iar dacă nu este de această formă atunci orbita va fi haotică, dar limitată la o coroană circulară în jurul originii.

Vom ignora cazul matricelor care nu au două valori proprii distincte, deoarece acestea nu apar în aplicații fiind cazuri 
instabile, deci printr-o perturbare mică devin diagonalizabile. 	

## Exemple

Revenim la modelul urșilor negri. Avem 
```math
\begin{aligned}
J_{n+1} & = & 0.65 J_n & + & 0.5A_n \\
A_{n+1} & = & 0.25 J_n & + & 0.9A_n.
\end{aligned}
```
Matricea sistemului este
```math
M=\left(\begin{array}{cc}0,65 & 0,5 \\ 0,25 & 0,9\end{array}\right).
```
Calculăm valorile proprii și obținem $1,15$ și $0.4$. Cele două valori proprii au $|\lambda_1| > 1$ și $|\lambda_2| < 1$. Originea este deci un punct șa.

Calculăm vectorii proprii corespunzători. Prima direcție proprie satisface $A=J$, iar cea de-a doua $A=-J/2.$
"""

# ╔═╡ 4cb73eca-7829-48bc-b72c-f2d065adb5e4
md"""
Deci originea este un punct șa, direcția corespunzătoare valorii proprii $0,4$ va fi atrasă de origine. 

Modelul pentru un an prost este 
```math
\begin{aligned}
J_{n+1} & = & 0.5 J_n & + & 0.4A_n \\
A_{n+1} & = & 0.1 J_n & + & 0.8A_n.
\end{aligned}
```

Matricea sistemului este 
```math
M=\left(\begin{array}{cc}0,5 & 0,4 \\ 0,1 & 0,8\end{array}\right)
```
cu valorile proprii $0,4$ și $0,9$. Deci originea este un nod stabil. Direcțiile proprii sînt $A=-J/4$ și respectiv $A=J$.
"""

# ╔═╡ 1d1f06b8-b994-4c25-9b32-8941072ee6c2
begin
	let
M=[0.5 0.4;0.1 0.8]
println(eigvals(M))
println(eigvecs(M))
init1 = [10,5]
init2 = [10,-5]
init3 = [5,20]
init4 = [10,30]
function orbits(init,iters)
    xs = []
    ys = []
    val = init
    push!(xs,val[1])
    push!(ys,val[2])
    for i in 1:iters
        val = M*val
        push!(xs,val[1])
        push!(ys,val[2])
    end
    xs,ys
end
n = 15
xs1,ys1 = orbits(init1,n)
xs2,ys2 = orbits(init2,n)
xs3,ys3 = orbits(init3,n)
xs4,ys4 = orbits(init4,n)
fig4 = Figure(size=(500,500))
ax1 = Axis(fig4[1,1])
scatter!(ax1,xs1,ys1)
scatter!(ax1,xs2,ys2)
scatter!(ax1,xs3,ys3)
scatter!(ax1,xs4,ys4)
fig4
	end
end

# ╔═╡ 4a345661-4192-4773-9ff9-328b87b36621
md"""
Un al treilea exemplu va fi de comportament oscilatoriu. 
```math
\begin{aligned}
J_{n+1} & = & 0.1 J_n & + & 1.4A_n \\
A_{n+1} & = & 0.4 J_n & + & 0.2A_n.
\end{aligned}
```

Matricea sistemului este 
```math
M=\left(\begin{array}{cc}0,1 & 1,4 \\ 0,4 & 0,2\end{array}\right)
```
cu valorile proprii $-0,6$ și $0,9$. De asemenea echilibrul este stabil, dar drumul spre echilibru este oscilant. Direcțiile proprii sînt $A=-J/2$ și $A=7J/4$.
"""

# ╔═╡ 9e254e84-2cd2-4ebf-9f2b-7bfe33693f4f
begin
	let
		M=[0.1 1.4;0.4 0.2]
println(eigvals(M))
println(eigvecs(M))
init1 = [10,5]
init2 = [10,-5]
init3 = [5,20]
init4 = [10,30]
function orbits(init,iters)
    xs = []
    ys = []
    val = init
    push!(xs,val[1])
    push!(ys,val[2])
    for i in 1:iters
        val = M*val
        push!(xs,val[1])
        push!(ys,val[2])
    end
    xs,ys
end
n = 15
xs1,ys1 = orbits(init1,n)
xs2,ys2 = orbits(init2,n)
xs3,ys3 = orbits(init3,n)
xs4,ys4 = orbits(init4,n)
fig5= Figure(size=(500,500))
ax1 = Axis(fig5[1,1])
scatter!(ax1,xs1,ys1)
scatter!(ax1,xs2,ys2)
scatter!(ax1,xs3,ys3)
scatter!(ax1,xs4,ys4)
fig5
	end
end

# ╔═╡ 7fe65162-6639-4753-ad86-c2c5711dd0ae
md"""
## Echilibre neutre
Modelul $S-I$ pe care l-am studiat mai demult era de forma
```math
\begin{aligned}
\left(\begin{array}{cc} S_{n+1} \\ I_{n+1}\end{array}\right)=\left(\begin{array}{cc} 1-a & b \\ a & 1-b\end{array}\right)
\left(\begin{array}{cc} S_{n} \\ I_{n}\end{array}\right).
\end{aligned}
```
Ecuația sa caracteristică este
```math
\lambda^2-(2-a-b)\lambda+(1-a-b)=0.
```

Rădăcinile sale sînt $\lambda_1=1$ și $\lambda_2=1-a-b$. Direcțiile proprii sînt $S=bI/a$ și respectiv $S=-I$. Luăm spre exemplu 
$a=0,1$ și $b=0,2$.

Se vede că iterațiile sînt atrase de direcția corespunzătoare lui $1$, paralel cu cealaltă direcție proprie.
"""

# ╔═╡ db341cee-38e7-4967-9d91-64ad9c189521
begin
	let
		M=[0.9 0.2;0.1 0.8]
println(eigvals(M))
println(eigvecs(M))
init1 = [10,5]
init2 = [10,-5]
init3 = [5,20]
init4 = [10,30]
function orbits(init,iters)
    xs = []
    ys = []
    val = init
    push!(xs,val[1])
    push!(ys,val[2])
    for i in 1:iters
        val = M*val
        push!(xs,val[1])
        push!(ys,val[2])
    end
    xs,ys
end
n = 15
xs1,ys1 = orbits(init1,n)
xs2,ys2 = orbits(init2,n)
xs3,ys3 = orbits(init3,n)
xs4,ys4 = orbits(init4,n)
fig6= Figure(size=(500,500))
ax1 = Axis(fig6[1,1])
scatter!(ax1,xs1,ys1)
scatter!(ax1,xs2,ys2)
scatter!(ax1,xs3,ys3)
scatter!(ax1,xs4,ys4)
fig6
	end
end

# ╔═╡ 7339dbea-acc2-412b-81f5-4ca36ff7bcba
md"""
# Sisteme de ecuații diferențiale liniare
``\newcommand{\vecx}{\overrightarrow{x}}``
Am văzut mai demult că pentru un sistem liniar decuplat
```math
\begin{aligned}
x' & = & ax \\
y' & = & dy
\end{aligned}
```
soluția este $\left(x(t),y(t)\right)=\left(x_0e^{at},y_0e^{dt}\right).$
Semnul coeficienților dă stabilitatea:
  - ``a > 0`` și $d > 0$ atunci originea este un echilibru instabil;
  - ``a < 0`` și $d > 0$ sau invers, atunci originea este un punct șa;
  - ``a < 0`` și $d < 0$ atunci originea este un echilibru stabil.
  
Trecem acum la cazul general. Pentru început presupunem că matricea sistemului este
```math
A=\left(\begin{array}{cc} a & b \\ c & d \end{array}\right)
```
cu două valori proprii distincte $\lambda_1$ și $\lambda_2$ cu vectorii proprii $\vecv_1$, respectiv $\vecv_2$. Orice soluție 
a sistemului de ecuații diferențiale este de forma
```math
\left(x(t),y(t)\right)=c_1e^{\lambda_1t}\vecv_1+c_2e^{\lambda_2t}\vecv_2,
```
constantele fiind determinate de valorile inițiale. Să presupunem că valoarea proprie dominantă este $\lambda_1$. Aunci soluția 
poate fi rescrisă
```math
\left(x(t),y(t)\right)=e^{\lambda_1t}\left(c_1\vecv_1+c_2e^{(\lambda_2-\lambda_1)t}\vecv_2\right).
```
Coeficientul lui $\vecv_2$ tinde la $0$. Deci comportamentul la infinit va fi determinat de $\lambda_1$. În rezumat:
 - Dacă $\lambda_1 < 0$ atunci originea este un echilibru stabil;
 - Dacă $\lambda_1 > 0$ și $\lambda_2 < 0$, atunci avem un punct șa;
 - Dacă ambele valori proprii sînt pozitive, originea este pur instabilă (sursă).

Studiem acum cazul a două valori proprii complexe $\lambda_{1,2}=\alpha\pm\beta i$. Atunci am văzut că există o bază (formată
din părțile reală și imaginară a vectorului propriu corespunzător uneia din valorile proprii) în raport cu care matricea 
sistemului este
```math
M=\left(\begin{array}{cc}\alpha & \beta \\ -\beta & \alpha\end{array}\right).
```
Din nou ne gîndim că matricea $M:\mathbb{C}^2 \to \mathbb{C}^2$. Soluția generală a sistemului este
```math
\vecx(t)=e^{(\alpha+i\beta)t}\vecv=e^{(\alpha+i\beta)t}(\vecv_1+i\vecv_2)=e^{\alpha t}\left(\cos(\beta t) + i\sin(\beta t)\right)
(\vecv_1+i\vecv_2).
```
Va trebui să revenim în spațiul real. O soluție a sistemului va fi o combinație liniară a părților reală și imaginară din 
expresia de mai sus:
```math
\left(x(t),y(t)\right)=e^{\alpha t}\left[c_1\left(\cos(\beta t)\vecv_1-\sin(\beta t)\vecv_2\right)+c_2\left(\sin(\beta t)\vecv_1
+\cos(\beta t)\vecv_2\right)\right].
```

În concluzie ambele valori proprii au aceeași parte reală $\alpha$ care determină comportamentul sistemului:
 - Dacă $\alpha < 0$ orbitele sistemului sînt spirale care se duc în origine;
 - Dacă $\alpha > 0$ orbitele sistemului sînt spirale care se îndepărtează de origine;
 - Dacă $\alpha=0$ orbitele sînt închise în jurul originii, spunem că originea este un centru.
"""

# ╔═╡ d4689b6d-94b0-417c-a3d6-d58cbc86769e
md""" a1 = $(@bind a1 NumberField(-2.0:0.1:1.5))"""

# ╔═╡ a1caceb6-ad12-4c03-b40d-80d5f066a6ed
md""" b1 = $(@bind b1 NumberField(-2.0:0.1:1.5))"""

# ╔═╡ 2f3c3a0b-850c-4caa-a664-2d8d9aa1e50c
md""" c1 = $(@bind c1 NumberField(-2.0:0.1:1.5))"""

# ╔═╡ fc4597d9-1455-4c8b-9c7e-329ae0668e70
md""" d1= $(@bind d1 NumberField(-2.0:0.1:1.5))"""

# ╔═╡ b6c0500a-9951-4186-996e-be687aaf223d
begin
	function gen!(du,u,p,n)
		a,b,c,d = p
		du[1] = a*u[1]+b*u[2]
		du[2] = c*u[1]+d*u[2]
		return nothing
	end
	initv = [3.0,8.0]
	p1 = [a1,b1,c1,d1]
	sisdiff = CoupledODEs(gen!, initv, p1)
	total_time = 40.5
    sampling_time = 0.02
    Y, ti = trajectory(sisdiff, total_time; Ttr = 2.2, Δt = sampling_time)
	fig2=Figure(size=(1500,500))
	ax2=Axis(fig2[1,1])
	ax20=Axis(fig2[1,2])
	scatter!(ax2,Y)
	for col in columns(Y)
		scatter!(ax20,col)
	end
	fig2
end

# ╔═╡ Cell order:
# ╟─e2576524-c9f3-11f0-180f-1b51b8a5b421
# ╟─71b64ca1-62e8-431b-bdf3-58d12810d900
# ╟─d92c50aa-19c6-4486-9612-4a04ed932d85
# ╟─fced5079-d314-45d4-a956-f765df74b99c
# ╟─d3721a4f-9d81-4220-86dd-82d70b181b64
# ╟─aac19b27-6864-444f-81fd-f2317eaab47b
# ╟─d00845bc-f88b-4afa-a7ab-fec17b35c033
# ╟─6391e7e3-c7a7-43b4-b6b7-4a5ca4d6b398
# ╟─252de602-5ce3-48e7-8321-3da20f88a59a
# ╟─5c4a5a48-8cd4-4444-865d-1b8c5517a4fd
# ╟─2a08f8d0-a81c-40f1-86c0-9984f1a27cea
# ╟─4cb73eca-7829-48bc-b72c-f2d065adb5e4
# ╟─1d1f06b8-b994-4c25-9b32-8941072ee6c2
# ╟─4a345661-4192-4773-9ff9-328b87b36621
# ╟─9e254e84-2cd2-4ebf-9f2b-7bfe33693f4f
# ╟─7fe65162-6639-4753-ad86-c2c5711dd0ae
# ╟─db341cee-38e7-4967-9d91-64ad9c189521
# ╟─7339dbea-acc2-412b-81f5-4ca36ff7bcba
# ╟─d4689b6d-94b0-417c-a3d6-d58cbc86769e
# ╟─a1caceb6-ad12-4c03-b40d-80d5f066a6ed
# ╟─2f3c3a0b-850c-4caa-a664-2d8d9aa1e50c
# ╟─fc4597d9-1455-4c8b-9c7e-329ae0668e70
# ╠═b6c0500a-9951-4186-996e-be687aaf223d
