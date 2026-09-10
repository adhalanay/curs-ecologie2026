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

# ╔═╡ e93b92b0-b81a-11f0-25f3-f5adb9a75a31
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

# ╔═╡ 9b1341dd-1725-4e58-b63f-c9c78452c8f4
begin
	using Roots
a1=0.01
f(x)=(a1+x^2)/(1+x^2)
g(x)=0.4*x
h(x)=f(x)-g(x)
ts=range(0.0,2.5,length=50)
ys3=[f(t) for t ∈ ts]
ys4=[g(t) for t ∈ ts]
fig8=Figure()
ax8=Axis(fig8[1,1])
lines!(ts,ys3,linewidth=5)
lines!(ts,ys4,color=:red,linewidth=5)
sol1=find_zero(h,(0,0.2))
sol2=find_zero(h,(0.3,0.6))
sol3=find_zero(h,(1,3))
scatter!([sol1,sol2,sol3],[f(sol1),f(sol2),f(sol3)])
println(sol1)
println(sol2)
println(sol3)
fig8
end

# ╔═╡ f9ad9e24-7143-43f7-a5aa-1133914fbe29
WidthOverDocs()

# ╔═╡ ac2735d2-078c-4970-8f5f-68d930dd4a1d
md"""
# Echilibre în mai multe dimensiuni

## Sisteme liniare

 Am văzut că în cazul ecuațiilor scalare avem doar două posibilități: soluțiile cresc la nesfîrșit sau se stabilizează către o valoare de echilibru. Pentru sisteme multidimensionale pot apărea mai multe situații după cum vom vedea.
  
  Pentru sisteme de ecuații diferențiale punctul de echilibru este un punct în care toate componentele sînt staționare. Echivalent, toate funcțiile din membrul drept trebuie să se anuleze. 
  
  Fie sistemul de ecuații diferențiale:
```math
\begin{aligned}
x_1' & = & f_1(x_1,\dots,x_n) \\
x_2' & = & f_2(x_1,\dots,x_n) \\
\vdots & \vdots & \vdots \\
x_n' & = & f_n(x_1,\dots,x_n)
\end{aligned}
```

**Un punct de echilibru pentru sistemul de ecuații diferențiale este un punct** $(x_1^*,\dots,x_n^*)$ **pentru care** $f_1(x_1^*,\dots,x_n^*)=0,\dots,f_n(x_1^*,\dots,x_n^*)=0.$

Primele exemple de echilibre in 2 dimensiuni se obțin din echilibre în dimensiune 1. După cum am văzut din două spații 1 dimensionale $V$ și $W$ putem forma un spațiu 2 dimensional $V \times W$ format din perechile $(x,y)$ cu $x \in V$ și $y \in W$.
Geometric putem considera $V$ și $W$ ca cele 2 axe ale planului.

Să presupunem că avem ecuațiile scalare $x'=x$ și $y'=y$. Acestea au amîndouă cîte un echilibru instabil în $0$. Compunerea lor duce la un echilibru instabil în $(0,0)$ numit **sursă** sau **nod instabil**.
"""

# ╔═╡ cdfb91f4-3ee9-41bb-80e5-376096d9afc7
begin
	f1(x) = Point2f(x[1],x[2])
	fig1 = Figure()
	ax1 = Axis(fig1[1,1],title="x\'=x,y'=y")
	xs1 = LinRange(-2, 2, 10)
    ys1 = LinRange(-2, 2, 10)
	us1 = [-y for x ∈ xs1, y ∈ ys1]
	vs1 = [x for x ∈ xs1, y ∈ ys1]
	strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax1,f1,xs1,ys1,colormap = :magma,linewidth = 1)
	fig1
end

# ╔═╡ 37160ed0-8527-43da-99ba-c202e7caacf3
md"""
Combinînd două echilibre stabile obținem un echilibru stabil în punctul $(0,0)$ numit **scurgere** sau **nod stabil**.
"""

# ╔═╡ beacfa71-a080-4d5e-ad7a-2d25dd435eab
begin
	f2(x) = Point2f(-x[1],-x[2])
	fig2 = Figure()
	ax2 = Axis(fig2[1,1],title="x\'=x,y\'=-y")
	# xs1 = LinRange(-2, 2, 10)
    # ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax2,f2,xs1,ys1,colormap = :magma,linewidth = 1)
	fig2
end

# ╔═╡ 3fea14fd-b04a-4e73-a47c-a894d92c0d01
md"""
În cazul dimensiunii $2$ apare o situație nouă, atunci cînd combinăm un echilibru stabil cu unul instabil. Atunci se obține un punct **șa**. Pentru un astfel de punct avem o direcție stabilă și una instabilă. Un punct se va mișca în direcția axei instabile și se va îndepărta de cea stabilă. Deci singura posibilitate să atingă echilibrul este ca punctul să se plece exact de pe axa stabilă.
"""

# ╔═╡ a7fff0c6-3ac9-4053-a024-8abb5b273efd
begin
	f3(x) = Point2f(x[1],-x[2])
	fig3 = Figure()
	ax3 = Axis(fig3[1,1],title="x\'=x,y\'=-y ")
	# xs1 = LinRange(-2, 2, 10)
    # ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax3,f3,xs1,ys1,colormap = :magma,linewidth = 1)
	fig3
end

# ╔═╡ a545145a-8c51-4243-88b2-5dbdc672a804
md"""
Un alt fenomen nou care apare în dimensiune 2 este unul care nu provine din compunerea a două sisteme $1-$dimensionale. Aceste puncte de echilibru implică mișcări de rotație în jurul lor. Ne reamintim sistemul arcului:
```math
\begin{aligned}
x' & = & v \\
v' & = & -x -v. 
\end{aligned}
```
Echilibrul obținut se numește **spirală stabilă**.
"""

# ╔═╡ 65fa6491-541c-4728-8b52-9072c26b5e0d
begin
	f4(x) = Point2f(x[2],-x[1]-x[2])
	fig4 = Figure()
	ax4 = Axis(fig4[1,1],title="x\'=y,y\'=-x-y ")
	# xs1 = LinRange(-2, 2, 10)
    # ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax4,f4,xs1,ys1,colormap = :magma,linewidth = 1)
	fig4
end

# ╔═╡ d575e59c-3647-41a0-9503-7a8ffd366572
md"""
Dacă modificăm ecuația arcului astfel:
```math
\begin{aligned}
x' & = & v \\
v' & = & -x +v. 
\end{aligned}
```
Obținem o **spirală instabilă**:
"""

# ╔═╡ 5f384a87-31b3-4f33-96e6-9f1ceeeb3715
begin
	f5(x) = Point2f(x[2],-x[1]+x[2])
	fig5 = Figure()
	ax5 = Axis(fig5[1,1],title="x\'=y,y\'=-x+y ")
	# xs1 = LinRange(-2, 2, 10)
    # ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax5,f5,xs1,ys1,colormap = :magma,linewidth = 1)
	fig5
end

# ╔═╡ fc4212e8-ee8d-497f-8672-5efd1c9ae6c1
md"""
În sfîrșit mai există un tip de echilibru în plan, pe care l-am mai întîlnit atunci cînd am studiat arcul fără frecare:
```math
\begin{aligned}
x' & = & v \\
v' & = & -x.
\end{aligned}
```
Acest tip de echilibru se numește **centru**.
"""

# ╔═╡ 924ffee8-d3b3-4de5-b0c2-1b7fa4970e53
begin
	f6(x) = Point2f(x[2],-x[1])
	fig6 = Figure()
	ax6 = Axis(fig6[1,1],title="x\'=y,y\'=-x ")
	# xs1 = LinRange(-2, 2, 10)
    # ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax6,f6,xs1,ys1,colormap = :magma,linewidth = 1)
	fig6
end

# ╔═╡ 1938dc9f-56d5-45e6-a419-ec73dd0218ae
md"""
## Sisteme neliniare

Considerăm acum exemple de sisteme neliniare, adică sisteme de ecuații diferențiale pentru care membrul drept este neliniar. Acestea au de obicei mai multe echilibre care au adesea comportamente diferite. Să luăm ca exemplu un model cu două populații (căprioare și elani), notate cu $D$ și $M$ care concurează pentru resurse. Presupunem că în absența concurenței populația de căprioare va avea o rată de creștere de $3$, iar cea de elani va avea o rată de creștere de $2$. Competiția are loc în cadrul aceleiași specii producînd termenii de aglomerare $-D^2$ și $-M^2$ și de asemenea concurează cu cealaltă specie, dar impactul căprioarelor asupra elanilor este $0,5$, iar a elanilor asupra căprioarelor este mai puternică, de $1$. Ecuațiile care se obțin sînt
```math
\begin{aligned}
D' & = & 3D & - & D^2 & - & MD \\
M' & = & 2M & - & M^2 & - & 0,5 MD.
\end{aligned}
```
Echilibrele sale $(D^*,M^*)$ sînt $(0,0)$, numit și **echilibrul trivial**, $(0,2)$, $(3,0)$ și $(2,1)$.

Pentru studiul stabilității acestor echilibre nu avem deocamdată nici un instrument matematic, dar putem să ne facem o idee din cimpul de vectori asociat.
"""

# ╔═╡ edc58e8e-dec7-44e9-8423-0de390cab9db
begin
	f7(x) = Point2f(3*x[1]-x[1]^2-x[1]*x[2],2*x[2]-x[2]^2-0.5*x[1]*x[2])
	fig7 = Figure()
	ax7 = Axis(fig7[1,1],title=L"x=3x-x^2-xy',y=2y-y^2-0,5xy' ",limits = (-0.1,5,-0.1,5))
	xs2 = LinRange(0, 5, 20)
    ys2 = LinRange(0, 5, 20)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax7,f7,xs2,ys2,colormap = :magma,linewidth = 1)
	scatter!(ax7,(0,2),markersize=20)
	scatter!(ax7,(3,0),markersize=20)
	scatter!(ax7,(0,0),markersize=20)
	scatter!(ax7,(2,1),markersize=20)
	fig7
end

# ╔═╡ 9622372d-2da7-4775-af21-f543086d7594
md"""
Se vede că sînt două puncte șa ($(0,2)$ și $(3,0)$), unul este instabil $(0,0)$ și unul este un punct echilibru stabil ($(2,1)$). Vom vedea că toate echilibrele sistemelor neliniare sînt de tipul celor pe care le-am văzut mai înainte.

## Metode empirice pentru studiul stabilității echilibrelor

O metodă foarte eficientă este metoda **null-clineolor.** Considerăm sistemul de ecuații
```math
\begin{aligned}
x' & = & f(x,y) \\
y' & = & g(x,y)
\end{aligned}
```
 Pentru
Luăm acum curbele definite de ``x'=0``, respectiv ``y'=0``. Adică ``f(x,y)=0`` și ``g(x,y) = 0``.Pentru sistemul nostru avem ``D-``null-clina ``3-D-M=0`` și ``M-``null-clina ``2-M-0,5D=0``. Alături de ele mai avem și ``D=0``, respectiv ``M=0``.
"""

# ╔═╡ b75f3d68-5c68-49c0-b4b6-3b3bc20a6882
begin
	n1(x) = 3-x
	n2(x) = 2-0.5*x
	ds = [n1(x) for x in xs2]
	ms = [n2(x) for x in xs2]
	lines!(ax7,[(0,0),(0,5)],linewidth=4)
	lines!(ax7,[(0,0),(5,0)],linewidth=4)
	lines!(ax7,xs2,ds,linewidth=4)
	lines!(ax7,xs2,ms,linewidth=4)
	fig7
end

# ╔═╡ b420830b-ccfe-41e1-9898-69f642a2b9a5
md"""
Aceste drepte împart planul în patru sectoare. În fiecare sector toți vectorii care dau evoluția sistemului au vîrful în aceeași direcție. Spre exemplu vedem că echilibrul din centru este stabil (săgețiile din toate sectoarele arată spre el). Toate celelalte
echilibre au vectori care arată dinspre ele, deci sînt instabile. 

Putem folosi null-clinele și pentru a studia comportamentul cîmpului de vectori. Începem prin a desena null-clinele. Apoi pe fiecare null-clină vom lua cîte un punct de fiecare parte a punctului de echilibru și vom evalua ecuația null-clinei în aceste puncte,
dacă obținem o valoare pozitivă atunci vectorul are vîrful în sus, iar dacă e negativă atunci îl are în jos.

Folosirea null-clinelor își arată valoarea atunci cînd lucrăm cu sisteme pentru care valorile parametrilor sînt necunoscute. Chiar și în această situație putem să avem o imagine aproximativă a comportamentului cîmpurilor de vectori. 

## Bazine de atracție

Să presupunem că avem un sistem cu mai multe puncte de echilibru. Pentru echilibrele stabile există o regiune din plan pentru care toate valorile inițiale aflate acolo se apropie de echilibru. Acesta este **bazinul de atracție** al respectivului echilibru. Spre exemplu pentru ecuația cu efectul Allee bazinul de atracție al lui $X=0$ este format din toate populațiile mai mici decît $a$, iar pentru $X=k$ toate populațiile mai mari decît $a$. Numele provine din geografie, fiind inspirat de noțiunea de bazin hidrografic. 

### Exemplu : operonul lactozei (lac operon)

În multe procese biologice conceptul de comutator și cel de prag joacă un rol foarte mare. De exemplu producția unei enzime sau hormon se activează atunci cînd anumite semnale trec un anumit prag; celulele care se dezvoltă după un anumit punct devin au un "destin" ireversibil; activitatea electrică (tensiunea) din neuroni sau celule cardiace este stabilă pînă cînd un anumit stimul o face să treacă pragul.

Un exemplu clasic este dat de bacteria *Escherichia coli*. Această bacterie poate transforma lactoza în energie, iar pentru introducerea lactozei în celulă este nevoie de o proteină, numită permeaza de lactoză. Producția de permează este costisitoare pentru celulă, deci este produsă doar atunci cînd concentrația de lactoză este suficient de mare. 

Vom modela acest proces printr-o ecuație diferențială. Fie $X$ concentrația de lactoză în celulă. Atunci X'=lactoză importată - lactoză metabolizată. Rata lactozei importate este proporțională cu cantitatea de permează. Presupunem că există un nivel constant de bază al producției de permează $a.$ Cînd concentrația de lactoză crește și producția de permează crește, dar la niveluri foarte mari de lactoză, se stabilizează la un anumit nivel:
O astfel de funcție este 
```math
f(x)=\frac{a+x^2}{1+x^2},
```
numită **sigmoid.**
Pentru degradarea lactozei vom folosi o funcție liniară $g(x)=kx$. Grupînd obținem
```math
x'=\frac{a+x^2}{1+x^2}-kx.
```
Ca să studiem echilibrele acestei ecuații vom folosi o metodă grafică. Anume desenăm graficele celor două funcții $f$ și $g$ și vom vedea unde se intersectează, apoi vom vedea semnul acestora.
"""

# ╔═╡ 01015c51-df72-4f44-8ea7-3df352bf63fe
md"""
Acolo unde graficul albastru este deasupra celui roșu $x' > 0$, iar unde cel roșu este deasupra celui albastru $x'<0$. Se vede de aici că cel mai mic și cel mai mare echilbru sînt stabile, iar echilibrul din mijloc este instabil. Echilibrul din mijloc joacă rolul de întrerupător: pentru valori mai mici decît el, concentrația de lactază din celulă scade, dar secreția de permează nu crește. Pentru valori mai mari decît acest echilibru, crește producția de permează, și deci și cantitatea de lactoză absorbită.
"""

# ╔═╡ 24f3153b-9fb7-469c-8722-6a4184b12a59
md"""
## Sisteme discrete liniare

Considerăm sistemul format din doua sisteme liniare cu echilbru instabil
```math
\begin{aligned}
x_{t+1} & = & 2 x_t \\
y_{t+1} & = & 2 y_t.
\end{aligned}
```
Punctul de echlibru este $(0,0)$ si este un echlibru instabil.
"""

# ╔═╡ 4f8714f2-b709-43ff-9378-e5c3da1d7015
md"""
Consideram acum sistemul format din doua sisteme liniare cu echilbru stabil
```math
\begin{aligned}
x_{t+1} & = & 0.2 x_t \\
y_{t+1} & = & 0.8 y_t.
\end{aligned}
```
Punctul de echlibru este $(0,0)$ si este un echlibru stabil.
"""

# ╔═╡ eaaa9fb7-d6f9-4643-89b5-063811cf7808
md"""
Daca avem sistemul 
```math
\begin{aligned}
x_{t+1} & = & 1.1 x_t \\
y_{t+1} & = & 0.9 y_t.
\end{aligned}
```

atunci echilibrul $(0,0)$ a fi un punct sa.

Alte posibilități se văd mai jos:
"""

# ╔═╡ efc40691-5a54-4171-8640-bda2980cb5c7
md""" a = $(@bind a NumberField(0:0.1:1.5))"""

# ╔═╡ 9c1bd8f9-4442-41c7-889e-0758c14f2dc9
md""" b = $(@bind b NumberField(0:0.1:1.5))"""

# ╔═╡ 9ba25a30-26af-41a2-9833-9a0e2ba4753c
md""" c = $(@bind c NumberField(0:0.1:1.5))"""

# ╔═╡ e13c93f7-8914-4431-8a55-299440f94593
md""" d = $(@bind d NumberField(0:0.1:1.5))"""

# ╔═╡ 83548534-8a28-4b58-bf3b-fa8ac9a98bca
md""" numărul de pași: $(@bind timp PlutoUI.Slider(1:200,show_value=true))"""

# ╔═╡ 07ac56fe-c8c2-4af4-91eb-4c87e0ce8a21
begin
	using DynamicalSystems
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
	fig9 = Figure(size=(1000,500))
	ax9 = Axis(fig9[1,1],title="Sisteme liniare cu parametrii $a, $b, $c, $d")
	scatter!(ax9,u0[1],u0[2],markersize=15,color=:red)
	scatter!(ax9,X)
	ax10 = Axis(fig9[1,2], title="Seriile temporale")
	for v in columns(X)
		scatter!(ax10,v)
	end
	fig9
end


# ╔═╡ db397d50-fe1d-48c4-8276-d255a572448a
md"""
Vedem că un sistem poate avea diferite alte comportamente: spirale stabile sau instabile sau orbite periodice. Spre exemplu pentru orice valoare inițială nenulă sistemul 
```math
\begin{aligned}
x_{n+1} & = & -y_n \\
y_{n+1} & = & x_n.
\end{aligned}
```
are orbită periodică de perioadă ``4``. De asemenea pentru anumite valori ale parametrilor fiecare componentă se va stabiliza la cîte un punct fix. Ce nu apare în cazul sistemelor discrete liniare este comportamentul haotic. 

## Modele Leslie

În cazul modelelor biologice, sistemele liniare se mai numesc modele **matriciale**. 

O clasă importantă de astfel de modele sînt modelele **Leslie** sau modele compartimentale pentru dinamica populațiilor. Astfel împărțim populația în diferite clase (de exemplu stadii de dezvoltare) și descriem procentul de indivizi care trec dintr-o clasă în alta ținînd cont de natalitate și mortalitate.

Să presupunem că avem o populație de urși negrii americani (Ursus americanus) cu două stadii de viață : juvenili și adulți. Matricea sistemului este de forma
```math
M=\left(\begin{array}{cc} 0.65 & 0.5 \\ 0.25 & 0.9\end{array}\right),
```
pe diagonală este procentul de indivizi care rămîn ăn clasă, iar celelalte valori sînt cele ale indivizilor care păresesc clasa. Explicit avem
```math
\begin{aligned}
J_{n+1} & = & 0.65 J_n & + & 0.5A_n \\
A_{n+1} & = & 0.25 J_n & + & 0.9 A_n,
\end{aligned}
```
adică ``65\%`` din juvenili rămîn juvenili și rata de reproducere este ``50\%.``Pentru adulți avem că ``25\%`` din juvenili devin adulți și avem o mortalitate de ``10\%`` în rîndul adulților deci ``90\%`` rămîn pentru anul următor. 

Vedem că numărul de urși va crește cu numărul de juvenili apropiindu-se de cel de adulți.

Dacă presupunem acum că avem un an prost (cu natalitate mică și mortalitate mare) avem acum matricea
```math
M=\left(\begin{array}{cc} 0.5 & 0.4 \\ 0.1 & 0.8\end{array}\right).
```
Modelul prezice că ambele clase se vor duce la $0$.

Pentru a treia situație să presupunem că avem natalitate mare, dar și mortalitate mare:
```math
M=\left(\begin{array}{cc} 0.1 & 1.4 \\ 0.4 & 0.2\end{array}\right).
```
mai precis ``10\%`` din juvenili rămîn juvenili, ``40\%`` devin adulți și restul mor. Natalitatea este de ``1.4`` și doar ``20\%`` din adulți supraviețuiesc pentru anul următor. Dinnou numărul de indivizi se duce spre ``0``, dar de o manieră oscilantă. 

Un caz special este atunci cînd suma elementelor de pe coloane este ``1``. Spunem că matricea ``M`` este **matrice Markov.** În acest caz se vede că valorile tind către un punct de echilibru, dar acesta este **diferit** în funcție de punctul inițial.
"""

# ╔═╡ Cell order:
# ╟─e93b92b0-b81a-11f0-25f3-f5adb9a75a31
# ╠═f9ad9e24-7143-43f7-a5aa-1133914fbe29
# ╟─ac2735d2-078c-4970-8f5f-68d930dd4a1d
# ╟─cdfb91f4-3ee9-41bb-80e5-376096d9afc7
# ╟─37160ed0-8527-43da-99ba-c202e7caacf3
# ╟─beacfa71-a080-4d5e-ad7a-2d25dd435eab
# ╟─3fea14fd-b04a-4e73-a47c-a894d92c0d01
# ╟─a7fff0c6-3ac9-4053-a024-8abb5b273efd
# ╟─a545145a-8c51-4243-88b2-5dbdc672a804
# ╟─65fa6491-541c-4728-8b52-9072c26b5e0d
# ╟─d575e59c-3647-41a0-9503-7a8ffd366572
# ╟─5f384a87-31b3-4f33-96e6-9f1ceeeb3715
# ╟─fc4212e8-ee8d-497f-8672-5efd1c9ae6c1
# ╟─924ffee8-d3b3-4de5-b0c2-1b7fa4970e53
# ╟─1938dc9f-56d5-45e6-a419-ec73dd0218ae
# ╠═edc58e8e-dec7-44e9-8423-0de390cab9db
# ╟─9622372d-2da7-4775-af21-f543086d7594
# ╟─b75f3d68-5c68-49c0-b4b6-3b3bc20a6882
# ╟─b420830b-ccfe-41e1-9898-69f642a2b9a5
# ╟─9b1341dd-1725-4e58-b63f-c9c78452c8f4
# ╟─01015c51-df72-4f44-8ea7-3df352bf63fe
# ╟─24f3153b-9fb7-469c-8722-6a4184b12a59
# ╟─4f8714f2-b709-43ff-9378-e5c3da1d7015
# ╟─eaaa9fb7-d6f9-4643-89b5-063811cf7808
# ╠═efc40691-5a54-4171-8640-bda2980cb5c7
# ╠═9c1bd8f9-4442-41c7-889e-0758c14f2dc9
# ╠═9ba25a30-26af-41a2-9833-9a0e2ba4753c
# ╠═e13c93f7-8914-4431-8a55-299440f94593
# ╠═83548534-8a28-4b58-bf3b-fa8ac9a98bca
# ╠═07ac56fe-c8c2-4af4-91eb-4c87e0ce8a21
# ╟─db397d50-fe1d-48c4-8276-d255a572448a
