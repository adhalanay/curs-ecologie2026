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

# ╔═╡ e0fdf9e0-b175-11f0-2b32-11f3f8375fcd
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

# ╔═╡ f8248b11-863f-4fb8-b49c-9d23018269d2
WidthOverDocs()

# ╔═╡ 56090739-0808-499f-b740-9d0910230637
md"""
# Sisteme discrete

În foarte multe situații fenomenele se petrec la anumite momente de timp. Spre exemplu speciile **semelpare** se reproduc o singură dată în viață, apoi mor. Astfel de animale sînt unele specii din genul Onchorhyncus (somon de Pacific), speciile genului Magicicada ale căror membrii ies din hibernare o dată la 13-17 ani pentru a se reproduce și muri. În regnul vegetal, există gențiana verde (Frasera speciosa) sau planta secolului (Agavi americana) care înfloresc după 20 de ani de viață, apoi  produc semințe și mor. 

Pentru a înțelege evoluția acestui tip de specii avem nevoie de sisteme discrete. Dinamica este de tipul 
```math
x_{t+1}=f(x_t),
```
unde $x_t$ reprezintă starea sistemului la momentul $t$. În acest caz $t$ ia doar valori discrete (ex. $0,1,2,\dots$). Expresia se numește *sistem dinamic discret* sau *ecuație cu diferențe finite*.

O soluție a unei astfel de ecuație este o expresie care ne dă valorile lui $x_t$ pentru orice $t$. De exemplu dacă $x_{t+1}=rx_t$, atunci $x_1=rx_0$, 
$x_2=rx_1=r^2x_0$ și în general $x_t=r^t x_0$. Spre deosebire de cazul diferențial, acum dinamica este scrisă explicit.



## Sisteme exponențiale

Să presupunem că avem o populație de căprioare care crește cu $5\%$ pe an. Dacă notăm cu $x_n$ populația din anul $n$, atunci $x_{n+1}=x_n+0.05x_n=1.05x_n.$ Deci este de tip exponențial. 

Evoluția acestui tip de sisteme dinamice se vede mai jos. Dacă ``r<1`` atunci valorile tind către ``0,`` iar pentru ``r > 1`` valorile se duc spre infinit.
 """

# ╔═╡ c25c69b1-70ae-4603-be76-742b509f63ac
md""" r $(@bind r PlutoUI.Slider(0.01:0.01:1.5,show_value=true))"""

# ╔═╡ b91646bb-6af8-4946-b49b-20cc231f8c24
begin
x0=1.0
f_r(x)=r*x
xsr=[x0]

for t ∈ range(1,50,length=50)
    xr=f_r(last(xsr))
    push!(xsr,xr)
end
fig1=Figure()
ax1 = Axis(fig1[1,1])
scatter!(xsr,color=:red)
fig1
end

# ╔═╡ 56f7a412-d62e-4aa0-a43b-aa3a2a48732d
md"""
Dacă ne situăm în cazul ``r > 1`` vedem că valorile cresc la început lent, apoi din ce în ce mai repede. Ca să vedem cît de repede, calculăm
```math
\Delta_{n+1} x=x_{n+2}-x_{n+1}=r(x_{n+1}-x_n)=r\Delta_n x.
```
Deci rata de modificare (viteza) crește cu rata $r.$ Cînd populația este mică rata de creștere este mică, dar pe măsură ce populația crește, ea crește din ce în ce mai repede.

## Modelul logistic

Să studiem acum modelul logistic discret introdus de May în 1976. Presupunem că avem o populație de insecte care trăiesc un an, depun ouă, apoi mor. Populația la  momentul $n+1$ va fi $x_{n+1}=f(x_n).$ Presupunem că mediul suportă maximum $k$ insecte. Dacă populația este $x_n$ atunci ea va consuma $\frac{x_n}{k}$ resurse, deci resursele neutilizate vor fi $1-\frac{x_n}{k}.$ Aceste resurse neutilizate sînt disponibile pentru a susține noi nașteri, deci rata nașterilor este proporțională cu $1-\frac{x_n}{k}$ cu un factor de proporționalitate $r$. Prin urmare avem
```math
f(x)=r\cdot x\cdot \left(1-\frac{x}{k}\right).
```
Ca de obicei presupunem $k=1$. Obținem sistemul dinamic
```math
x_{n+1}=rx_n(1-x_n).
```
Punctele sale de echilibru sînt ``0`` și ``\frac{r-1}{r}.``Să simulăm acum sistemul pentru diferite valori ale parametrului ``r`` și pentru diferite valori inițiale. 
"""

# ╔═╡ bfcb1cb7-a349-4d77-87cf-1fe973076c1c
md""" valoarea inițială $(@bind x1 PlutoUI.Slider(0.000:0.001:1,show_value=true))"""

# ╔═╡ 425c759d-268f-4bcd-afe9-e737fa55c74d
md""" r $(@bind r1 PlutoUI.Slider(0.000:0.0005:4.0,show_value=true))"""

# ╔═╡ 4aed83f0-30a9-4f9d-80ad-a22b84c83a21
begin
f_l(x)=r1*x*(1-x)
xs1 = [x1]
ts = collect(range(0,400))

for t ∈ range(0,400)
    x_c=f_l(last(xs1))
    push!(xs1,x_c)
end
fig2=Figure(size=(1500,600))
ax2 = Axis(fig2[1,1])
	scatter!(xs1,color=:red)
fig2
end

# ╔═╡ 78d927ea-1b59-4fbf-8f40-c1102aac0a94
md"""
Ce se observă:
  - pentru $0\leq r<1$ orbita traiectoria oricărei valori inițiale se duce la ``0;``
  - pentru ``1\leq r \leq 2`` traiectoria tinde către ``\frac{r-1}{r}`` pentru orice punct inițial (înafară de ``0`` și ``1``);
  - pentru ``2 \leq r \leq 3`` orbitele vor tinde tot la ``\frac{r-1}{r}`` dar mai întîi oscilează în jurul acestuia;
  - pentru ``3 \leq r \leq 1+\sqrt{6} \simeq`` $(1+sqrt(6)) orbitele vor lua alternativ ``2`` valori. Avem o orbită periodică de perioadă 2;
  - pentru ``1+\sqrt{6} \leq r \leq 3.44949``(aproximativ) orbita va avea perioadă ``4`` pentru aproape toate valorile inițiale;
  - pentru `` 3.44949< r \leq 3.56995`` perioada se dublează rapid (avem orbite de perioadă ``8,16,32,`` etc.);
  - pentru ``3.56995 < r \leq 4`` se instalează haosul: orbitele evoluează aparent impredictibil pentru aproape toate valorile lui ``r``. Există și excepții, de exemplu pentru ``r=1+\sqrt{8} \simeq`` $(1+sqrt(8)) traiectoria se stabilizează la una periodică de perioadă ``3``.

Sistemul logistic prezintă toate comportamentele posibile pentru un sistem discret. De asemenea se observă câ înainte de a se stabiliza la un anumit tip de orbită traiectoria are un comportament tranzient. Acesta poate dura un termen relativ lung, deci este important să studiem riguros 
teoria calitativă.
"""

# ╔═╡ 0dae07cc-86b6-4255-989d-a08da7bfe946
begin
f3(x)=(1+sqrt(8))*x*(1-x)
xs3 = [x1]

for t ∈ range(0,6000)
    x3=f3(last(xs3))
    push!(xs3,x3)
end
fig3=Figure(size=(1500,300))
ax3 = Axis(fig3[1,1])
	scatter!(xs3,color=:red)
fig3
end

# ╔═╡ d441ba55-03d2-4cb8-ade3-f97f2d7b42b0
md"""

Ca și în cazul ecuațiilor diferențiale, în general este imposibil să determinăm explicit soluțiile unui sistem discret, așa că sîntem nevoiți să le studiem din punct de vedere calitativ. Ca și în cazul continuu putem considera **puncte de echilibru**. În acest caz este vorba de valori ``\hat{x}`` astfel încît 
```math
\hat{x}=f(\hat{x}).
```

Continuînd analogia cu ecuațiile diferențiale punctele de echilibru pot fi *stabile* sau *instabile*.

Criteriul de stabilitate a echilibrului este 

**Teoremă** Fie sistemul discret $x_{t+1}=f(x_t)$ cu punctul de echilibru $\hat{x}$. Echilibrul este
  - stabil dacă $|f'(\hat{x})| < 1;$
  - instabil dacă $| f'(\hat{x})| > 1.$

Pentru modelul logistic, avem ``f(x)=rx(1-x)``. Am văzut că punctele sale de echilibru, cele care satisfc $x=f(x)$ sînt ``\hat{x}=0`` și ``\hat{x}=\frac{r-1}{r}``. Derivata lui ``f`` este 
```math
f'(x) = r(1-2x).
```
Avem ``f'(0)=r`` și ``f'\left(\frac{r-1}{r}\right)=1-r.``Deci dacă ``0 < r < 1``, ``0`` este echilibru stabil. Dacă $1<r<2$ echilibrul ``\frac{r-1}{r}`` este stabil. Apoi pentru ``r \geq 2`` ambele echilibre sînt instabile.

## Haos

Se vede că sistemul are o comportare neregulată, dar totuși pentru valori apropiate de ``0``, variațiile sînt mici și la fel pentru valori apropiate de ``0.75``. 

Comportamentul sistemului logistic discret este un exemplu de comportament **haotic.**

Noțiunea de haos este foarte folosită cu diferite înțelesuri. Uneori se vorbește de sisteme haotice, deși haosul este un comportament, același sistem putînd avea un comportament haotic sau nu în funcție de valorile parametrilor sau valorile inițiale.

Comportarea unui sistem este haotică dacă este
  - Deterministă
  - Mărginită
  - Neregulată
  - Sensibilitate față de valorile inițiale

Comportarea este deterministă dacă starea sistemului la un moment dat este complet determinată doar de stările precedente.

Mărginirea se referă la faptul că orbitele sistemului rămîn între anumite limite.

Neregularitatea înseamnă că orbitele nu sînt periodice sau nu tind la un punct fix. Adică traiectoria nu se repetă **exact**. Datorită caracterului determinist dacă o orbită se repetă adică se întoarce exact în același punct, atunci obligatoriu este o orbită periodică.

Sensibilitatea se referă la faptul că la mici variații ale valorilor inițiale corespund variații mari ale orbitelor. Pe scurt: orbitele care încep foarte aproape una de alta se îndepărtează foarte mult. 

Riguros sensibilitatea se măsoară astfel: presupunem că avem două orbite ``M_t`` și ``N_t`` corespunzînd valorilor inițiale ``M_0`` și ``N_0.`` Atunci sistemul este sensibil la condițiile inițiale, dacă 
```math
dist(M_t,N_t)=e^{λt}dist(M_0,N_0).
```
Numărul ``λ`` se numește **coeficientul Liapunov** al sistemului.

Schimbarea comportamentului în funcție de valorile parametrului se numește **bifurcație.** În cazul nostru avem de-a face cu bifurcația  de dublare a perioadei. Aceasta este unul din drumurile spre haos.

Mai jos vedem o altă reprezentare a comportamentului sistemului, numită diagramă de bifurcație:
"""

# ╔═╡ 32913023-6499-402a-881f-b368d19282f3
@bind rstart PlutoUI.Slider(0.0:0.0001:4.0,show_value=true)

# ╔═╡ d4b65f9c-3749-4941-ac91-042f0128485b
begin
    let
        rs = LinRange(rstart,4.0, 4000)
        x1 = 0.5
        maxiter = 1000 # maximum iterations
        x = zeros(length(rs), maxiter)
        for k in eachindex(rs)
            x[k,1] = x1 # initial condition
            for j = 1 : maxiter-1
                x[k, j+1] = rs[k] * x[k, j] * (1 - x[k,j])
            end
        end
		fig4=Figure(size=(1500,1000),aspect=DataAspect())
		ax4=Axis(fig4[1,1],title="Diagrama de bifurcație a modelului logistic")
		for j in range(50,0,step=-1)
		    scatter!(rs,x[:,end-j],markersize = 1.5)
		end
		fig4
		# x[:,end-50:end]
    end
end

# ╔═╡ f860cbbe-469c-40a4-8b97-e95555699f92
begin
    let
      rs = LinRange(rstart,4.0, 4000)
      x1 = 0.5
      maxiter = 1000 # maximum iterations
        xric = zeros(length(rs), maxiter)
        for k in eachindex(rs)
            xric[k,1] = x1 # initial condition
            for j = 1 : maxiter-1
                xric[k, j+1] = xric[k, j] * exp(rs[k]*(1-xric[k,j]))
            end
        end
		fig5=Figure(size=(1500,1000),aspect=DataAspect())
		ax5=Axis(fig5[1,1],title="Diagrama de bifurcație a modelului Rickert")
		for j in range(50,0,step=-1)
		    scatter!(rs,xric[:,end-j],markersize = 1.5)
		end
		fig5
		# x[:,end-50:end]
    end
end

# ╔═╡ f7f72587-c1bd-4b2d-8be0-d21158b7f20e
md"""
## Modelele Rickert și Beverton-Holt
Un alt sistem discret care prezintă un comportament haotic este sistemul **Ricker**. Acesta fost introdus în contextul management-ului crescătoriilor de pești. Acest sistem ia în considerare mortalitatea, care în cazul acesta este proporțională cu densitatea indivizilor. Evoluția sa este
```math
x_{n+1} = x_n e^{r\left(1-\frac{x_n}{K}\right)},
```
``K`` este capacitatea de încăracare a mediului. ``e^r`` este rata constantă de reproducere, iar factorul ``e^{-r\frac{x}{K}}`` reprezintă mortalitatea care este proporțională cu densitatea. Cu cît densitatea este mai mare, cu atît mortalitatea este mai mare. 

Se observă că nivelul populației poate depăși ``K``. O versiune simplificată este dată de modelul **Beverton-Holt** în care într-adevăr ``K`` este limita maximă a populației. Acest model este
```math
x_{n+1} = \frac{rx_n}{1+\frac{r-1}{K}x_n}.
```
Se observă că are tot două puncte de echilibru, anume ``0`` și ``K``. În funcție de ``r`` avem
 - dacă ``r > 1`` atunci ``K`` este stabil, iar ``0`` este instabil;
 - dacă `` r < 1`` atunci este invers ``0`` este stabil, iar ``K`` este instabil.

Aceste sisteme sînt de forma ``x_{n+1} = f(x_n)`` graficele funcțiilor (pentru ``K=1``) și simularea sistemului Beverton-Holt sînt mai jos. Se vede că sistemul Beverton-Holt nu are un comportament haotic pentru orice valoare a lui ``r``.
"""

# ╔═╡ 3de65238-17f6-41ed-94ca-986aed400382
@bind r0 PlutoUI.Slider(0:0.1:4.0,show_value=true)

# ╔═╡ c50cb750-8971-4933-a7eb-d61168d4710c
begin
	f1(x) = r0*x*(1-x)
	f2(x) = exp(r0*(1-x))*x
	f4(x) = r0*x/(1+(r0-1)*x)
	bv = [x1]
	xs = LinRange(0.0,1.0,500)
	xs4 = LinRange(0.0,4.0,500)
	xs5 = LinRange(0.0, 10.0, 1000)
	ys1 = [f1(x) for x in xs]
	ys2 = [f2(x) for x in xs4]
	ys3 = [f4(x) for x in xs5]
	for t ∈ range(0,6000)
       bv0=f4(last(bv))
       push!(bv,bv0)
    end
	fig5 = Figure(size=(1500,500))
	ax51 = Axis(fig5[1,1][1,1],title = "Verhulst")
		scatter!(xs,ys1)
	ax52 = Axis(fig5[1,1][1,2], title = "Ricker")
		scatter!(xs4,ys2)
	ax53 = Axis(fig5[1,1][1,3], title = "Beverton-Holt")
		scatter!(xs5,ys3)
	ax54 = Axis(fig5[2,1])
	    scatter!(bv,color = :red)
	fig5
end

# ╔═╡ a9f4d36b-7183-4cf3-afef-96e18bf4c8f3
md"""
Se observă că sistemele cu comportament haotic funcțiile care definesc sistemul au un punct de maxim, iar pentru Beverton-Holt graficul admite asimptotă orizontală.
"""

# ╔═╡ Cell order:
# ╟─e0fdf9e0-b175-11f0-2b32-11f3f8375fcd
# ╟─f8248b11-863f-4fb8-b49c-9d23018269d2
# ╟─56090739-0808-499f-b740-9d0910230637
# ╠═c25c69b1-70ae-4603-be76-742b509f63ac
# ╟─b91646bb-6af8-4946-b49b-20cc231f8c24
# ╟─56f7a412-d62e-4aa0-a43b-aa3a2a48732d
# ╟─bfcb1cb7-a349-4d77-87cf-1fe973076c1c
# ╟─425c759d-268f-4bcd-afe9-e737fa55c74d
# ╟─4aed83f0-30a9-4f9d-80ad-a22b84c83a21
# ╟─78d927ea-1b59-4fbf-8f40-c1102aac0a94
# ╟─0dae07cc-86b6-4255-989d-a08da7bfe946
# ╟─d441ba55-03d2-4cb8-ade3-f97f2d7b42b0
# ╟─d4b65f9c-3749-4941-ac91-042f0128485b
# ╟─32913023-6499-402a-881f-b368d19282f3
# ╟─f860cbbe-469c-40a4-8b97-e95555699f92
# ╟─f7f72587-c1bd-4b2d-8be0-d21158b7f20e
# ╠═3de65238-17f6-41ed-94ca-986aed400382
# ╟─c50cb750-8971-4933-a7eb-d61168d4710c
# ╟─a9f4d36b-7183-4cf3-afef-96e18bf4c8f3
