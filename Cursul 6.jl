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

# ╔═╡ c571d396-bd45-11f0-2305-732890aee08a
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

# ╔═╡ e2ea1bbb-5698-4d53-a294-477a238a6ff8
WidthOverDocs()

# ╔═╡ 30cad42d-f396-4efd-939c-ad86e605eda1
md"""
# Comportament periodic și bifurcații

Multă vreme s-a crezut că toate sistemele din natura evoluează către o stare de echilibru. Aceasta a făcut ca modelarea să se concentreze pe găsirea de modele (sisteme discrete sau de ecuații diferențiale) care prezintă echilibre stabile. S-a observat însă că foarte multe fenomene au un comportment oscilant permanent. Astfel sînt sistemul reprezentat de interacțiunea între rîși și iepuri sau reacția Belousov-Jabotinski din chimia anorganică. Alte fenomene oscilatorii sînt: evoluția concentrației unui hormon în corp, expresia genelor este de asemenea oscilatorie cu cicluri de ore sau zile. 

În comportamentul unui sistem trebuie să facem distincția între comportamente **tranziente** (pe termen scurt) și cele **asimptotice** (pe termen lung). Pentru a studia comportamentul considerăm noțiunea de **atractor**. 

Dacă ``X`` este spațiul fazelor, o mulțime ``A \subset X`` se numește atractor dacă există o vecinătate de puncte ințiale ale căror orbite se apropie de ``A``. Mai exact ``A`` este atractor dacă există o submulțime ``X_0 \subset X`` astfel că pentru orice
valoare inițială din ``X_0`` avem că 
```math
\text{dist}(X(t),A)) \rightarrow 0.
```

Am văzut deja exemple de atractori, anume punctele de echilibru stabile. O altă clasă de atractori este reprezentată de **oscilațiile stabile**. Oscilațiile pe care le-am întîlnit deja în cazul sistemului rîși-iepuri sau pendulul fără frecare nu sînt stabile, ci depind de valorile ințiale. 

Oscilația poate fi privită din două puncte de vedere: o variabilă de stare ``X(t)`` este periodică dacă există un ``P`` astfel încît ``X(t+P)=X(t)`` pentru orice ``t.`` În spațiul fazelor a traiectorie este o oscilație dacă este o curbă închisă. Ne interesează
doar oscilațiile stabile, adică orbitele periodice care sînt și atractori.Acestea se mai numesc **cicluri limită stabile.**

## Clarinetul lui Raileigh

Un prim exemplu este dat de modelului clarinetului lui Raleigh. Acesta este de fapt o modificare a pendulului cu frecare:
```math
\begin{aligned}
x' & = & v \\
v' & = & - x-(v^3-v)
\end{aligned}
```
Portretul de fază este:
"""

# ╔═╡ 35b51444-672b-43ff-8377-896d33a60c18
begin
	f1(x) = Point2f(x[2],-x[1]-(x[2]^3-x[2]))
	fig1 = Figure()
	ax1 = Axis(fig1[1,1])
	xs1 = LinRange(-2, 2, 10)
    ys1 = LinRange(-2, 2, 10)
	# us1 = [-y for x ∈ xs1, y ∈ ys1]
	# vs1 = [x for x ∈ xs1, y ∈ ys1]
	# strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax1,f1,xs1,ys1,colormap = :magma,linewidth = 1)
	fig1
end

# ╔═╡ b57c6a08-6a72-4542-90ea-4637c13b6da3
md"""
Cîteva orbite se văd mai jos și seriile temporale corespunzătoare sînt reprezentate mai jos. Se observă că orbita roșie este periodică și atrage toate celelalte orbite.
"""

# ╔═╡ d961200f-1dd6-42a4-9bc0-c4cb17373e09
begin
function f_raileigh!(du,u,p,t)
    du[1]=u[2]
    du[2]=-u[1]-(u[2]^3-u[2])
	return nothing
end
t_final=30.0
s_t=0.01
init_vals = [[0.5,0.5],[0.0,0.1],[1.0,0.9]]
probs = [] 
	for u in init_vals 
		raileigh = CoupledODEs(f_raileigh!,u)
		push!(probs,raileigh)
	end
fig30 = Figure(size=(2000,1000))
	ax30 = Axis(fig30[1,1],title="Portretul de fază")
	let i = 2
	for prob in probs
		sol,t = trajectory(prob,t_final;Ttr=2.2, Δt=s_t)
		scatter!(ax30,sol,colormap=:thermal)
		ax = Axis(fig30[1,i],title="Serii temporale")
		for col in columns(sol)
			lines!(ax,t,col,linewidth=4)
		end
		i+=1
	end
	fig30
end
end

# ╔═╡ 9056c9e6-7594-4f3d-97a3-95593061431a
md"""

## Reglarea hormonală

Unul din mecanismele care dau naștere la orbite periodice stabile este feddback-ul negativ. Un model pentru secreția de hormoni este următorul: considerăm secreția unui hormon secretat de o gonadă. Notăm cu ``G`` cantitate de hormon. Secreția sa este controlată de doi hormoni secretați de glanda hipofiză: hormonul luteinizant (LH) și cel foliculostimulant (FSH). Ambii hormoni stimulează gonadele. Pentru simplitate vom considera un singur hormon notat cu ``P``. La rîndul ei hipofiza este controlată de 
creier (mai precis hipotalamus) via hormonul eliberator de gonadotropine ``H``. Cînd cantitatea de hormon sexual este prea mare hipotalamus scade producția de hormon eliberator și deci cantitatea celorlalți doi hormoni scade și ea. Modelul pe care îl obținem este
```math
\begin{aligned}
H' & = & \frac{1}{1+G^n} & -  k_1 H \\
P' & = & H & -  k_2 P \\
G' & = & P & -  k_3 G.
\end{aligned}
```
Chiar dacă vorbim de feedback negativ funcția de feedback este pozitivă. Tot ceea ce face este săs cadă producția de ``H``, deoarece hipotalamusul nu poate absorbi la loc hormonul. O simulare a sistemului este 
"""

# ╔═╡ 1145c760-e750-4c1b-b6c1-6c178ae6ef3c
md""" ``n=`` $(@bind n PlutoUI.Slider(2:10,show_value=true)) """

# ╔═╡ 58b70a0f-8180-49f7-b5b9-2e7147cf6b00
begin
# tspan=(0.0,30.0)
function f!(du,u,p,t)
    du[1]=1/(u[3]^n+1) -p[1]*u[1]
    du[2]=u[1]-p[2]*u[2]
	du[3]=u[2]-p[3]*u[3]
	return nothing
end
p = [0.2,0.2,0.2]
u_init=[0.2,0.2,0.2]
horm = CoupledODEs(f!,u_init,p)
total_time = 400.0
sampling_time = 0.02
Y, t = trajectory(horm, total_time; Ttr = 2.2, Δt = sampling_time)
fig3 = Figure(size=(1700,500))
ax3 = Axis(fig3[1,1])
	for col in columns(Y)
		lines!(t,col)
	end
ax33=Axis3(fig3[1,2])
lines!(ax33,Y[:,1],Y[:,2],Y[:,3])
fig3
end

# ╔═╡ 6fad810e-9c3b-40dd-851b-145ddbc8ad23
md"""
Se vede că pentru ``n < 8`` sistemul are un echilibru stabil, apoi cînd ``n`` trece de 8 acesta devine instabil și apare un ciclu stabil. Intuitiv putem spune că un feedbak prea puternic duce la apariția unei orbite periodice stabile, dar dacă studiem sistemul
din care am eliminat ecuația din mijloc vom pierde periodicitatea. Deci mai există un mecanism implicat, anume întîrzierea determinată de acțiunea lui ``P``. 

## Expresii genetice oscilante

Factorul de trasncripție ``Hes1`` are un comportament oscilant cu perioada de aproximativ 2 ore. Modelul este următorul: ARN-ul mesager al lui Hes1 (Y) este convertit în proteina ``Hes1`` (X) cu rata ``B``. Se presupune că există un factor de degradare ``Z`` care combinat cu proteina ``Hes1``. Proteina ``Hes`` inhibă propria sa transcripție ceea ce se modelează printr-o funcție sigmoid . De asemenea proteina inhibă și factorul de interacțiune. Sistemul de ecuații este 
```math
\begin{aligned}
X' & = & -AXZ & + & BY & - & CX \\
Y' & = & \frac{E}{1+X^2} & - & DY \\
Z' & = & -AXZ & + & \frac{F}{1+X^2} & - & GZ
\end{aligned}
```
O altă cale spre apariția orbitelor periodice este dată de anumite bifurcații pe care le vom discuta mai jos.
"""

# ╔═╡ acea996d-1913-4488-93fb-f40e880ea59d
begin
# tspan=(0.0,30.0)
function hes1!(du,u,p,t)
	A,B,C,D,E,F,G = p
    du[1] = -A*u[1]*u[3]+B*u[2]-C*u[1]
    du[2] = E/(1+u[1]^2)-D*u[2] 
	du[3] = -A*u[1]*u[3]+F/(1+u[1]^2)-G*u[3]
	return nothing
end
p_hes = [0.022,0.3,0.031,0.028,0.5,20.0,0.3]
u0=[0.2,0.2,0.2]
hes = CoupledODEs(hes1!,u0,p_hes)
total_time0 = 3000.0
# sampling_time = 0.02
Y_hes, t_hes = trajectory(hes, total_time0; Ttr = 2.2, Δt = sampling_time)
fig4 = Figure(size=(1500,1000))
ax4 = Axis(fig4[1,1])
	for col in columns(Y_hes)
		lines!(t_hes,col)
	end
ax41 = Axis3(fig4[1,2])
	lines!(ax41,Y_hes[:,1],Y_hes[:,2],Y_hes[:,3])
fig4
end

# ╔═╡ dd8ad90a-b58e-4a6b-a59f-e9deb5e46e0b
md"""

# Bifurcații 

Am văzut deja că anumite sisteme dinamice au comportamente diferite în funcție de valorile unui parametru. Un astfel de sistem este sistemul logistic.

## Bifurcații ale punctelor de echilibru

Atunci cînd stabilitatea echilibrelor se schimbă o dată cu parametrii spunem că avem o **bifurcație locală.** Spre exemplu în cazul sistemului căprioare-elani scăderea competiției interspecifice duce la stabilizarea la coexistență. Deci studiul bifucațiilor ne
indică ce și cum să schimbăm parametrii pentru a obține rezultatul dorit.

### Bifurcația transcritică 

Să vedem evoluția unei populații cu efectul Allee
```math
x' = 0.1 x \left(1-\frac{x}{k}\right)\left(\frac{x}{a}-1\right),
```
unde ``k`` este capacitatea de suport a mediului și ``a`` este mărimea minimă a populației care asigură creșterea. Echilibrele sînt ``0,k`` și ``a``. Cînd ``a < k`` atunci ``k`` este stabil și ``a`` este instabil. Invers cînd ``a>k`` ``a`` este stabil și ``k`` este instabil. În acest caz cele două puncte de echilibru se ciocnesc și schimbă stabilitatea între ele. 
"""

# ╔═╡ 43e87692-0176-458c-aae8-afe94db20b42
md""" a=$(@bind a PlutoUI.Slider(300:1500,show_value=true))"""

# ╔═╡ ce3937dc-d5f2-4b8c-a24a-b66cf33dba5b
begin
tspan=(0.0,300.0)
function alle!(du,u,p,t)
	a,k = p
    du[1] = 0.1*u[1]*(1-u[1]/k)*(u[1]/a-1)
end
k0=800.0
p_alle = [a,k0]
init_values=[[200.0],[400.0],[600.0],[800.0],[1000.0],[1200.0],[1500.0],[a]]
#total_time0 = 3000.0
# sampling_time = 0.02
	probs1 = [] 
	for u in init_values 
		alle_pr = ODEProblem(alle!,u,tspan,p_alle)
		sol = solve(alle_pr, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
		push!(probs1,sol)
	end
fig5 = Figure(size=(1000,500))
	ax5 = Axis(fig5[1,1],title="Efectul Allee cu k=$k0")
	for sol in probs1
		series!(sol.t,hcat(sol.u...),linewidth=5)
	end
fig5
end


# ╔═╡ 50dd9248-8970-46c1-9f70-5fa6edbf0c3c
md"""
### Bifurcația șa-nod

Pentru lac operan avem ecuația
```math
x' = \frac{a+x^2}{1+x^2}-rx.
```

"""

# ╔═╡ a4268686-3a7e-46e0-8b5b-d41a2dd2e878
md""" r=$(@bind r PlutoUI.Slider(4.0:-0.01:0.01, show_value=true))"""

# ╔═╡ 7d4f375e-1603-484a-9d9c-6d321d54f0e6
begin
	using Roots
function lac!(du,u,p,t)
	a,r = p
    du[1] = (a+u[1]^2)/(1+u[1]^2)-r*u[1]
end
a0=0.01
lac0(x)=(a0+x^2)-r*x*(1+x^2)
pars_lac = [a0,r]
init_values_lac=[[0.2],[0.4],[0.6],[0.8],[1.0],[1.2],[1.5]]
	xs = LinRange(0.0,2.0,1000)
	ys = [lac0(x) for x in xs]
	zs = [0 for _ in xs]
#total_time0 = 3000.0
# sampling_time = 0.02
	probs2 = [] 
	for u in init_values_lac 
		prob = ODEProblem(lac!,u,tspan,pars_lac)
		sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
		push!(probs2,sol)
	end
fig6 = Figure(size=(1000,500))
	ax6 = Axis(fig6[1,1],title="Lac operon cu a=$a0")
	for sol in probs2
		series!(ax6, sol.t,hcat(sol.u...),linewidth=5)
	end
	ax61 = Axis(fig6[1,2],title="Punctele de echilibru pentru r=$r")
	    scatter!(ax61, xs,ys)
	    scatter!(ax61,xs,zs)
fig6
	
end

# ╔═╡ 190fe814-24fe-4e0e-a98c-e50de4532cfb
md"""
Se vede că pentru valori mari ale lui ``r`` avem un singur echilibru care este stabil. Apoi pe măsură ce ``r`` scade apar încă două echilibre, unul stabil și unul instabil. Acesta se apropie din ce în ce mai mult de primul echilibrul stabil, pînă se ciocnesc
și dispar, rămînînd doar noul echilibru stabil. O bifurcație unde apare o pereche de echilibre se numește bifurcație **șa-nod.**

Studiul bifurcaților este foarte important pentru că ne oferă o imagine calitativă despre comportamentul sistemelor și dependența acestuia de valorile parametrilor. 

## Bifurcația Hopf

În celulele de drojdie s-au observat pentru prima oară oscilațiile glicolitice. Reacția de glicoliză este guvernată de enzima fosfofructochinază(PFK). La încput glucoza este convertită fructoză-6-fosofat (F6F) apoi PFK convertește F6F în fructoză 1,6-bifosfat (FBP). FBP intră în metabolismul celular și produce ATP (adenozin trifosfat). Reacția PFK are nevoie de o moleculă de ATP și produce două molecule de ADP, mai puțin utilă. Avem reacțiile 
```math
\begin{aligned}
\emptyset \xrightarrow{V_0} F6P \\
F6P + 2 ADP \xrightarrow[PFK]{c} 3 ADP \\
ADP \xrightarrow{k} \emptyset
\end{aligned}
```

sistemul de ecuații corespunzător este 
```math
\begin{aligned}
S' & = & V_0 & - & cSP^2 \\
P' & = & cSP^2 & - & kP,
\end{aligned}
```
unde ``S=[F6P]`` și ``P=[ADP].``
""" 


# ╔═╡ 8533a47c-964d-41eb-adb4-db15b8025641
md""" c = $(@bind c PlutoUI.Slider(1.20:-0.01:0.80,show_value=true)) """

# ╔═╡ 3da8634d-9bf4-4dbc-8494-1e7e0974d8f0
begin
# tspan=(0.0,30.0)
function glico!(du,u,p,t)
	c,k = p
    du[1] = 1-c*u[1]*u[2]^2
    du[2] = c*u[1]*u[2]^2-k*u[2]
	return nothing
end
p_glico = [c, 1.0]
init_vals0 = [1.2, 1.2]
glic = CoupledODEs(glico!,init_vals0,p_glico)
total_time1 = 350.0
# sampling_time = 0.02
Y_glico, t_glico = trajectory(glic, total_time1; Ttr = 2.2, Δt = sampling_time)
fig7 = Figure(size=(1500,1000))
ax7 = Axis(fig7[1,1])
		scatter!(ax7,Y_glico,markersize=3)
ax71 = Axis(fig7[1,2])
	lines!(ax71,t_glico,Y_glico[:,1])
	lines!(ax71,t_glico,Y_glico[:,2])
fig7
end

# ╔═╡ 575d6d45-b5cd-44be-a09b-809769bc7d60
md"""
Se vede că atunci cînd parametrul ``c`` scade echilibrul stabil dispare și este înlocuit cu un ciclu limită stabil. Acest tip de bifurcație se numește bifurcația Poincare-Andronov-Hopf. 

Un alt model apare în ecologie. Acesta este un model de pradă-prădător mai realist decît Lotka-Volterra, numit modelul Holling-Tanner. Populația de pradă va fi ``N``, iar prădătorii ``P``. Presupunem că populația prăzii crește logistic. Pentru prădătorii
hrănirea este ``f(N)P``. De obicei
```math
f(N)=\frac{C_{max}N}{N=h},
```
unde ``C_[max]`` este capacitatea maximă de absorbție, iar ``h`` este densitatea de pradă pentru care consumul este la jumătate de rata maximă. Sistemul de ecuații va fi
```math
\begin{aligned}
N' & = & r_1N\left(1-\frac{N}{k}\right) & - &\frac{wN}{N+d}P \\
P' & = & r_2P\left(1-\frac{jP}{N}\right)
\end{aligned}
```

"""

# ╔═╡ abf156c0-412b-4263-9233-a08cf05c6a56
md""" w=$(@bind w PlutoUI.Slider(0.3:0.01:1.0,show_value=true))"""

# ╔═╡ aeed174d-f167-45be-b6a8-60892b036052
begin
# tspan=(0.0,30.0)
function holling!(du,u,p,t)
	r1,r2,k,d,j,w = p
    du[1] = r1*u[1]*(1-u[1]/k)-w*u[1]*u[2]/(u[1]+d)
    du[2] = r2*u[2]*(1-j*u[2]/u[1])
	return nothing
end
p_holling = [1.0,0.1,7.0,1.0,1.0,w]
init_vals1 = [1.2, 1.2]
holling = CoupledODEs(holling!,init_vals1,p_holling)
total_time2 = 550.0
# sampling_time = 0.02
Y_holling, t_holling = trajectory(holling, total_time2; Ttr = 2.2, Δt = sampling_time)
fig8 = Figure(size=(1500,1000))
ax8 = Axis(fig8[1,1])
		scatter!(ax8,Y_holling,markersize=3)
ax81 = Axis(fig8[1,2])
	lines!(ax81,t_holling,Y_holling[:,1])
	lines!(ax81,t_holling,Y_holling[:,2])
fig8
end

# ╔═╡ dc19464e-292c-4964-976e-8353b746b1eb
md"""
Se observă că avem de asemenea o bifurcație Hopf în funcție de parametrul ``w=C_{max}.``
"""

# ╔═╡ Cell order:
# ╟─c571d396-bd45-11f0-2305-732890aee08a
# ╟─e2ea1bbb-5698-4d53-a294-477a238a6ff8
# ╟─30cad42d-f396-4efd-939c-ad86e605eda1
# ╟─35b51444-672b-43ff-8377-896d33a60c18
# ╟─b57c6a08-6a72-4542-90ea-4637c13b6da3
# ╠═d961200f-1dd6-42a4-9bc0-c4cb17373e09
# ╟─9056c9e6-7594-4f3d-97a3-95593061431a
# ╟─1145c760-e750-4c1b-b6c1-6c178ae6ef3c
# ╟─58b70a0f-8180-49f7-b5b9-2e7147cf6b00
# ╟─6fad810e-9c3b-40dd-851b-145ddbc8ad23
# ╟─acea996d-1913-4488-93fb-f40e880ea59d
# ╟─dd8ad90a-b58e-4a6b-a59f-e9deb5e46e0b
# ╟─43e87692-0176-458c-aae8-afe94db20b42
# ╟─ce3937dc-d5f2-4b8c-a24a-b66cf33dba5b
# ╟─50dd9248-8970-46c1-9f70-5fa6edbf0c3c
# ╟─a4268686-3a7e-46e0-8b5b-d41a2dd2e878
# ╟─7d4f375e-1603-484a-9d9c-6d321d54f0e6
# ╟─190fe814-24fe-4e0e-a98c-e50de4532cfb
# ╟─8533a47c-964d-41eb-adb4-db15b8025641
# ╟─3da8634d-9bf4-4dbc-8494-1e7e0974d8f0
# ╟─575d6d45-b5cd-44be-a09b-809769bc7d60
# ╟─abf156c0-412b-4263-9233-a08cf05c6a56
# ╟─aeed174d-f167-45be-b6a8-60892b036052
# ╟─dc19464e-292c-4964-976e-8353b746b1eb
