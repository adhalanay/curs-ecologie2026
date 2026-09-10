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

# ╔═╡ 96931d62-a744-11f0-296e-dd13032fbfff
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using DifferentialEquations
	TableOfContents()
end

# ╔═╡ 34cc7eec-c950-4491-bc7f-8ec29ffef73e
begin
	using QuadGK
	t2 = range(0,10,length=150)
	V(t)=2*sqrt(1000-t^3)
	y2 = [V(t) for t ∈ t2]
	integral, error = quadgk(t -> V(t), 0, 10)
	print(integral)
	fig2=Figure()
	ax2=Axis(fig2[1,1])
	lines!(t2,y2)
	fig2
end

# ╔═╡ 34ee0b2b-12f1-4d5e-9739-be3e8f9de99b
html"""
<style>
	@media screen {
		main {
			margin: 0 auto;
			max-width: 2000px;
    		padding-left: max(283px, 10%);
    		padding-right: max(383px, 10%); 
            # 383px to accomodate TableOfContents(aside=true)
		}
	}
</style>
"""

# ╔═╡ 4ed380d4-314a-426f-9ae3-7597edc47808
md"""
# Cursul III-Integrala și Ecuații Diferențiale

## Integrala

După cum am văzut în cursul precedent, este destul de ușor să calculăm derivata unei funcții pornind de la cîteva formule și aplicînd cîteva reguli de calcul. Dar ce putem inversa procesul? Adică să de determinăm $f(x)$ dacă știm doar $f'(x)$. Procedura care ne permite să facem aceasta se numește **integrare**.

Pentru funcții simple este nu este foarte greu: de exemplu dacă $f'(x)=2x$ atunci orice funcție de forma $f(x)=x^2+c$ are derivata $2x$ pentru orice constantă $c$. De asemenea dacă $f'(x)=e^x$, atunci $f(x)=e^x+c$, pentru orice constantă $c$. În general însă, este mult mai dificil să calculăm **primitiva** decît derivata.

Se pune întrebarea de ce este nevoie să calculăm primitiva unei funcții. Adesea cunoaștem derivata funcției fără să știm funcția însăși: să presupunem că avem acces la înregistrările vitezometrului unei mașini și vrem să determinăm distanța parcursă în fiecare moment $t$. Matematic aceasta înseamnă că știm $x'(t)$ (viteza). Într-o altă situație știm viteza cu care un pacient primește medicament intra-venos (via un debitmetru atașat perfuziei), dar ne interesează cantitatea administrată. Dacă putem găsi o expresie analitică pentru viteză sau debit atunci putem încerca să calculăm primitiva (deplasarea sau cantitatea de medicament în cazul al doilea). De cele mai multe ori însă avem doar un grafic de forma:
"""

# ╔═╡ 95bf4b62-d4c3-4714-bcb2-c108cfafa53b
begin
	f(x) = x-cos(x)*sin(x)
	t1 = range(0,10,length=15)
	y1 = [f(t) for t ∈ t1]
	fig1=Figure()
	ax=Axis(fig1[1,1])
	scatter!(t1,y1)
	lines!(t1,y1)
	fig1
end

# ╔═╡ 2a0ac124-34d7-455c-bccc-1f875083981a
md"""
Atunci avem nevoie să calculăm folosind metode de aproximare. Ne amintim de metoda Euler de calcul aproximativ al derivatei
```math
X'(t)\approx \frac{X(t+\Delta t)-X(t)}{\Delta t}.
```
Atunci avem că $X(t+\Delta t) \approx X(t)+\Delta t\cdot X'(t)$. Acum repetăm construcția pentru $t=0$,$t=\Delta t$, $t=2\Delta t$, etc. Avem acum
```math
\begin{aligned}
X(\Delta t)=X(0) + \Delta t\cdot X'(0) \\
X(2\Delta t) = X(\Delta t)+\Delta t\cdot X'(\Delta t)=X(0)+\Delta t\cdot X'(0)+\Delta t\cdot X'(\Delta t) \\
X(3\Delta t)=X(0)+X'(0)\Delta t+X'(\Delta t)\Delta t+X'(2\Delta t)\Delta t
\end{aligned}
```
Putem reface calculul pentru distanța parcursă în funcție de viteză: 
```math
X(t)\approx X(0)+V(0)\Delta t+V(\Delta t)\Delta t+V(2\Delta t)\Delta t+\dots V(n\Delta t)\Delta t
```
cu $n$ numărul de pași pe care îi facem.  Pe scurt scriem
```math
X(t) \approx X(0)+\sum_{k=0}^n X'(k\Delta t)\cdot\Delta t.
```
Ca să obținem egalitate trecem la limită cu $\Delta t \to 0$. Obținem 
```math
X(t) = X(0)+\lim_{\Delta t \to 0} \sum_{k=0}^n X'(k\Delta t)\cdot\Delta t=:X(0)+\int_0^t X'(t)dt.
```
Aceast formulă se numește teorema fundamentală a analizei sau formula Leibniz-Newton. 

Exemplu: presupunem că viteza unei mașini este $V(t)=2\sqrt{1000-t^3}.$ Mașina pornește cu aproximativ $63 m/s=234 km/h$ și se oprește după $10$s. Primitiva funcției $V$ nu poate fi exprimată prin funcții elementare, deci va trebui să calculăm sume Riemann. 
"""

# ╔═╡ 5300c0d9-36b7-421f-ba05-af81b6c786ed
md"""
Presupunem că $X(0)=0$ și $\Delta t=0.1$. Avem următoarele calcule:
"""

# ╔═╡ c6560f36-88ef-4160-b63a-12a96a0d726c
begin
Δt=0.1
X=0
Xs=[0.0]
print("---------------------------------------------------\n")
print("t               |     V(t)              | X(t+1)         \n")
print("---------------------------------------------------\n")
for k ∈ range(0,100,length=100) 
    X = X+V(k*Δt)*Δt
    append!(Xs,X)
    print(k*Δt,'|', V(k*Δt),'|', X,'\n')
    end
print(integral)
end

# ╔═╡ 4124bec3-78dd-4f8c-b3a6-ab7046bdbd19
begin
	t3 = range(0,10,length=101)
    fig3=Figure()
    ax3=Axis(fig3[1,1],title="Graficul deplasării")
    lines!(t3,Xs)
    fig3
end

# ╔═╡ f1728f88-5926-4ccb-b1b2-0424d0c3027e
md"""
Ce am făcut de fapt pentru a calcula deplasarea:
 - am descompus intervalul $(0,t)$ în segmente de lungime egală $\Delta t$;
 - am presupus că $V(t)$ este constantă pe fiecare sub interval;
 - am calculat distanța parcursă pe fiecare interval cu formula distanța=viteza $\times$ timpul;
 - am adunat lungimile obținute.

O observație simplă este că produsul $V(k\cdot\Delta t)\cdot\Delta t$ este aria dreptunghiului cu laturile respective. Deci suma Riemann aproximează aria aflată sub graficul funcției $V'$, iar integrala este egală exact cu aria de sub grafic.

Dacă vrem să calculăm aria unei porțiuni de sub grafic între două valori ale lui $t$, să spunem $a < b$. Atunci aria este
```math
\int_a^b V(t)dt=\int_0^b V(t)dt-\int_0^a V(t)dt.
```

Dacă $F(t)$ este primitiva lui $f(t)$, atunci
```math
\int_a^b f(t)dt=F(b)-F(a).
```

### Formula de medie
Valoarea medie a unei funcții $f:[a,b] \to \mathbb{R}$ este 

```math
f(c) = \frac{1}{b-a} \int_a^b f(t) dt.
```

## Ecuații diferențiale

În cîteva cazuri fericite putem folosi primitiva pentru a găsi soluții explicite pentru ecuații diferențiale. În general ecuațiile sînt de forma $x'=f(x)$ sau echivalent  $\frac{dx}{dt}(t)=f(x(t))$. Cele mai simple cazuri sînt ecuațiile
$x'=kx$ și $x'=-kx$ cu $k > 0$. 

Presupunem că o populație $X$ are o rată de nașteri $b$, o rată de deces $d$, o rată de imigrare $i$ și o rată de emigrare $e$, toate presupuse constante. Atunci variația populației este 
```math
X'=rX,
```
unde $r=b+i-d-e$, pe care îl presupunem pozitivă. Care este comportarea lui $X$? Putem aplica metoda lui Euler, dar în cazul nostru, noi am văzut că
```math
\frac{de^{kt}}{dt}=ke^{kt}.
```
Deci funcția $X(t)=e^{kt}$ este soluția ecuației. Creșterea dată de această regulă se numește **exponențială**. Soluția generală a ecuației este $$X(t)=X(0)e^{rt}.$$

Să presupunem acum că $r<0$, spre exemplu dacă populația este în scădere. Atunci ecuația noastră va fi $X'=-kX$ cu soluția
$X(t)=X(0)\cdot e^{-kt}.$

Ecuațiile de acest tip modelează fenomene în care o fracție constantă $k$ este extrasă în fiecare moment. Spre exemplu o substanță care se descompune sau
o populație care scade cu o fracție constantă etc. Graficul soluției este
"""

# ╔═╡ ae8a653a-48fa-4247-859b-d2f4dd195966
begin
	f1(t) = 0.3*exp(-0.5*t)
    y4=[f1(t) for t ∈ t2]
	fig4=Figure()
	ax4=Axis(fig4[1,1])
	lines!(t2,y4)
	fig4
end

# ╔═╡ 47bff75c-3d9a-4a71-b596-37155dfbb2a3
md"""
Curba aceasta se numește scădere exponențială. 

Dacă $r >0$ atunci avem creștere exponențială.
"""

# ╔═╡ 18b60c22-1365-4e86-90c5-466eef2c856d
begin
	f2(t) = 0.3*exp(0.5*t)
	y5=[f2(t) for t ∈ t2]
	fig5=Figure()
	ax5=Axis(fig5[1,1])
	lines!(t2,y5)
	fig5
end

# ╔═╡ 894ad5c0-b051-45c7-b465-98a14659dd9f
md"""
Modelarea cu ajutorul unei ecuații diferențiale de tip exponențial este valabilă doar pe un termen foarte scurt, deoarece prezice că populația nu se oprește din crescut, ceea ce este evident imposibil. Deci modelul nostru trebuie modificat pentru a lua în considerare fenomenul de limitare a creșterii.

Același lucru se întîmplă și pentru scăderea exponențală. Foarte repede ajungem la valori foarte mici, în practică $0$.

Atunci cînd formulăm modele matematice, trebuie să verificăm întotdeauna dacă rezultatele obținute se potrivesc cu observațiile concrete și cu principiile biologice.

## Sisteme dinamice. Puncte de echilibru

Modelele deterministe care studiază evoluția în timp a unor variabile de stare se numesc sisteme dinamice. Dacă presupunem că timpul ia doar valori discrete (numere întregi) atunci avem de-a face cu **sisteme discrete**. Dacă lăsăm timpul să ia toate valorile dintr-un interval atunci avem **sisteme continue**. Principala problemă cu care ne confruntăm este că sînt cunoscute doar legile după care variază variabilele și nu formule după care evoluează variabilele. În general nu putem scrie explicit orbitele unui astfel de sistem, așa că recurgem la **teoria calitativă** a acestora. Un prim pas este studiul **punctelor de echilibru.**

Pentru o ecuație diferențială $X'=f(X)$, punctele de echilibru sînt punctele $X_0$ pentru care $f(X_0)=0$, adică în aceste puncte sistemul are "viteza" nulă

Pînă acum am studiat ecuațiile diferențiale folosind metoda Euler pentru simulări numerice, dar acest lucru are mari limitări. Fie ecuația 
```math
x'=r\left(1-\frac{x}{k}\right)\left(\frac{x}{a}-1\right)
```
foarte asemănătoare cu ecuația logistică. Să simulăm comportamentul ei pentru cîteva valori ale parametrilor și ale valorii inițiale $x(0).$ Anume luăm $r=0.1,a=5,k=100$ și valorile inițiale $10,20,150$, respectiv $4$.
"""

# ╔═╡ 608f4fa5-37e6-465e-886b-5412fa1b3e5c
begin

r=0.1;a=5;k=100
f(u,p,t)=r*(1-u/k)*(u/a-1)
tspan = (0.0,800.0)
u0=10.0
prob1 = ODEProblem(f,10.0,tspan)
sol1 = solve(prob1, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
prob2 = ODEProblem(f,20.0,tspan)
sol2 = solve(prob2, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
prob3 = ODEProblem(f,150.0,tspan)
sol3 = solve(prob3, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
#prob4 = ODEProblem(f,4.0,tspan)
#sol4 = solve(prob4, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
fig6=Figure()
ax6=Axis(fig6[1,1])
series!(sol1.t,hcat(sol1.u...),linewidth=5)
series!(sol2.t,hcat(sol2.u...),color=:reds,linewidth=5)
series!(sol3.t,hcat(sol3.u...),color=:vik,linewidth=5)
# WGLMakie.series!(sol4.t,hcat(sol4.u...))
fig6
end

# ╔═╡ df77753d-59db-4378-97ad-33551d0973de
begin
fig7=Figure()
ax7 = Axis(fig7[1,1])
prob4 = ODEProblem(f,4.0,tspan)
sol4 = solve(prob4, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
series!(sol4.t,hcat(sol4.u...),color= :roma, linewidth=5)
fig7
end

# ╔═╡ e25589cd-579f-47f5-be8c-7573e76cc2de
@bind u1 html"<input type='range' style='width: 1000px' min='0' max='200' step='0.05'>"

# ╔═╡ 8217cecf-19c9-43c1-93b4-940c1cdd15d0
u1

# ╔═╡ 916800aa-7986-48ad-af69-934e115e5972
begin


prob_var = ODEProblem(f,u1,tspan)
sol_var = solve(prob_var, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
fig_var=Figure()
ax_var=Axis(fig_var[1,1])
series!(sol_var.t,hcat(sol_var.u...),linewidth=5)
fig_var
end

# ╔═╡ f171e206-6e46-4525-9ecb-cf4b3a982b63
md"""
Ca să vedem mai bine de ce se întîmplă vom mai lua două valori inițiale, pentru care $f=0$, adică $u0=k$ și $u0=a$.
"""

# ╔═╡ b0eacf7f-e6a9-4f98-9287-a58d87881397
begin

prob5 = ODEProblem(f,a,tspan)
sol5 = solve(prob5, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
prob6 = ODEProblem(f,k,tspan)
sol6 = solve(prob6, Tsit5(), reltol = 1e-8, abstol = 1e-8,saveat=0.01)
fig8=Figure()
ax8=Axis(fig8[1,1])
series!(sol5.t,hcat(sol5.u...),linewidth=5)
series!(sol6.t,hcat(sol6.u...),color=:reds,linewidth=5)
fig8
end

# ╔═╡ 7a5c562f-8c8e-41a1-88b2-dfaf2153bbdc
md"""
Se vede că în primele trei cazuri, comportarea este asemănătoare cu cea a ecuației logistice: dacă pornim cu o valoare mai mică decît $k$, atunci soluția va crește către el, iar dacă pornim cu o valoare mai mare, soluția va scădea către $k$. În cazul al patrulea se întimplă ceva diferit, soluția nu crește cum ne-am aștepta ci scade. 

Comportamentul nou provine din faptul că am luat $x(0)$ mai mică decît $a$ în acest caz, în vreme ce în primele trei cazuri ele au fost toate mai mari decît $a$.

Cum nu putem simula pentru toate valorile inițiale este important să găsim o metodă să determinăm comportamentul soluțiilor. în cazul dimensiunii $1$, o astfel de metodă există. Pentru cazul multi-dimensional această metodă oferă informații importante, chiar dacă nu complete. 

Primul pas constă în a determina punctele în care sistemul nu se schimbă, numite **punctele de echilibru**. Aceast sînt puncte pentru care $x'=0$ (deoarece am văzut că derivata unei constante este $0$, și reciproc doar funcțiile constante au această proprietate). Practic pentru ecuația $x'=f(x)$, trebuie să găsim valorile lui $x$ pentru care $f(x)=0$. De exemplu pentru ecuația noastră punctele de echilibru satisfac
```math
\begin{aligned}
1-\frac{x}{k}=0 & \Longrightarrow & x=k \text{ și } \\
1-\frac{x}{a}=0 & \Longrightarrow & x=a. 
\end{aligned}
```

## Stabilitatea echilibrelor - Metoda punctului intermediar

După ce am determinat punctele de echilibru este important să decidem ce se întimplă în jurul acestui punct: dacă pornim din apropierea echilibrului răminem în apropiere sau ne îndepărtăm? De asemenea soluțiile care pornesc dintr-un punct oarecare se apropie de una dintre soluțiile constante? 

Matematic vorbim de **stabilitatea** echilibrelor. Un echilbru este stabil dacă pornind din apropierea lui rămînem în apropiere și este instabil dacă pornind din apropierea echilibrului ne îndepărtăm.

Pentru a studia stabilitatea unui echilibru trebuie să vedem semnul funcției $f$. Să luăm ca exemplu ecuația logistică 
```math
x'=rx\left(1-\frac{x}{k}\right).
```
Punctele sale de echibru sînt $x_0=0$ și $x_1=k>0$. Pentru valori ale lui $x$ între $0$ și $k$ $f(x) > 0,$ iar pentru $x > k,$
$f(x) < 0$. Din desen devine imediată stabilitatea: dacă vectorii îndepărtează punctul de punctul de echilibru, atunci acesta este instabil, iar dacă îl apropie este stabil.


"""

# ╔═╡ 6db1dbae-20f2-44c5-9a91-4ee14cdf744f
begin
fig9=Figure()
ax9=Axis(fig9[1,1])
hlines!(0;xmin=0.0,xmax=10.0)
p1=Point2f(0.0,0.0)
p2=Point2f(10.0,0.0)
p3=Point2f(5.0,0.0)
p4=Point2f(14.0,0.0)
p5=Point2f(-4.0,0.0)
scatter!(p1)
text!(p1,text="x=0",align = (:center, :bottom),offset=(0,10))
arrows!([p1],[p3],linewidth=3)
scatter!(p2)
text!(p2,text="x=k",align = (:center, :bottom),offset=(0,10))
arrows!([p4],[p5],linewidth=3)
xlims!(ax9,-0.5,20.0)
ylims!(ax9,-0.5,10.0)
fig9
end

# ╔═╡ 3ea69d1b-21ee-4da3-b10f-2cd5e4e858cd
md"""
Nu este întotdeauna la fel de ușor ca în cazul ecuației logistice să determinăm semnul funcției care definește ecuația diferențială. Să luăm următoarea ecuație, bazată pe modelul logistic, unde o parte din populație este scoasă spre exemplu prin vînătoare sau pescuit:
```math
x'=0,2 x\left(1-\frac{x}{1000}\right)-0,1x.
```
Se poate cu ușurință vedea că cele două echilibre sînt $x=0$ și $x=500$. Pentru a găsi semnul lui $x'$ pe fiecare interval $(0,500)$, respectiv $(500,\infty)$ este suficient să luăm cîte un punct în fiecare interval și să vedem semnul funcției în acest punct.  Astfel avem o figură similară cu cea precedentă. 

Am folosit implicit faptul că o funcție continuă nu poate schimba semnul fără să treacă prin . Aceasta este teorema valorii intermediare. Noi cunoaștem toate zerourile, funcția nu poate schimba semnul fără să treacă prin echilibru.

## Stabilitatea echilibrelor: analiza liniară a stabilității

Metoda grafică are succes în dimensiune $1$, dar în dimensiune mai mare nu funcționează așa de bine. Studiem dinnou ecuația logistică 
```math
x'=x\left(1-\frac{x}{k}\right).
```
Facem graficul funcției $f(x)=x\left(1-\frac{x}{k}\right)$ și vedem că punctele de echilibru sînt cele în care graficul taie axa $Ox$. Pentru $x=0$ funcția trece de la negativ la pozitiv, deci crește, iar pentru $x=k$ din contră trece de la pozitiv la negativ, decit scade. Prin urmare panta tangentei la grafic în primul echilibru este pozitivă, iar în cel de-al doilea panta este negativă.
"""

# ╔═╡ 9654f5b8-c980-48ca-8dc9-99da68f331a9
begin
k1=5
f3(x)=x*(1-x/k1)
ts = range(0.0,20.0,length=100)
ys = [f3(t) for t ∈ ts]
fig10=Figure()
ax10=Axis(fig10[1,1])
hlines!(0;xmin=0.0,xmax=10.0)
p11=Point2f(0.0,0.0)
p21=Point2f(5.0,0.0)
p31=Point2f(1.0,0.0)
p41=Point2f(7.0,0.0)
p51=Point2f(-1.0,0.0)
scatter!(p11)
text!(p11,text="x=0",align = (:center, :bottom),offset=(0,10))
arrows!([p11],[p31],linewidth=3)
scatter!(p21)
text!(p21,text="x=k",align = (:center, :bottom),offset=(0,10))
arrows!([p41],[p51],linewidth=3)
lines!(ts,ys)
xlims!(ax10,-0.5,8.0)
ylims!(ax10,-3.0,3.0)
fig10
end

# ╔═╡ a5291a15-baac-4a5c-9334-e90dcd8468e7
md"""
Calculăm derivata funcției:
```math
f'(x)=1-\frac{2x}{k}.
```
Deci $f'(0)=1 > 0$ și $f'(k)=-1 < 0$. Semnul derivatei ne indică stabilitatea unui echilibru $x^*$:
 - dacă $f'(x^*) > 0,$ atunci echilibrul este instabil;
 - dacă $f'(x^*) < 0,$ atunci echilibrul este stabil. 
 
 Aceasta este o primă aplicației a **teoremei Hartman-Grobman**: comportarea unei ecuații diferențiale în jurul unui punct de echilibru este determinată de aproximarea sa liniară. 

## Efectul Allee

În anumite populații s-a constat următorul fenomen destul de neintuitiv: pentru anumite populații supraviețuirea nu are loc decît pentru un anumit număr minim de indivizi. Dacă numărul de indivizi scade sub un anumit prag, atunci populația intră în declin, putînd chiar dispărea. Acesta este *efectul Alle*.

Efectul Allee se modelează modificînd ecuația logistică într-o formă pe care am mai întilnit-o:
```math
x'=rx\left(1-\frac{x}{k}\right)\left(\frac{x}{a}-1\right).
```

Avem trei puncte de echilibru, anume $x=0$, $x=k$ și $x=a$. Pentru a studia stabilitatea cu ajutorul metodei punctului intermediar luăm cîte un punct în intervalele $(0,a)$, $(a,k)$ și $(k,\infty)$ (am presupus $a < k$). 

Privind graficul vedem că echilibrul $x=0$ este stabil, $x=a$ este instabil și $x=k$ este iar stabil. Putem verifica acest lucru și prin calcule. Avem
```math
f(x)=-\frac{r}{ak}x^3+\frac{r}{a}x^2+\frac{r}{k}x^2-rx,
```
deci
```math
f'(x)=-\frac{3r}{ak}x^2+\frac{2r}{a}x+\frac{2r}{k}x-r.
```
Evaluăm acum derivata în punctele de echilibru:
```math
\begin{aligned}
\left.\frac{df}{dx}\right|_{x=0} & = & -r & < & 0 \\
\left.\frac{df}{dx}\right|_{x=a} & = & r\left(1-\frac{a}{k}\right) & > & 0 \\
\left.\frac{df}{dx}\right|_{x=k} & = & r\left(1-\frac{k}{a}\right) & < & 0 
\end{aligned}
```
De aici rezultă stabilitatea echilibrelor.
"""

# ╔═╡ c5bd9a17-9455-4a64-a7e6-d7ce4dae29b6
begin
k2=7.0;a2=3.0;r2=5.0
f4(x)=r2*x*(1-x/k2)*(x/a2-1)
ts4 = range(0.0,20.0,length=100)
ys4 = [f4(t) for t ∈ ts4]
fig11=Figure()
ax11=Axis(fig11[1,1])
hlines!(0;xmin=0.0,xmax=10.0)
p12=Point2f(0.0,0.0)
p22=Point2f(k2,0.0)
p32=Point2f(a2,0.0)
p42=Point2f(1.0,0.0)
scatter!(p12)
text!(p12,text="x=0",align = (:center, :bottom),offset=(0,10))
scatter!(p22)
text!(p22,text="x=k",align = (:center, :bottom),offset=(0,10))
scatter!(p32)
text!(p32,text="x=a",align = (:center, :bottom),offset=(0,10))
lines!(ts4,ys4)
xlims!(ax11,-0.5,8.0)
ylims!(ax11,-4.0,7.0)
fig11
end

# ╔═╡ Cell order:
# ╟─96931d62-a744-11f0-296e-dd13032fbfff
# ╟─34ee0b2b-12f1-4d5e-9739-be3e8f9de99b
# ╟─4ed380d4-314a-426f-9ae3-7597edc47808
# ╟─95bf4b62-d4c3-4714-bcb2-c108cfafa53b
# ╟─2a0ac124-34d7-455c-bccc-1f875083981a
# ╟─34cc7eec-c950-4491-bc7f-8ec29ffef73e
# ╟─5300c0d9-36b7-421f-ba05-af81b6c786ed
# ╟─c6560f36-88ef-4160-b63a-12a96a0d726c
# ╟─4124bec3-78dd-4f8c-b3a6-ab7046bdbd19
# ╟─f1728f88-5926-4ccb-b1b2-0424d0c3027e
# ╟─ae8a653a-48fa-4247-859b-d2f4dd195966
# ╟─47bff75c-3d9a-4a71-b596-37155dfbb2a3
# ╟─18b60c22-1365-4e86-90c5-466eef2c856d
# ╟─894ad5c0-b051-45c7-b465-98a14659dd9f
# ╠═608f4fa5-37e6-465e-886b-5412fa1b3e5c
# ╟─df77753d-59db-4378-97ad-33551d0973de
# ╠═e25589cd-579f-47f5-be8c-7573e76cc2de
# ╠═8217cecf-19c9-43c1-93b4-940c1cdd15d0
# ╠═916800aa-7986-48ad-af69-934e115e5972
# ╟─f171e206-6e46-4525-9ecb-cf4b3a982b63
# ╟─b0eacf7f-e6a9-4f98-9287-a58d87881397
# ╟─7a5c562f-8c8e-41a1-88b2-dfaf2153bbdc
# ╟─6db1dbae-20f2-44c5-9a91-4ee14cdf744f
# ╟─3ea69d1b-21ee-4da3-b10f-2cd5e4e858cd
# ╟─9654f5b8-c980-48ca-8dc9-99da68f331a9
# ╟─a5291a15-baac-4a5c-9334-e90dcd8468e7
# ╟─c5bd9a17-9455-4a64-a7e6-d7ce4dae29b6
