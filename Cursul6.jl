### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> title = "Cursul VI"
#> description = "Comportament periodic și bifurcații"
#> 
#>     [[frontmatter.author]]
#>     name = "Andrei Halanay"

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

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000000
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie, WGLMakie
	using PlutoUI
	using PlutoTeachingTools
	using Roots
	using Printf

	set_theme!(theme_minimal())
	TableOfContents()
end

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ 329b5e95-328d-41d9-8d69-c41a0fc4b4b7
begin
	# un pas al metodei Runge–Kutta de ordinul 4 pentru u' = F(u) (Cursul 5);
	# u poate fi un număr sau un tuplu cu oricîte componente
	function rk4(F, u, h)
		k1 = F(u)
		k2 = F(u .+ h / 2 .* k1)
		k3 = F(u .+ h / 2 .* k2)
		k4 = F(u .+ h .* k3)
		u .+ h / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
	end

	# soluția sistemului u' = F(u) pe [0, T]: vectorul timpilor și cîte un vector pentru
	# fiecare componentă, de exemplu  t, (X, Y) = solutie(F, (1, 0), 20)
	function solutie(F, u0::Tuple, T; h = 0.01)
		u = float.(u0)
		t, U = [0.0], [u]
		for k in 1:round(Int, T / h)
			u = rk4(F, u, h)
			push!(t, k * h); push!(U, u)
		end
		t, [getindex.(U, i) for i in eachindex(u)]
	end

	# varianta pentru o singură ecuație x' = f(x):  t, X = solutie(f, 0.5, 20)
	function solutie(f, x0::Real, T; h = 0.01)
		x = float(x0)
		t, X = [0.0], [x]
		for k in 1:round(Int, T / h)
			x = rk4(f, x, h)
			push!(t, k * h); push!(X, x)
		end
		t, X
	end

	# minimul, maximul și perioada mărimii X după regimul tranzitoriu (t ≥ t0);
	# perioada este NaN dacă X nu oscilează
	function oscilatie(t, X, t0)
		τ, v = t[t .>= t0], X[t .>= t0]
		lo, hi = extrema(v)
		hi - lo < 1e-3 * max(abs(hi), 1) && return lo, hi, NaN
		m = (lo + hi) / 2
		urcari = [τ[k] for k in 2:length(v) if v[k-1] < m <= v[k]]
		lo, hi, length(urcari) >= 2 ? (urcari[end] - urcari[1]) / (length(urcari) - 1) : NaN
	end

	# cîmpul de vectori, desenat prin linii de curent (mai închise unde mișcarea este mai rapidă)
	cimp!(ax, F, xr, yr; densitate = 1.0) = streamplot!(ax, p -> Point2f(F(p)...),
		xr[1] .. xr[2], yr[1] .. yr[2], colormap = [:gray85, :gray35],
		linewidth = 1, arrow_size = 9, density = densitate)

	# punct de echilibru: negru dacă este stabil, alb dacă este instabil, gri în cazurile limită
	echilibru!(ax, x, y, stabil) = scatter!(ax, [x], [y], markersize = 16,
		strokewidth = 2.5, strokecolor = :black,
		color = isnothing(stabil) ? :gray60 : stabil ? :black : :white)

	# graficul lui f și linia de fază a ecuației x' = f(x) (Cursul 3): săgeți verzi unde
	# x crește, roșii unde x scade; stabilitatea echilibrelor ech din semnul lui f'
	function linie_de_faza!(ax, f, ech, xmin, xmax)
		xs = range(xmin, xmax, length = 400)
		hlines!(ax, 0, color = :gray50)
		lines!(ax, xs, f.(xs), color = :royalblue, linewidth = 3)
		repere = sort([xmin; filter(e -> xmin < e < xmax, ech); xmax])
		for (p, q) in zip(repere[1:end-1], repere[2:end])
			s = sign(f((p + q) / 2))
			s == 0 && continue
			L = 0.5 * (q - p)
			arrows2d!(ax, [(p + q) / 2 - s * L / 2], [0.0], [s * L], [0.0],
				color = s > 0 ? :seagreen : :crimson, shaftwidth = 3, tipwidth = 14, tiplength = 12)
		end
		for e in filter(e -> xmin <= e <= xmax, ech)
			d = (f(e + 1e-6) - f(e - 1e-6)) / 2e-6
			echilibru!(ax, e, 0.0, abs(d) < 1e-8 ? nothing : d < 0)
		end
	end
end

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 📈 Cursul 6 — Comportament periodic și bifurcații
"""

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000002
md"""
## 🎯 Obiectivele lecției
1. Să deosebim un **centru** (oscilații a căror amplitudine depinde de starea inițială) de un **ciclu limită** (oscilații robuste, cu amplitudine și perioadă proprii).
2. Să înțelegem două mecanisme care produc oscilații: **pomparea de energie** la amplitudini mici (clarinetul lui Rayleigh) și **feedback-ul negativ cu întîrziere** (hormoni, gene).
3. Să definim **bifurcația** și să citim o **diagramă de bifurcație** pentru un model care depinde de un parametru.
4. Să recunoaștem bifurcațiile **transcritică** și **șa-nod** și fenomenul de **histerezis** (pescuitul, comutatorul genetic lac).
5. Să vedem cum se naște un ciclu limită prin **bifurcația Hopf** (clarinetul, paradoxul îmbogățirii).
"""

# ╔═╡ bfbc972c-4703-4eb1-9405-cddb964d5bc2
md"""
## 🔁 0. Ce știm deja
* **Cursul 3:** pentru o singură variabilă, ``x' = f(x)``, soluțiile sînt **monotone**: nu pot traversa un echilibru, deci **nu pot oscila**. Stabilitatea echilibrelor se citește pe linia de fază sau din semnul lui ``f'(x^*)``.
* **Cursul 4:** în timp discret chiar și o singură variabilă poate oscila. În **diagrama de bifurcație** a modelului logistic discret am văzut cum, pe măsură ce ``r`` crește, echilibrul este înlocuit de cicluri de perioadă ``2, 4, 8, \dots`` și apoi de haos.
* **Cursul 5:** în plan traiectoriile se pot roti. Modelul Lotka–Volterra al rîșilor și iepurilor are **curbe închise** în jurul unui **centru**. La arcul cu frecare ``x' = v``, ``v' = -x - c\,v``, valoarea ``c = 0`` separă spiralele stabile (``c > 0``) de cele instabile (``c < 0``).

**Două întrebări pentru azi:**
1. 🔄 Bătăile inimii, ritmul circadian, secreția pulsatilă a hormonilor, ciclurile de aproximativ ``10`` ani ale rîșilor și iepurilor (Cursul 1) sînt oscilații care **se repetă cu aceeași amplitudine și perioadă**, deși sînt perturbate permanent. Ce fel de modele produc astfel de oscilații?
2. 🎚️ Ce se întîmplă cu echilibrele și cu oscilațiile unui model cînd un parametru se schimbă încet (pescuitul devine tot mai intens, mediul tot mai bogat)? De ce uneori schimbarea este **treptată și reversibilă**, iar alteori **bruscă și ireversibilă**?
"""

# ╔═╡ ad37b4e7-848d-41d4-b35b-995ad51c7e4d
md"""
## 🔄 1. Centre și cicluri limită

### 1.1. Două feluri de oscilații
Comparăm două sisteme din plan. Primul este modelul **rîși–iepuri** din Cursul 5, cu iepurii ``I`` pe axa orizontală și rîșii ``R`` pe cea verticală:
```math
I' = I - R\,I, \qquad R' = R\,I - R .
```
Al doilea este **clarinetul lui Rayleigh**, propus de Lord Rayleigh în secolul al XIX-lea pentru a explica sunetul instrumentelor de suflat. Aici ``x`` este deplasarea **ancii** (lamela care vibrează în muștiucul clarinetului), iar ``v`` este viteza ei:
```math
\begin{aligned}
x' &= v \\
v' &= -x + v - v^3 .
\end{aligned}
```
Ambele sisteme au cîte un singur echilibru în regiunea desenată: ``(1, 1)``, respectiv ``(0, 0)``. Pornim din puncte aflate la distanța ``\delta`` de echilibru și urmărim ce se întîmplă.

### 🎛️ Experimentați
Curbele gri sînt traiectorii care pornesc din alte puncte, pentru comparație.

Distanța de la echilibru ``\delta``: $(@bind δ_osc PlutoUI.Slider(0.05:0.05:2.5, show_value=true, default=0.5))
"""

# ╔═╡ 1e41598f-baa1-4e36-89bd-3e52818e5ec9
let
	lv(u) = (u[1] - u[1] * u[2], u[1] * u[2] - u[2])      # u = (I, R)
	ray(u) = (u[2], -u[1] + u[2] - u[2]^3)                 # u = (x, v)
	T = 40.0
	t, (I, R) = solutie(lv, (1 + δ_osc, 1.0), T)
	_, (X, V) = solutie(ray, (δ_osc, 0.0), T)
	loI, hiI, _ = oscilatie(t, I, 0)
	loX, hiX, perX = oscilatie(t, X, 15)

	fig = Figure(size = (1200, 850))
	ax1 = Axis(fig[1, 1], title = "Rîși și iepuri (Lotka–Volterra): un centru",
		subtitle = @sprintf("iepurii oscilează între %.2f și %.2f", loI, hiI),
		xlabel = "Iepuri (I)", ylabel = "Rîși (R)", aspect = 1, limits = (0, 4.5, 0, 4.5))
	ax2 = Axis(fig[1, 2], title = "Clarinetul lui Rayleigh: un ciclu limită",
		subtitle = @sprintf("după regimul tranzitoriu, x oscilează între %.2f și %.2f, cu perioada %.2f",
			loX, hiX, perX),
		xlabel = "Deplasarea ancii (x)", ylabel = "Viteza (v)", aspect = 1, limits = (-3, 3, -3, 3))
	cimp!(ax1, lv, (0, 4.5), (0, 4.5))
	cimp!(ax2, ray, (-3, 3), (-3, 3))
	for d in (0.25, 0.75, 1.5, 2.25)
		_, (Ig, Rg) = solutie(lv, (1 + d, 1.0), T)
		lines!(ax1, Ig, Rg, color = (:gray40, 0.6), linewidth = 1.5)
		_, (Xg, Vg) = solutie(ray, (d, 0.0), T)
		lines!(ax2, Xg, Vg, color = (:gray40, 0.6), linewidth = 1.5)
	end
	lines!(ax1, I, R, color = :crimson, linewidth = 3)
	lines!(ax2, X, V, color = :purple, linewidth = 3)
	scatter!(ax1, [1 + δ_osc], [1.0], color = :crimson, marker = :star5, markersize = 18)
	scatter!(ax2, [δ_osc], [0.0], color = :purple, marker = :star5, markersize = 18)
	echilibru!(ax1, 1, 1, nothing)
	echilibru!(ax2, 0, 0, false)

	ax3 = Axis(fig[2, 1], title = "Evoluția în timp", xlabel = "Timp (t)", ylabel = "Populație",
		limits = (0, T, 0, 1.4 * max(maximum(I), maximum(R))))
	lines!(ax3, t, I, color = :royalblue, linewidth = 2.5, label = "iepuri I(t)")
	lines!(ax3, t, R, color = :crimson, linewidth = 2.5, label = "rîși R(t)")
	axislegend(ax3, position = :rt, orientation = :horizontal, backgroundcolor = (:white, 0.85))
	ax4 = Axis(fig[2, 2], title = "Evoluția în timp", xlabel = "Timp (t)", ylabel = "x(t)",
		limits = (0, T, -3, 3))
	lines!(ax4, t, X, color = :purple, linewidth = 2.5)
	rowsize!(fig.layout, 2, Relative(0.3))
	fig
end

# ╔═╡ 7ff60562-f93a-4471-bda3-0bec3962c417
md"""
### 🤔 Ce am observat?
* **Rîși și iepuri:** fiecare stare inițială are **propria** curbă închisă, iar amplitudinea oscilațiilor este dictată de starea inițială. O perturbare (o iarnă grea, o epidemie) mută sistemul pe altă curbă, iar sistemul nu mai revine: oscilațiile „țin minte” toate perturbările. Mai mult, o modificare oricît de mică a ecuațiilor distruge curbele închise (secțiunea 4.2). Un astfel de comportament este **fragil** (despre sistemele de acest tip, numite conservative, vom vorbi în Cursul 10).
* **Clarinetul:** pentru orice ``\delta > 0`` traiectoria se apropie de **aceeași** curbă închisă: din interior spiralează spre exterior, din exterior spre interior. Amplitudinea și perioada oscilațiilor sînt proprietăți ale **sistemului**, nu ale stării inițiale. După o perturbare, sistemul revine la aceleași oscilații.

> 🔁 Un **ciclu limită** este o traiectorie închisă **izolată**: în vecinătatea ei nu mai există alte traiectorii închise. Dacă traiectoriile vecine se apropie de el, ciclul limită este **stabil**: este un **atractor**, la fel ca un echilibru stabil. Oscilațiile robuste din biologie (inima, ritmul circadian, ciclurile hormonale) se modelează prin cicluri limită stabile.

### 1.2. De ce oscilează clarinetul?
Scriem ecuația vitezei sub forma arcului cu frecare din Cursul 5:
```math
v' = -x - c(v)\,v, \qquad \text{unde } c(v) = v^2 - 1 .
```
Coeficientul de frecare depinde acum de viteză:
* pentru ``\lvert v\rvert < 1`` avem ``c(v) < 0``, adică **frecare negativă**: ancia primește energie de la aerul suflat de muzician, ca un leagăn împins la fiecare trecere. Lîngă echilibru avem o **spirală instabilă**: oscilațiile mici cresc.
* pentru ``\lvert v\rvert > 1`` avem ``c(v) > 0``, o frecare obișnuită, care disipă energia și împiedică oscilațiile să crească la nesfîrșit.

Ciclul limită este curba pe care, de-a lungul unei perioade, energia primită este egală cu energia pierdută.

### ✅ Verificare rapidă
Energia ancii este ``E = \frac{1}{2}\left(x^2 + v^2\right)``. Calculați ``E'`` cu regula lanțului, folosind ecuațiile clarinetului. Unde crește energia și unde scade? Poate ciclul limită să rămînă în întregime în banda ``\lvert v\rvert < 1``?
"""

# ╔═╡ 672a0d81-5c6c-4906-a828-8c847dd4b80b
hint(md"``E' = x\,x' + v\,v' = xv + v\,(-x + v - v^3) = v^2\,(1 - v^2)``. Energia crește cînd ``\lvert v\rvert < 1`` și scade cînd ``\lvert v\rvert > 1``. După o perioadă energia revine la valoarea inițială, deci de-a lungul ciclului ``E'`` trebuie să fie uneori pozitivă și alteori negativă: ciclul limită traversează dreptele ``v = \pm 1``. Pe figură, ``v`` ajunge pînă la aproximativ ``\pm 1.16``.")

# ╔═╡ 39bf943b-72f1-42b4-ab31-a1a8840af5d9
md"""
> 🧭 **Teorema Poincaré–Bendixson.** În plan, dacă o regiune mărginită care nu conține echilibre „prinde” traiectoriile (ele intră în regiune, dar nu mai ies), atunci în ea există un ciclu limită. Pentru clarinet, regiunea este un inel în jurul originii: lîngă origine traiectoriile pleacă spre exterior, iar departe de ea se întorc spre interior. O consecință surprinzătoare: un sistem continuu cu **două** variabile **nu poate fi haotic**. Pentru haos în timp continuu sînt necesare cel puțin trei variabile (E. Lorenz, 1963), spre deosebire de timpul discret (Cursul 4), unde ajunge una singură.

Putem acum completa lista comportamentelor pe termen lung (**atractori**) întîlnite pînă acum:

| Atractor | Comportament pe termen lung | Exemple |
| :--- | :--- | :--- |
| 📍 **echilibru stabil** | stare constantă | populația logistică ajunge la capacitatea de suport (Cursul 3) |
| 🔄 **ciclu limită stabil** | oscilații periodice robuste | clarinetul, bătăile inimii, ciclurile hormonale |
| 🌪️ **atractor haotic** | oscilații neregulate, mărginite, extrem de sensibile la condițiile inițiale | modelul logistic discret cu ``r = 3.9`` (Cursul 4) |
"""

# ╔═╡ f1287f65-d265-45d6-8eb6-e8f1cae80852
md"""
## 🧬 2. Feedback negativ cu întîrziere

În Cursul 1 am văzut că **feedback-ul negativ** stabilizează: termostatul ține temperatura constantă, insulina readuce glicemia la normal. Acest lucru este adevărat cînd reacția este **rapidă**. Dacă reacția vine **cu întîrziere**, sistemul corectează prea tîrziu și prea mult, iar corecția depășește ținta. 🚿 Gîndiți-vă la un duș la care apa caldă ajunge cu întîrziere: apa este rece, rotim robinetul spre cald, apa este tot rece, îl mai rotim, apoi apa devine fierbinte, rotim înapoi... iar temperatura oscilează.

### 2.1. Un model hormonal
Hormonul ``H`` secretat de hipotalamus (GnRH) stimulează hipofiza să secrete hormonul ``P`` (LH), care stimulează gonadele să secrete hormonul ``G`` (testosteronul). Bucla se închide: un nivel ridicat de ``G`` **inhibă** secreția lui ``H``.
```math
\begin{aligned}
H' &= \frac{1}{1 + G^n} - k\,H \\
P' &= H - k\,P \\
G' &= P - k\,G
\end{aligned}
\qquad\qquad
H \xrightarrow{\ +\ } P \xrightarrow{\ +\ } G \xrightarrow{\ -\ } H
```
* fiecare hormon este eliminat din sînge cu aceeași rată ``k = 0.2`` (termenii ``-kH``, ``-kP``, ``-kG``);
* termenul ``1/(1 + G^n)`` este o **funcție Hill** descrescătoare (Cursul 1): cînd ``G`` este mare, secreția lui ``H`` scade. Exponentul ``n``, numit **cooperativitate**, măsoară cît de bruscă este inhibiția: pentru ``n`` mare, ``1/(1 + G^n)`` se comportă ca un întrerupător;
* **întîrzierea** vine din lanțul cu trei etape: ``G`` „află” că ``H`` a crescut abia după ce ``H`` a stimulat ``P``, iar ``P`` a stimulat ``G``.

Modelele de acest tip se numesc **oscilatoare Goodwin** (B. Goodwin, 1965).

### 🎛️ Experimentați
Cooperativitatea ``n``: $(@bind n_horm PlutoUI.Slider(1:0.5:14, show_value=true, default=10))
"""

# ╔═╡ ecd0a0e6-78ee-435a-9c71-fb9c679a4d7a
let
	k = 0.2
	F(u) = (1 / (1 + u[3]^n_horm) - k * u[1], u[1] - k * u[2], u[2] - k * u[3])
	T = 300.0
	t, (H, P, G) = solutie(F, (0.2, 0.2, 0.2), T; h = 0.02)
	# echilibrul: G*(1 + G*^n) = 1/k³, P* = kG*, H* = k²G*; el este instabil (și apare un
	# ciclu limită) dacă panta buclei de feedback depășește 8k³ (liniarizare, Cursul 9)
	Gs = find_zero(G -> G * (1 + G^n_horm) - 1 / k^3, (0.0, 1 / k^3))
	instabil = n_horm * Gs^(n_horm - 1) / (1 + Gs^n_horm)^2 > 8k^3
	lo, hi, per = oscilatie(t, G, 200)
	sub = instabil ?
		@sprintf("echilibrul este instabil ⟹ ciclu limită: G oscilează între %.2f și %.2f, cu perioada %.1f",
			lo, hi, per) :
		@sprintf("echilibrul G* ≈ %.2f este stabil: oscilațiile se sting", Gs)

	# primul vîrf, datorat stării inițiale depărtate de echilibru, iese din grafic
	dupa = t .>= 50
	fig = Figure(size = (1200, 500))
	ax1 = Axis(fig[1, 1], title = "Concentrațiile hormonilor (n = $(n_horm))", subtitle = sub,
		xlabel = "Timp (t)", ylabel = "Concentrație", limits = (0, T, 0, 1.6 * maximum(G[dupa])))
	lines!(ax1, t, H, color = :royalblue, linewidth = 2, label = "H (hipotalamus)")
	lines!(ax1, t, P, color = :darkorange, linewidth = 2, label = "P (hipofiză)")
	lines!(ax1, t, G, color = :crimson, linewidth = 2.5, label = "G (gonade)")
	axislegend(ax1, position = :rt, orientation = :horizontal, backgroundcolor = (:white, 0.85))

	ax2 = Axis3(fig[1, 2], title = "Spațiul stărilor, după t = 50 (culoarea arată timpul)",
		xlabel = "H", ylabel = "P", zlabel = "G")
	lines!(ax2, H[dupa], P[dupa], G[dupa], color = t[dupa], colormap = :viridis, linewidth = 2)
	scatter!(ax2, [k^2 * Gs], [k * Gs], [Gs], color = instabil ? :white : :black,
		strokewidth = 2, strokecolor = :black, markersize = 14)
	fig
end

# ╔═╡ 6cf2aef4-e72e-422d-b9d9-188d863172be
md"""
### 🤔 Ce am observat?
* Pentru ``n`` mic, oscilațiile se sting: traiectoria spiralează spre echilibru, iar nivelurile hormonilor devin constante.
* Pentru ``n`` mai mare decît aproximativ ``8.1``, echilibrul devine instabil și apare un **ciclu limită**: secreția hormonilor devine **pulsatilă**. Valoarea critică se poate calcula prin liniarizare (Cursul 9).
* Maximele se succed în ordinea ``H``, ``P``, ``G``: fiecare hormon rămîne în urma celui care îl stimulează, așa cum, în Cursul 5, maximul rîșilor urma maximului iepurilor.
* Trecerea de la echilibru la oscilații cînd ``n`` depășește o valoare critică este o **bifurcație** (secțiunea 3), mai precis o **bifurcație Hopf** (secțiunea 4).

### 2.2. Cît de multă întîrziere este necesară?
Generalizăm modelul la un lanț cu ``m`` etape: ``x_1`` stimulează ``x_2``, ..., ``x_{m-1}`` stimulează ``x_m``, iar ``x_m`` inhibă ``x_1``:
```math
x_1' = \frac{1}{1 + x_m^n} - k\,x_1, \qquad x_i' = x_{i-1} - k\,x_i \quad (i = 2, \dots, m).
```
Pentru fiecare ``m`` și ``n`` simulăm sistemul mult timp și măsurăm cît de mari sînt oscilațiile care rămîn după regimul tranzitoriu.
"""

# ╔═╡ 866002fb-4d84-45ac-843d-22a888db5be7
# amplitudinea relativă (max − min)/(max + min) a oscilațiilor variabilei x_m, după regimul
# tranzitoriu, pentru m = 2, 3, 4, 5 etape și n între 1 și 14
date_lant = let
	k = 0.2
	function amplitudine(::Val{m}, n) where {m}
		F(u) = ntuple(i -> i == 1 ? 1 / (1 + u[m]^n) - k * u[1] : u[i-1] - k * u[i], Val(m))
		u = ntuple(_ -> 0.5, Val(m))
		for _ in 1:30000
			u = rk4(F, u, 0.05)
		end
		lo, hi = Inf, -Inf
		for _ in 1:8000
			u = rk4(F, u, 0.05)
			lo, hi = min(lo, u[m]), max(hi, u[m])
		end
		(hi - lo) / (hi + lo)
	end
	ns, ms = 1:0.25:14, 2:5
	(; ns, ms, A = [amplitudine(Val(m), n) for m in ms, n in ns])
end

# ╔═╡ 11a7f979-4026-486e-92c8-834077357a1f
let
	(; ns, ms, A) = date_lant
	culori = [:gray50, :crimson, :royalblue, :seagreen]
	fig = Figure(size = (1000, 450))
	ax = Axis(fig[1, 1], title = "Oscilațiile unei bucle de feedback negativ cu m etape",
		subtitle = "amplitudinea relativă (max − min)/(max + min) a variabilei xₘ, după regimul tranzitoriu",
		xlabel = "Cooperativitatea n", ylabel = "Amplitudinea relativă", limits = (1, 14, -0.02, 1))
	for (j, m) in enumerate(ms)
		lines!(ax, ns, A[j, :], color = culori[j], linewidth = 3, label = "m = $m etape")
	end
	vlines!(ax, n_horm, color = :black, linestyle = :dash, label = "n ales în secțiunea 2.1")
	axislegend(ax, position = :lt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ af3eceb0-6556-49d5-b0c7-acbfec2480ff
md"""
**Ce observăm:**
* Cu **două** etape nu apar niciodată oscilații susținute, oricît de mare ar fi ``n``: doar oscilații amortizate (o spirală stabilă).
* Cu trei etape oscilațiile apar pentru ``n > 8.1``, cu patru pentru ``n > 4``, cu cinci pentru ``n > 2.9``.
* **Cu cît bucla de feedback negativ are mai multe etape (deci mai multă întîrziere), cu atît oscilațiile apar mai ușor.** Pentru rate de eliminare egale, pragul este aproximativ ``\left(\cos\frac{\pi}{m}\right)^{-m}``: ``8`` pentru ``m = 3``, ``4`` pentru ``m = 4``, ``2.9`` pentru ``m = 5``.

### 2.3. Ceasul segmentării: gena Hes1
În embrionul vertebratelor, **somitele** (din care se vor forma vertebrele și mușchii spatelui) apar una cîte una, în ritmul unui „ceas” molecular. Una dintre piesele acestui ceas este gena **Hes1**, a cărei expresie oscilează cu o perioadă de aproximativ ``2`` ore (Hirata și colaboratorii, 2002). Mecanismul este din nou un feedback negativ cu întîrziere: proteina Hes1 **inhibă transcripția propriei gene**, iar între transcripție și apariția proteinei trec zeci de minute (transcripția, transportul ARN-ului mesager din nucleu, translația).

Un model cu trei variabile: proteina ``X``, ARN-ul mesager ``Y`` și un al treilea factor ``Z``, care se leagă de proteină (termenul ``-AXZ``), cu timpul măsurat în minute:
```math
\begin{aligned}
X' &= -A\,X Z + B\,Y - C\,X \\
Y' &= \frac{E}{1 + X^2} - D\,Y \\
Z' &= -A\,X Z + \frac{F}{1 + X^2} - G\,Z .
\end{aligned}
```
Termenii ``E/(1 + X^2)`` și ``F/(1 + X^2)`` descriu inhibiția genelor de către proteină.
"""

# ╔═╡ 4d4ae524-a518-46d9-952c-d27415be1d79
let
	A, B, C, D, E, F, G = 0.022, 0.3, 0.031, 0.028, 0.5, 20.0, 0.3
	hes1(u) = (-A * u[1] * u[3] + B * u[2] - C * u[1],
		E / (1 + u[1]^2) - D * u[2],
		-A * u[1] * u[3] + F / (1 + u[1]^2) - G * u[3])
	t, (X, Y, Z) = solutie(hes1, (0.2, 0.2, 0.2), 1200.0; h = 0.05)
	_, _, per = oscilatie(t, X, 600)

	fig = Figure(size = (1200, 480))
	ax1 = Axis(fig[1, 1], title = "Modelul Hes1",
		subtitle = @sprintf("perioada oscilațiilor ≈ %.0f minute ≈ %.1f ore", per, per / 60),
		xlabel = "Timp (ore)", ylabel = "Concentrație")
	lines!(ax1, t ./ 60, X, color = :crimson, linewidth = 2.5, label = "proteina X")
	lines!(ax1, t ./ 60, Y, color = :royalblue, linewidth = 2.5, label = "ARNm Y")
	lines!(ax1, t ./ 60, Z, color = :gray50, linewidth = 2, label = "factorul Z")
	axislegend(ax1, position = :rt, backgroundcolor = (:white, 0.85))
	ax2 = Axis3(fig[1, 2], title = "Ciclul limită (după primele 5 ore)",
		xlabel = "X", ylabel = "Y", zlabel = "Z")
	dupa = t .> 300
	lines!(ax2, X[dupa], Y[dupa], Z[dupa], color = :purple, linewidth = 2)
	fig
end

# ╔═╡ aea3cccd-b20e-4b7b-b5d8-a76b71c1386b
md"""
**Ce observăm:** după un scurt regim tranzitoriu, toate cele trei mărimi oscilează cu o perioadă de aproximativ ``2`` ore, ca în experimente. Bucla este: multă proteină ``\Rightarrow`` gena este inhibată ``\Rightarrow`` ARNm scade ``\Rightarrow`` proteina scade ``\Rightarrow`` gena este eliberată ``\Rightarrow`` ARNm crește ``\Rightarrow`` ... Maximul ARNm precede maximul proteinei.

### ✅ Verificare rapidă
1. Un medicament reduce cooperativitatea inhibiției în modelul hormonal de la ``n = 10`` la ``n = 6``. Ce se întîmplă cu secreția hormonilor?
2. De ce ecuația ``x' = \dfrac{1}{1 + x^n} - k\,x`` (o buclă de feedback negativ cu o singură etapă) nu poate oscila, oricît de mare ar fi ``n``?
"""

# ╔═╡ 2a217a54-02cf-4500-a29e-9375a75eb4d9
hint(md"1. Deoarece ``6 < 8.1``, echilibrul redevine stabil: pulsurile se sting, iar nivelurile hormonilor devin constante (verificați cu slider-ul din secțiunea 2.1). 2. Este o ecuație cu o singură variabilă, deci soluțiile ei sînt monotone (Cursul 3). Biologic: fără întîrziere, corecția vine la timp și nu depășește ținta.")

# ╔═╡ a7fb6f10-1998-4827-a6d4-9af07fbe67d1
md"""
## 🦋 3. Bifurcații

Toate modelele depind de **parametri**: efortul de pescuit, cooperativitatea ``n``, rata de degradare ``k``... Cînd un parametru ``\mu`` variază puțin, echilibrele se deplasează puțin și, de obicei, nimic esențial nu se schimbă. La anumite valori ale parametrului însă **numărul** echilibrelor sau **stabilitatea** lor se schimbă, ori apare un ciclu limită. Spunem atunci că are loc o **bifurcație**, iar valoarea respectivă a lui ``\mu`` se numește **punct de bifurcație**.

**Diagrama de bifurcație** reprezintă echilibrele ``x^*`` în funcție de parametrul ``\mu``: cu linie continuă echilibrele stabile, cu linie punctată cele instabile. Este analogul diagramei de bifurcație din Cursul 4: pe fiecare verticală ``\mu = `` const. citim linia de fază pentru acea valoare a parametrului.

### 3.1. Bifurcația transcritică: pescuit cu efort constant
Reluăm exemplul din Cursul 3: o populație de pești (``r = 0.2`` pe an, ``K = 1000`` de tone) din care se pescuiește în fiecare an fracțiunea ``E`` din stoc (**efortul** de pescuit):
```math
x' = 0.2\,x\left(1 - \frac{x}{1000}\right) - E\,x .
```
Echilibrele sînt ``x^* = 0`` și ``x^* = 1000\left(1 - \dfrac{E}{0.2}\right)``. Cînd ``E`` crește, al doilea echilibru coboară, iar pentru ``E = 0.2`` el **se ciocnește** cu ``0``. Pentru ``E > 0.2`` cele două echilibre trec unul prin celălalt și **își schimbă stabilitatea**: ``0`` devine stabil, iar populația dispare. Echilibrele negative nu au sens biologic (zona gri), dar ne ajută să vedem ce se întîmplă.

Efortul de pescuit ``E``: $(@bind E_pesc PlutoUI.Slider(0.0:0.01:0.3, show_value=true, default=0.1))
"""

# ╔═╡ 4cb83737-82f7-46f2-99e2-1f49222ada4a
let
	r = 0.2
	f(x) = r * x * (1 - x / 1000) - E_pesc * x
	x_ech = 1000 * (1 - E_pesc / r)
	critic = isapprox(E_pesc, r)

	fig = Figure(size = (1300, 440))
	ax1 = Axis(fig[1, 1], title = "Linia de fază (E = $(E_pesc))", xlabel = "Stocul x (tone)",
		ylabel = "x'", limits = (-500, 1100, -60, 60))
	vspan!(ax1, -500, 0, color = (:gray, 0.15))
	linie_de_faza!(ax1, f, [0.0, x_ech], -500, 1100)

	ax2 = Axis(fig[1, 2], title = "Diagrama de bifurcație",
		subtitle = "continuu: echilibru stabil;  punctat: instabil",
		xlabel = "Efortul de pescuit E", ylabel = "Echilibrele x*", limits = (0, 0.3, -500, 1100))
	hspan!(ax2, -500, 0, color = (:gray, 0.15))
	lines!(ax2, [0, r], [1000, 0], color = :black, linewidth = 3)
	lines!(ax2, [r, 0.3], [0, -500], color = :black, linewidth = 2, linestyle = :dash)
	lines!(ax2, [0, r], [0, 0], color = :black, linewidth = 2, linestyle = :dash)
	lines!(ax2, [r, 0.3], [0, 0], color = :black, linewidth = 3)
	text!(ax2, r, 0, text = "E = 0.2", align = (:left, :bottom), offset = (6, 6))
	vlines!(ax2, E_pesc, color = :darkorange, linewidth = 2)
	echilibru!(ax2, E_pesc, 0.0, critic ? nothing : E_pesc > r)
	echilibru!(ax2, E_pesc, x_ech, critic ? nothing : E_pesc < r)

	ax3 = Axis(fig[1, 3], title = "Stocul în timp",
		subtitle = @sprintf("captura la echilibru: E·x* = %.1f tone pe an", E_pesc * max(x_ech, 0)),
		xlabel = "Timp (ani)", ylabel = "x(t) (tone)", limits = (0, 100, 0, 1100))
	hlines!(ax3, max(x_ech, 0), color = :gray50, linestyle = :dash)
	for (x0, c) in ((1000.0, :crimson), (200.0, :royalblue))
		t, X = solutie(f, x0, 100; h = 0.05)
		lines!(ax3, t, X, color = c, linewidth = 2.5, label = "x₀ = $(Int(x0))")
	end
	axislegend(ax3, position = :rt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ 53cd73af-6c92-4bb4-97ac-ef4a3798e1f3
md"""
**Ce observăm:**
* Pentru ``E < 0.2`` stocul se stabilizează la ``x^* = 1000(1 - 5E)``, iar captura anuală este ``E\,x^*``. Ea este maximă pentru ``E = 0.1`` (``50`` de tone pe an), exact exemplul din Cursul 3 (Exercițiul 1).
* Cînd ``E`` se apropie de ``0.2``, echilibrul coboară **treptat** spre ``0``: declinul stocului se vede din timp și poate fi corectat reducînd efortul. În plus, revenirea la echilibru devine tot mai lentă, deoarece ``f'(x^*) = E - 0.2`` se apropie de ``0`` (fenomenul se numește **încetinire critică** și este un semnal de avertizare).
* Pentru ``E \ge 0.2`` populația dispare.

### 3.2. Bifurcația șa-nod: pescuit cu cotă fixă
În Exercițiul 4 din Cursul 3 se pescuia o cantitate **fixă** ``H`` (cota) în fiecare an:
```math
x' = 0.5\,x\left(1 - \frac{x}{100}\right) - H .
```
Echilibrele sînt rădăcinile ecuației ``0.005\,x^2 - 0.5\,x + H = 0``:
```math
x^*_{\pm} = 50 \pm \sqrt{2500 - 200\,H}, \qquad \text{definite doar pentru } H \le 12.5 .
```
Echilibrul ``x^*_+`` este stabil (stocul „sănătos”), iar ``x^*_-`` este instabil: un **prag** sub care populația se prăbușește, ca la efectul Allee. Cînd ``H`` crește, cele două echilibre se apropie, pentru ``H = 12.5`` se contopesc în ``x = 50``, iar pentru ``H > 12.5`` **dispar** amîndouă. O astfel de ciocnire, în care un echilibru stabil (un **nod**) și unul instabil (în plan, un **punct șa**) se anihilează, se numește **bifurcație șa-nod**.

Cota anuală ``H``: $(@bind H_pesc PlutoUI.Slider(0.0:0.25:16.0, show_value=true, default=10.0))
"""

# ╔═╡ 30a7f680-64d8-4667-8288-80b0cbfe4e8d
let
	f(x) = 0.5 * x * (1 - x / 100) - H_pesc
	Δ = 2500 - 200 * H_pesc
	ech = Δ >= 0 ? [50 - sqrt(Δ), 50 + sqrt(Δ)] : Float64[]

	fig = Figure(size = (1300, 440))
	ax1 = Axis(fig[1, 1], title = "Linia de fază (H = $(H_pesc))", xlabel = "Stocul x (tone)",
		ylabel = "x'", limits = (0, 110, -17, 14))
	linie_de_faza!(ax1, f, ech, 0, 110)

	ax2 = Axis(fig[1, 2], title = "Diagrama de bifurcație",
		subtitle = "continuu: echilibru stabil;  punctat: instabil",
		xlabel = "Cota H (tone pe an)", ylabel = "Echilibrele x*", limits = (0, 16, 0, 110))
	Hs = range(0, 12.5, length = 200)
	lines!(ax2, Hs, 50 .+ sqrt.(max.(2500 .- 200 .* Hs, 0)), color = :black, linewidth = 3)
	lines!(ax2, Hs, 50 .- sqrt.(max.(2500 .- 200 .* Hs, 0)), color = :black, linewidth = 2,
		linestyle = :dash)
	text!(ax2, 12.5, 50, text = "H = 12.5", align = (:left, :center), offset = (8, 0))
	vlines!(ax2, H_pesc, color = :darkorange, linewidth = 2)
	for (i, e) in enumerate(ech)
		echilibru!(ax2, H_pesc, e, Δ == 0 ? nothing : i == 2)
	end

	ax3 = Axis(fig[1, 3], title = "Stocul în timp", xlabel = "Timp (ani)", ylabel = "x(t) (tone)",
		limits = (0, 80, 0, 110))
	for (x0, c) in ((100.0, :crimson), (60.0, :darkorange), (20.0, :royalblue))
		t, X = solutie(f, x0, 80; h = 0.02)
		n = something(findfirst(<(0), X), length(X))
		lines!(ax3, t[1:n], X[1:n], color = c, linewidth = 2.5, label = "x₀ = $(Int(x0))")
		if x0 == 100 && n < length(X)
			ax3.subtitle = @sprintf("pornind de la x₀ = 100, stocul dispare după %.0f ani", t[n])
		end
	end
	axislegend(ax3, position = :rt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ 551bf8e5-e871-41d7-b24a-5baf75052a78
md"""
**Ce observăm:**
* Pînă la ``H = 12.5`` stocul rămîne mare (cel puțin ``50`` de tone): **nivelul stocului nu ne avertizează** că ne apropiem de pragul critic.
* Pentru ``H > 12.5`` avem ``x' < 0`` peste tot: prăbușirea este **bruscă** și pornește de la orice nivel al stocului, oricît de sănătos ar părea.
* Pentru ``H`` puțin peste ``12.5`` (încercați ``H = 12.75``), populația zăbovește mult timp în jurul valorii ``x = 50``, unde ``x'`` este aproape zero („fantoma” echilibrelor dispărute), apoi se prăbușește. Această aparentă stabilitate dă o falsă senzație de siguranță.
* Valoarea critică ``12.5`` este exact maximul creșterii naturale ``0.5\,x(1 - x/100)``: nu putem pescui mai mult decît se poate reface populația.

> 🐟 **Tipul bifurcației contează.** Cu efort constant (bifurcație transcritică) supra-exploatarea duce la un declin **treptat**, vizibil și reversibil. Cu cotă fixă (bifurcație șa-nod) stocul pare sănătos pînă în momentul în care se prăbușește **brusc**. Un exemplu celebru este prăbușirea stocului de cod de lîngă Newfoundland, la începutul anilor '90, după care pescuitul a fost interzis în 1992.

### 3.3. Histerezis: comutatorul genetic lac
În Cursul 5 (secțiunea 4.1) am studiat operonul lac al bacteriei *E. coli*:
```math
x' = \frac{a + x^2}{1 + x^2} - k\,x, \qquad a = 0.01 ,
```
unde ``x`` este concentrația de permează. Pentru ``k = 0.4`` am găsit trei echilibre: „**oprit**” (``x^* \approx 0.027``, stabil), un **prag** (``x^* \approx 0.465``, instabil) și „**pornit**” (``x^* \approx 2.008``, stabil). Acum variem rata de degradare ``k``. Un ``k`` mai mic înseamnă că producția cîștigă în fața degradării, ceea ce se întîmplă cînd în mediu crește concentrația de lactoză.

Echilibrele sînt intersecțiile dintre curba de producție și dreapta de degradare. Pentru diagrama de bifurcație este mai simplu să exprimăm ``k`` în funcție de echilibru: ``k = \dfrac{a + x^2}{x\,(1 + x^2)}``.

Rata de degradare ``k``: $(@bind k_lac PlutoUI.Slider(0.10:0.01:0.70, show_value=true, default=0.40))
"""

# ╔═╡ b7f5c448-b5c7-4d19-9d0d-8897e5708ab0
# diagrama de bifurcație a comutatorului lac: echilibrele satisfac k = (a + x²)/(x(1 + x²)),
# iar un echilibru este stabil acolo unde această funcție de x este descrescătoare;
# întoarce valorile lui k în care au loc cele două bifurcații șa-nod
function diagrama_lac!(ax; a = 0.01)
	g(x) = (a + x^2) / (x * (1 + x^2))
	x1, x2 = find_zeros(x -> (g(x + 1e-6) - g(x - 1e-6)) / 2e-6, 0.02, 3)
	for (xa, xb, stabil) in ((0.008, x1, true), (x1, x2, false), (x2, 12.0, true))
		xs = range(xa, xb, length = 300)
		lines!(ax, g.(xs), xs, color = :black, linewidth = stabil ? 3 : 2,
			linestyle = stabil ? :solid : :dash)
	end
	g(x1), g(x2)
end

# ╔═╡ 48ab4099-ca5b-47c6-a65b-6031fcb25efa
let
	a = 0.01
	productie(x) = (a + x^2) / (1 + x^2)
	f(x) = productie(x) - k_lac * x
	ech = find_zeros(f, 0, 15)
	stabil(e) = (f(e + 1e-6) - f(e - 1e-6)) < 0
	stare = length(ech) == 3 ? @sprintf("bistabil: „oprit” (x ≈ %.3f) sau „pornit” (x ≈ %.2f)", ech[1], ech[3]) :
		ech[end] < 0.5 ? @sprintf("un singur echilibru: „oprit” (x ≈ %.3f)", ech[1]) :
		@sprintf("un singur echilibru: „pornit” (x ≈ %.2f)", ech[end])

	fig = Figure(size = (1200, 480))
	xmax = max(3.0, 1.15 * maximum(ech))
	ax1 = Axis(fig[1, 1], title = "Producție și degradare (k = $(k_lac))", subtitle = stare,
		xlabel = "Concentrația de permează (x)", ylabel = "Rata", limits = (0, xmax, -0.02, 1.1))
	xs = range(0, xmax, length = 400)
	lines!(ax1, xs, productie.(xs), color = :royalblue, linewidth = 3, label = "producție (a + x²)/(1 + x²)")
	lines!(ax1, xs, k_lac .* xs, color = :crimson, linewidth = 3, label = "degradare k x")
	for e in ech
		echilibru!(ax1, e, productie(e), stabil(e))
	end
	axislegend(ax1, position = :rb, backgroundcolor = (:white, 0.85))

	ax2 = Axis(fig[1, 2], title = "Diagrama de bifurcație",
		subtitle = "continuu: echilibru stabil;  punctat: instabil",
		xlabel = "Rata de degradare k", ylabel = "Echilibrele x* (scară logaritmică)",
		yscale = log10, limits = (0.1, 0.7, 0.005, 12))
	k1, k2 = diagrama_lac!(ax2)
	vspan!(ax2, k1, k2, color = (:gold, 0.15))
	text!(ax2, (k1 + k2) / 2, 10, text = "zona bistabilă", align = (:center, :top))
	vlines!(ax2, k_lac, color = :darkorange, linewidth = 2)
	for e in ech
		echilibru!(ax2, k_lac, e, stabil(e))
	end
	fig
end

# ╔═╡ 25fdc3a9-caa3-49dc-89fc-b63ed025f11d
md"""
Diagrama are forma literei **S**, cu două puncte de întoarcere, ``k \approx 0.198`` și ``k \approx 0.505``: în fiecare dintre ele are loc o **bifurcație șa-nod**, în care pragul se ciocnește cu unul dintre echilibrele stabile. Între ele (zona galbenă) sistemul este **bistabil**.

**Un experiment.** Pornim cu operonul oprit și ``k = 0.7``. Scădem foarte încet ``k`` pînă la ``0.1``, apoi îl creștem la fel de încet înapoi pînă la ``0.7``. Urmărim starea celulei pe diagrama de bifurcație.
"""

# ╔═╡ 907565a1-48e0-4481-981e-63bcfbd3b2ac
let
	a = 0.01
	T = 30000.0
	# k scade liniar de la 0.7 la 0.1 în prima jumătate a timpului, apoi crește înapoi la 0.7
	k_t(t) = t < T / 2 ? 0.7 - 1.2 * t / T : 0.1 + 1.2 * (t - T / 2) / T
	# adăugăm timpul ca a doua variabilă de stare, ca să putem folosi solutie
	F(u) = ((a + u[1]^2) / (1 + u[1]^2) - k_t(u[2]) * u[1], 1.0)
	t, (X, _) = solutie(F, (0.015, 0.0), T; h = 0.2)
	K = k_t.(t)
	jos = t .<= T / 2

	fig = Figure(size = (1200, 500))
	ax1 = Axis(fig[1, 1], title = "Histerezis: drumul celulei pe diagrama de bifurcație",
		xlabel = "Rata de degradare k", ylabel = "Permeaza x (scară logaritmică)",
		yscale = log10, limits = (0.1, 0.7, 0.005, 12))
	k1, k2 = diagrama_lac!(ax1)
	lines!(ax1, K[jos], X[jos], color = :darkorange, linewidth = 4, label = "k scade")
	lines!(ax1, K[.!jos], X[.!jos], color = :teal, linewidth = 4, label = "k crește")
	text!(ax1, k1, 0.5, text = @sprintf("pornire bruscă\nla k ≈ %.2f", k1), align = (:left, :center),
		offset = (8, 0), color = :darkorange)
	text!(ax1, k2, 0.15, text = @sprintf("oprire bruscă\nla k ≈ %.2f", k2), align = (:right, :center),
		offset = (-8, 0), color = :teal)
	axislegend(ax1, position = :rt, backgroundcolor = (:white, 0.85))

	dreapta = fig[1, 2] = GridLayout()
	axk = Axis(dreapta[1, 1], title = "Același experiment în timp", ylabel = "k(t)",
		limits = (0, T, 0, 0.8))
	ax2 = Axis(dreapta[2, 1], xlabel = "Timp (t)", ylabel = "Permeaza x(t)", limits = (0, T, -0.3, 11),
		xticks = (0:5000:25000, string.(0:5000:25000)))
	linkxaxes!(axk, ax2)
	hidexdecorations!(axk, grid = false)
	rowsize!(dreapta, 1, Relative(0.3))
	lines!(axk, t[jos], K[jos], color = :darkorange, linewidth = 3)
	lines!(axk, t[.!jos], K[.!jos], color = :teal, linewidth = 3)
	lines!(ax2, t[jos], X[jos], color = :darkorange, linewidth = 3)
	lines!(ax2, t[.!jos], X[.!jos], color = :teal, linewidth = 3)
	fig
end

# ╔═╡ 6801d4b4-9f0f-452d-8fd1-84047b8a1b1b
md"""
**Ce observăm:**
* Cînd ``k`` scade, celula rămîne „oprită” mult după ce starea „pornit” a devenit posibilă (``k = 0.505``) și pornește brusc abia la ``k \approx 0.198``, cînd echilibrul „oprit” dispare printr-o bifurcație șa-nod.
* Cînd ``k`` crește înapoi, celula rămîne „pornită” pînă la ``k \approx 0.505``, unde dispare echilibrul „pornit”.
* Pentru aceeași valoare a parametrului, de exemplu ``k = 0.35``, starea celulei depinde de **drumul** pe care s-a ajuns acolo. Acest fenomen se numește **histerezis** (din greacă: „rămînere în urmă”): sistemul are **memorie**.

Experimental, histerezisul operonului lac a fost observat de A. Novick și M. Weiner (1957) și, mai precis, de E. Ozbudak și colaboratorii (2004): bacteriile care au fost expuse anterior la un inductor rămîn „pornite” la concentrații la care bacteriile neexpuse rămîn „oprite”. Același mecanism (două bifurcații șa-nod) explică schimbările bruște și greu reversibile din unele ecosisteme (Exercițiul 5).

### ✅ Verificare rapidă
Pentru ecuația ``x' = \mu - x^2`` găsiți echilibrele și stabilitatea lor în funcție de ``\mu``, apoi desenați diagrama de bifurcație. Ce tip de bifurcație are loc și pentru ce valoare a lui ``\mu``?
"""

# ╔═╡ 94c33587-fa82-4d98-bb41-de1070af33ba
hint(md"Pentru ``\mu < 0`` nu există echilibre, iar ``x`` scade la nesfîrșit. Pentru ``\mu = 0`` există un singur echilibru, ``x^* = 0``, semistabil. Pentru ``\mu > 0``: ``x^* = \pm\sqrt{\mu}``; cum ``f'(x) = -2x``, ``\sqrt{\mu}`` este stabil, iar ``-\sqrt{\mu}`` este instabil. Diagrama este parabola ``\mu = x^2``, culcată, cu ramura de sus continuă și cea de jos punctată: o **bifurcație șa-nod** în ``\mu = 0``. Aceasta este forma cea mai simplă a bifurcației șa-nod. Lîngă ``x = 50``, modelul cu cotă fixă arată la fel: notînd ``y = x - 50``, obținem ``y' = (12.5 - H) - 0.005\,y^2``.")

# ╔═╡ b0710d3e-1894-495c-9f37-7cfa1feeb686
md"""
## 💫 4. Bifurcația Hopf: nașterea unui ciclu limită

Bifurcațiile din secțiunea 3 privesc doar **echilibrele**. În sistemele cu cel puțin două variabile apare un tip nou de bifurcație: un echilibru stabil își pierde stabilitatea, iar **în jurul lui se naște un ciclu limită**. Aceasta este **bifurcația Hopf** (după matematicianul E. Hopf, 1942). Am întîlnit-o deja în modelul hormonal, cînd ``n`` a depășit aproximativ ``8.1``.

### 4.1. Cît de tare suflă clarinetistul?
Introducem în modelul clarinetului un parametru ``\mu``, care măsoară cîtă energie primește ancia de la aerul suflat de muzician:
```math
\begin{aligned}
x' &= v \\
v' &= -x + \mu\,v - v^3 .
\end{aligned}
```
* Lîngă echilibru termenul ``v^3`` este neglijabil și rămîne arcul cu frecare din Cursul 5, cu ``c = -\mu``. Prin urmare echilibrul este o **spirală stabilă** pentru ``\mu < 0`` și o **spirală instabilă** pentru ``\mu > 0``.
* Departe de echilibru termenul ``-v^3`` frînează mișcarea, deci oscilațiile nu pot crește la nesfîrșit.

### 🎛️ Experimentați
Traiectoria mov pornește din interior, aproape de echilibru, iar cea portocalie din exterior.

Energia primită de ancie ``\mu``: $(@bind μ_cl PlutoUI.Slider(-1.0:0.05:2.0, show_value=true, default=0.5))
"""

# ╔═╡ 0c280fcc-f6af-4789-abf6-c07ece80a1cd
# amplitudinea ciclului limită al clarinetului (maximul lui x) pentru 0 < μ ≤ 2
amplitudini_clarinet = let
	μs = 0.05:0.05:2.0
	amp = map(μs) do μ
		F(u) = (u[2], -u[1] + μ * u[2] - u[2]^3)
		t, (X, _) = solutie(F, (0.1, 0.0), 300.0; h = 0.02)
		maximum(X[t .> 200])
	end
	(; μs, amp)
end

# ╔═╡ bbf70a10-828d-4cc2-88af-9818e447263d
let
	F(u) = (u[2], -u[1] + μ_cl * u[2] - u[2]^3)
	T = 60.0
	t, (X1, V1) = solutie(F, (0.1, 0.0), T)
	_, (X2, V2) = solutie(F, (2.2, 0.0), T)
	lo, hi, per = oscilatie(t, X2, 40)
	sub = μ_cl <= 0 ? "echilibrul este stabil: oscilațiile se sting, nu se aude nimic" :
		@sprintf("ciclu limită: amplitudinea ≈ %.2f, perioada ≈ %.2f", hi, per)

	fig = Figure(size = (1200, 520))
	ax1 = Axis(fig[1, 1], title = "Clarinetul (μ = $(μ_cl))", subtitle = sub,
		xlabel = "Deplasarea ancii (x)", ylabel = "Viteza (v)", aspect = 1,
		limits = (-2.5, 2.5, -2.5, 2.5))
	cimp!(ax1, F, (-2.5, 2.5), (-2.5, 2.5))
	lines!(ax1, X1, V1, color = :purple, linewidth = 2.5)
	lines!(ax1, X2, V2, color = :darkorange, linewidth = 2.5)
	scatter!(ax1, [0.1, 2.2], [0.0, 0.0], color = [:purple, :darkorange], marker = :star5,
		markersize = 16)
	echilibru!(ax1, 0, 0, μ_cl <= 0)

	(; μs, amp) = amplitudini_clarinet
	ax2 = Axis(fig[1, 2], title = "Diagrama de bifurcație (Hopf)",
		subtitle = "echilibrul (negru) și extremele lui x pe ciclul limită (mov)",
		xlabel = "Parametrul μ", ylabel = "x", limits = (-1, 2, -2.3, 2.3))
	lines!(ax2, [-1, 0], [0, 0], color = :black, linewidth = 3)
	lines!(ax2, [0, 2], [0, 0], color = :black, linewidth = 2, linestyle = :dash)
	lines!(ax2, [0; μs], [0; amp], color = :purple, linewidth = 3, label = "ciclul limită")
	lines!(ax2, [0; μs], -[0; amp], color = :purple, linewidth = 3)
	μf = range(0, 2, length = 200)
	lines!(ax2, μf, 2 .* sqrt.(μf ./ 3), color = :gray50, linestyle = :dot, linewidth = 2,
		label = "aproximația 2√(μ/3)")
	lines!(ax2, μf, -2 .* sqrt.(μf ./ 3), color = :gray50, linestyle = :dot, linewidth = 2)
	vlines!(ax2, μ_cl, color = :darkorange, linewidth = 2)
	axislegend(ax2, position = :lt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ 9f93e844-d1a9-4e6c-9080-829b63067c6d
md"""
**Ce observăm:**
* Pentru ``\mu \le 0`` ambele traiectorii ajung la echilibru: ancia se oprește și nu se aude nimic.
* Pentru ``\mu > 0`` echilibrul devine instabil și se naște un ciclu limită. Amplitudinea lui pornește de la ``0`` și crește **continuu**, aproximativ ca ``2\sqrt{\mu/3}``: un clarinetist care suflă puțin peste prag produce un sunet foarte slab, care se întărește treptat. Aceasta este **bifurcația Hopf supercritică**. Există și varianta **subcritică**, în care oscilațiile apar brusc, cu amplitudine mare (Cursul 10).
* Lîngă bifurcație perioada este aproape ``2\pi``: **înălțimea** notei este dată de partea liniară a modelului (arcul), iar **intensitatea** ei de cît de tare se suflă.

În cursurile 8–10 vom vedea cum se recunoaște bifurcația Hopf din **valorile proprii** ale liniarizării: o pereche de valori proprii complexe își schimbă partea reală din negativă în pozitivă.

### 4.2. Paradoxul îmbogățirii
Modelul Lotka–Volterra are două defecte: iepurii fără rîși cresc exponențial, iar un rîs poate mînca oricît de mulți iepuri. Un model mai realist (M. Rosenzweig și R. MacArthur, 1963) corectează ambele defecte:
```math
\begin{aligned}
I' &= I\left(1 - \frac{I}{K}\right) - \frac{I}{1 + I}\,R \\
R' &= 0.5\,\frac{I}{1 + I}\,R - 0.2\,R .
\end{aligned}
```
* Fără rîși, iepurii cresc **logistic**, cu capacitatea de suport ``K``, care măsoară cît de bogată este vegetația.
* Un rîs mănîncă ``I/(1 + I)`` iepuri pe unitatea de timp: funcția Michaelis–Menten din Cursul 1, numită în ecologie **răspuns funcțional de tip II** (C. S. Holling, 1959). Cînd iepurii sînt mulți, rîsul este sătul și nu mai poate mînca mai mult.
* Rîșii transformă hrana în pui (coeficientul ``0.5``) și mor cu rata ``0.2``.

**Null-clinele** (Cursul 5):
* ``R' = 0 \iff R = 0`` sau ``0.5\,\dfrac{I}{1 + I} = 0.2 \iff I = \dfrac{2}{3}``: o dreaptă **verticală**;
* ``I' = 0 \iff I = 0`` sau ``R = \left(1 - \dfrac{I}{K}\right)(1 + I)``: o parabolă în formă de **cocoașă**, cu vîrful în ``I = \dfrac{K - 1}{2}``.

Echilibrul de coexistență este intersecția dintre dreapta verticală și cocoașă; el există pentru ``K > 2/3``.

> 📐 **Regula vîrfului** (Rosenzweig): echilibrul de coexistență este **stabil** dacă dreapta verticală taie cocoașa **la dreapta vîrfului** și **instabil**, înconjurat de un ciclu limită, dacă o taie **la stînga vîrfului**. Bifurcația Hopf are loc cînd dreapta trece prin vîrf: ``\dfrac{K - 1}{2} = \dfrac{2}{3}``, adică ``K = \dfrac{7}{3} \approx 2.33``.

### 🎛️ Experimentați
Capacitatea de suport a iepurilor ``K``: $(@bind K_rm PlutoUI.Slider(0.5:0.05:6.0, show_value=true, default=2.0))
"""

# ╔═╡ d993ab39-7b69-442e-8491-8f78b30daf20
# extremele numărului de rîși pe ciclul limită, pentru K > 7/3
cicluri_rm = let
	Ks = 2.35:0.05:6.0
	extreme = map(Ks) do K
		F(u) = (u[1] * (1 - u[1] / K) - u[1] * u[2] / (1 + u[1]), 0.5 * u[1] * u[2] / (1 + u[1]) - 0.2 * u[2])
		t, (_, R) = solutie(F, (1.0, 0.5), 1500.0; h = 0.1)
		extrema(R[t .> 1000])
	end
	(; Ks, Rmin = first.(extreme), Rmax = last.(extreme))
end

# ╔═╡ 2ccbaf4a-ef8c-4688-add0-33369b2bd0df
let
	K = K_rm
	Is = 2 / 3                         # null-clina verticală a rîșilor
	Rs(K) = (1 - Is / K) * (1 + Is)    # rîșii la echilibrul de coexistență
	F(u) = (u[1] * (1 - u[1] / K) - u[1] * u[2] / (1 + u[1]), 0.5 * u[1] * u[2] / (1 + u[1]) - 0.2 * u[2])
	T = 400.0
	t, (I, R) = solutie(F, (1.0, 0.5), T; h = 0.05)
	regim = K <= Is ? "prea puțină hrană: rîșii dispar" :
		K <= 7 / 3 ? "coexistență la echilibru (stabil)" :
		@sprintf("ciclu limită: iepurii coboară pînă la %.3f", minimum(I[t .> 200]))

	fig = Figure(size = (1350, 480))
	xmax, ymax = max(2.0, 1.1K), max(2.0, 1.1maximum(R))
	ax1 = Axis(fig[1, 1], title = "Null-cline și cîmpul de vectori (K = $(K))", subtitle = regim,
		xlabel = "Iepuri (I)", ylabel = "Rîși (R)", limits = (-0.02xmax, xmax, -0.02ymax, ymax))
	cimp!(ax1, F, (0, xmax), (0, ymax))
	Ig = range(0, K, length = 200)
	lines!(ax1, Ig, (1 .- Ig ./ K) .* (1 .+ Ig), color = :darkorange, linewidth = 3, label = "I' = 0")
	vlines!(ax1, 0, color = :darkorange, linewidth = 3)
	vlines!(ax1, Is, color = :teal, linewidth = 3, label = "R' = 0")
	hlines!(ax1, 0, color = :teal, linewidth = 3)
	if K > 1
		scatter!(ax1, [(K - 1) / 2], [(K + 1)^2 / (4K)], marker = :diamond, color = :darkorange,
			markersize = 14)
	end
	lines!(ax1, I, R, color = :crimson, linewidth = 2.5)
	scatter!(ax1, [1.0], [0.5], color = :crimson, marker = :star5, markersize = 18)
	echilibru!(ax1, 0, 0, false)
	echilibru!(ax1, K, 0, K < Is)
	K > Is && echilibru!(ax1, Is, Rs(K), K < 7 / 3)
	axislegend(ax1, position = :rt, backgroundcolor = (:white, 0.85))

	ax2 = Axis(fig[1, 2], title = "Evoluția în timp", xlabel = "Timp (t)", ylabel = "Populație",
		limits = (0, T, 0, nothing))
	lines!(ax2, t, I, color = :royalblue, linewidth = 2.5, label = "iepuri I(t)")
	lines!(ax2, t, R, color = :crimson, linewidth = 2.5, label = "rîși R(t)")
	axislegend(ax2, position = :rt, backgroundcolor = (:white, 0.85))

	(; Ks, Rmin, Rmax) = cicluri_rm
	ax3 = Axis(fig[1, 3], title = "Diagrama de bifurcație pentru rîși",
		subtitle = "transcritică în K = 2/3, Hopf în K = 7/3", xlabel = "Capacitatea de suport K",
		ylabel = "Rîși (R)", limits = (0.5, 6, -0.1, 4.2))
	lines!(ax3, [0.5, Is], [0, 0], color = :black, linewidth = 3)
	lines!(ax3, [Is, 6], [0, 0], color = :black, linewidth = 2, linestyle = :dash)
	K1, K2 = range(Is, 7 / 3, length = 100), range(7 / 3, 6, length = 100)
	lines!(ax3, K1, Rs.(K1), color = :black, linewidth = 3)
	lines!(ax3, K2, Rs.(K2), color = :black, linewidth = 2, linestyle = :dash)
	lines!(ax3, [7 / 3; Ks], [Rs(7 / 3); Rmax], color = :purple, linewidth = 3)
	lines!(ax3, [7 / 3; Ks], [Rs(7 / 3); Rmin], color = :purple, linewidth = 3)
	vlines!(ax3, K, color = :darkorange, linewidth = 2)
	fig
end

# ╔═╡ 7f4f3d54-78a9-4259-bea9-4522bd74618a
md"""
**Ce observăm:**
* ``K < 2/3``: chiar cu toți iepurii pe care îi poate susține vegetația, rîșii nu au destulă hrană și dispar. La ``K = 2/3`` echilibrul de coexistență se desprinde din axa ``R = 0``: o **bifurcație transcritică**, ca la pescuit. De acum rîșii pot invada.
* ``2/3 < K < 7/3``: vîrful cocoașei este la stînga dreptei verticale, iar populațiile se stabilizează (spiralînd) la echilibrul de coexistență.
* ``K > 7/3``: vîrful trece la dreapta, echilibrul devine instabil și apare un **ciclu limită** (bifurcație Hopf). Spre deosebire de Lotka–Volterra, oscilațiile nu depind de starea inițială: un astfel de model explică mai bine ciclurile regulate de aproximativ ``10`` ani ale rîșilor și iepurilor (Cursul 1).
* Cu cît ``K`` este mai mare, cu atît oscilațiile sînt mai ample, iar în fiecare ciclu iepurii (și apoi rîșii) ajung aproape de zero. În realitate, la numere atît de mici o întîmplare (o iarnă grea, o boală) poate duce la **dispariția** speciei.

> 🌾 **Paradoxul îmbogățirii** (M. Rosenzweig, 1971): îmbogățirea mediului (fertilizarea pășunii, hrana suplimentară) poate **destabiliza** un ecosistem în loc să-l ajute.

### ✅ Verificare rapidă
Un fermier fertilizează pășunea, iar capacitatea de suport a iepurilor crește de la ``K = 2`` la ``K = 4``. Folosind regula vîrfului, prevedeți ce se întîmplă. Care este cea mai mare valoare a lui ``K`` pentru care cele două specii coexistă la echilibru?
"""

# ╔═╡ 00613d90-c130-4b6c-996e-5aaf9ca125cd
hint(md"Pentru ``K = 2`` vîrful cocoașei este la ``I = 0.5 < 2/3``: dreapta verticală taie cocoașa la dreapta vîrfului, deci echilibrul ``(2/3,\ 10/9)`` este stabil. Pentru ``K = 4`` vîrful este la ``I = 1.5 > 2/3``: echilibrul devine instabil, iar populațiile oscilează pe un ciclu limită (verificați cu slider-ul). Valoarea maximă: ``(K - 1)/2 = 2/3``, adică ``K = 7/3 \approx 2.33``.")

# ╔═╡ 344d0337-37a3-4571-a7e3-36f83e480844
md"""
## 📝 5. Recapitulare și exerciții

### 🧠 Recapitulare
- 🔄 Un **centru** (Lotka–Volterra) are o familie de curbe închise, iar amplitudinea oscilațiilor depinde de starea inițială: comportament fragil. Un **ciclu limită** este o curbă închisă izolată; cînd este stabil, el atrage traiectoriile vecine: oscilații robuste, cu amplitudine și perioadă proprii.
- 🎷 **Primul mecanism:** energia este pompată în sistem la amplitudini mici (frecare negativă) și disipată la amplitudini mari (clarinetul lui Rayleigh).
- 🧬 **Al doilea mecanism:** **feedback negativ cu întîrziere** (hormonii, ceasul Hes1). Cu cît bucla are mai multe etape, cu atît oscilațiile apar mai ușor; cu o singură etapă (o variabilă) sau cu două etape nu apar oscilații susținute.
- 🦋 O **bifurcație** este o schimbare calitativă a comportamentului unui sistem cînd un parametru trece printr-o valoare critică; **diagrama de bifurcație** arată echilibrele (și ciclurile) în funcție de parametru.

| Bifurcația | Ce se întîmplă | Exemple din curs |
| :--- | :--- | :--- |
| **transcritică** | două echilibre se intersectează și își schimbă stabilitatea | pescuit cu efort constant (``E = 0.2``); invazia rîșilor (``K = 2/3``) |
| **șa-nod** | un echilibru stabil și unul instabil se ciocnesc și dispar | pescuit cu cotă fixă (``H = 12.5``); comutatorul lac |
| **Hopf** | un echilibru stabil devine instabil și se naște un ciclu limită | modelul hormonal (``n \approx 8.1``); clarinetul (``\mu = 0``); paradoxul îmbogățirii (``K = 7/3``) |

- 🔁 Două bifurcații șa-nod (diagrama în formă de S) produc **bistabilitate** și **histerezis**: salturi bruște și memorie (comutatorul lac).
- ⚠️ Tipul bifurcației decide dacă o schimbare a mediului produce un declin **treptat** și reversibil sau o prăbușire **bruscă** și greu reversibilă.

---
### ✍️ Exerciții
**Exercițiul 1 (captura maximă durabilă).** O populație crește logistic și este pescuită cu efort constant ``E``:
```math
x' = r\,x\left(1 - \frac{x}{K}\right) - E\,x .
```
a) Exprimați echilibrul nenul ``x^*(E)`` și captura la echilibru ``Y(E) = E\,x^*(E)``.

b) Găsiți efortul ``E`` pentru care captura este maximă și valoarea acestei capturi maxime.

c) Pentru ``r = 0.5`` și ``K = 100``, comparați captura maximă cu cota critică ``H = 12.5`` din secțiunea 3.2. De ce este periculos să se fixeze o cotă egală cu captura maximă?
"""

# ╔═╡ c3a69d29-5d2c-4c46-b372-3d205726705d
hint(md"a) ``x^* = K(1 - E/r)``, deci ``Y(E) = EK(1 - E/r)``. b) ``Y'(E) = K(1 - 2E/r) = 0`` pentru ``E = r/2``; atunci ``x^* = K/2`` și ``Y_{max} = rK/4``. c) ``rK/4 = 12.5``: cota egală cu captura maximă este **exact** punctul de bifurcație șa-nod. O mică supraestimare a lui ``r`` sau ``K``, sau un an prost, face ca ``x' < 0`` peste tot, iar stocul se prăbușește. Cu efortul ``E = r/2`` echilibrul rămîne stabil (``f'(x^*) = -r/2``), iar un an prost reduce stocul doar temporar.")

# ╔═╡ 2893cf23-e3f4-441d-b48e-f8513088facc
md"""
**Exercițiul 2 (bifurcația transcritică).**

a) Pentru ``x' = \mu x - x^2``, găsiți echilibrele și stabilitatea lor în funcție de ``\mu`` și desenați diagrama de bifurcație.

b) Arătați că modelul pescuitului cu efort constant, ``x' = r\,x\left(1 - \dfrac{x}{K}\right) - E\,x``, are aceeași formă. Cine joacă rolul lui ``\mu``?

c) În modelul cu efect Allee din Cursul 3, ``x' = r\,x\left(1 - \dfrac{x}{K}\right)\left(\dfrac{x}{a} - 1\right)``, ce se întîmplă cu echilibrele ``a`` și ``K`` cînd pragul ``a`` crește și trece prin ``K``?
"""

# ╔═╡ fe319fa8-862d-4e2c-a6c1-5d5450c569cb
hint(md"a) ``x^* = 0`` și ``x^* = \mu``; ``f'(x) = \mu - 2x``, deci ``f'(0) = \mu`` și ``f'(\mu) = -\mu``. Pentru ``\mu < 0``: ``0`` stabil, ``\mu`` instabil; pentru ``\mu > 0`` invers. Diagrama: dreptele ``x = 0`` și ``x = \mu``, care se intersectează în origine și își schimbă acolo stabilitatea (continuu ↔ punctat). b) ``x' = (r - E)\,x - \frac{r}{K}\,x^2``: rolul lui ``\mu`` îl joacă ``r - E``, iar bifurcația are loc pentru ``E = r``. c) ``f'(a) = r(1 - a/K)`` și ``f'(K) = r(1 - K/a)``. Pentru ``a < K``, ``a`` este instabil și ``K`` stabil; pentru ``a > K`` invers. La ``a = K`` echilibrele se ciocnesc și își schimbă stabilitatea: o bifurcație transcritică (biologic, cazul ``a > K`` nu are sens).")

# ╔═╡ b0201582-10bc-41bb-9421-23bf239f5609
md"""
**Exercițiul 3 (oscilații glicolitice).** În celulele de drojdie, concentrațiile intermediarilor glicolizei oscilează cu o perioadă de cîteva minute. Enzima PFK transformă fructoza-6-fosfat (F6P) în ADP și este **activată de ADP**: un feedback pozitiv. Modelul lui E. Sel'kov (1968), cu ``x`` = ADP și ``y`` = F6P, este
```math
\begin{aligned}
x' &= -x + a\,y + x^2 y \\
y' &= b - a\,y - x^2 y ,
\end{aligned}
\qquad a = 0.08 ,
```
unde ``b`` este rata cu care este furnizată F6P.

a) Arătați că echilibrul este ``x^* = b``, ``y^* = \dfrac{b}{a + b^2}``. (Indicație: adunați cele două ecuații.)

b) Rulați codul de mai jos pentru ``b = 0.2``, ``b = 0.6`` și ``b = 1``. Pentru care valori apare un ciclu limită?

c) Estimați intervalul de valori ale lui ``b`` pentru care glicoliza oscilează. Cîte bifurcații Hopf are modelul?
"""

# ╔═╡ 8d9323f5-4026-4fc8-a24b-aa17b7cac35f
hint(md"a) Adunînd ecuațiile, la echilibru ``b - x = 0``, deci ``x^* = b``; din a doua ecuație ``y^*(a + b^2) = b``. b) Pentru ``b = 0.2`` și ``b = 1`` oscilațiile se sting (echilibru stabil); pentru ``b = 0.6`` apare un ciclu limită. c) Aproximativ ``0.35 < b < 0.85``: modelul are **două** bifurcații Hopf, la capetele intervalului (exact, ``b^2 = \frac{1}{2}\left(1 - 2a \pm \sqrt{1 - 8a}\right)``, adică ``b \approx 0.346`` și ``b \approx 0.849``).")

# ╔═╡ e800ae38-a356-4488-9280-4f47e2663756
# EXERCIȚIUL 3 b): modificați b și rulați
let
	a, b = 0.08, 0.6
	F(u) = (-u[1] + a * u[2] + u[1]^2 * u[2], b - a * u[2] - u[1]^2 * u[2])
	t, (X, Y) = solutie(F, (0.5, 0.5), 100)
	fig = Figure(size = (1100, 420))
	ax1 = Axis(fig[1, 1], title = "Spațiul stărilor (b = $b)", xlabel = "ADP (x)", ylabel = "F6P (y)")
	lines!(ax1, X, Y, color = :teal, linewidth = 2)
	echilibru!(ax1, b, b / (a + b^2), nothing)   # completați: stabil (true) sau instabil (false)?
	ax2 = Axis(fig[1, 2], title = "Evoluția în timp", xlabel = "Timp (t)", ylabel = "Concentrație")
	lines!(ax2, t, X, color = :darkorange, linewidth = 2, label = "ADP x(t)")
	lines!(ax2, t, Y, color = :royalblue, linewidth = 2, label = "F6P y(t)")
	axislegend(ax2, position = :rt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ f27069f8-5f21-4121-b936-ab61edcd416f
md"""
**Exercițiul 4 (modelul Holling–Tanner).** O altă variantă a modelului prădător–pradă:
```math
\begin{aligned}
N' &= N\left(1 - \frac{N}{7}\right) - \frac{w\,N}{N + 1}\,P \\
P' &= 0.1\,P\left(1 - \frac{P}{N}\right) .
\end{aligned}
```
a) Interpretați termenii. Care este capacitatea de suport a prădătorilor și de ce depinde ea de ``N``? Ce măsoară parametrul ``w``?

b) Rulați codul de mai jos pentru mai multe valori ale lui ``w`` între ``0.5`` și ``2``. Găsiți aproximativ valoarea lui ``w`` la care are loc o bifurcație Hopf.

c) Interpretați biologic rezultatul: ce se întîmplă cînd prădătorii devin vînători mai eficienți?
"""

# ╔═╡ 898f5912-d2b9-47e3-9a14-074911103b8f
hint(md"a) Prada crește logistic (capacitatea de suport ``7``) și este consumată după un răspuns funcțional de tip II; ``w`` este numărul maxim de prăzi pe care le poate consuma un prădător pe unitatea de timp. Prădătorii cresc logistic, cu capacitatea de suport ``N``: fiecare prădător are nevoie de o pradă, deci numărul de prădători pe care îi poate susține mediul este proporțional cu prada disponibilă. b) Pentru ``w`` mai mic decît aproximativ ``0.9`` populațiile ajung la echilibru (cu oscilații amortizate); pentru ``w > 0.9`` apare un ciclu limită, de exemplu pentru ``w = 1`` prada oscilează între ``0.9`` și ``4.4``, cu perioada ``\approx 29``. c) Prădătorii mai eficienți destabilizează sistemul: apar oscilații ample, iar la minimul ciclului prada ajunge la valori mici, deci crește riscul de dispariție.")

# ╔═╡ d3b2bf23-66eb-4989-a3f9-233c899afbda
# EXERCIȚIUL 4 b): modificați w și rulați
let
	w = 0.5
	F(u) = (u[1] * (1 - u[1] / 7) - w * u[1] * u[2] / (u[1] + 1), 0.1 * u[2] * (1 - u[2] / u[1]))
	t, (N, P) = solutie(F, (2.0, 2.0), 400)
	fig = Figure(size = (1100, 420))
	ax1 = Axis(fig[1, 1], title = "Holling–Tanner (w = $w)", xlabel = "Pradă (N)", ylabel = "Prădători (P)")
	lines!(ax1, N, P, color = :teal, linewidth = 2)
	ax2 = Axis(fig[1, 2], title = "Evoluția în timp", xlabel = "Timp (t)", ylabel = "Densitate")
	lines!(ax2, t, N, color = :royalblue, linewidth = 2, label = "pradă N(t)")
	lines!(ax2, t, P, color = :crimson, linewidth = 2, label = "prădători P(t)")
	axislegend(ax2, position = :rt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ 6ae7fccb-569d-4ea6-9da6-9e757d6f83ad
md"""
**Exercițiul 5 (discuție: lacuri limpezi și lacuri tulburi).** 🏞️ Un lac puțin adînc poate fi în două stări: **limpede**, cu multe plante acvatice care fixează sedimentele și consumă nutrienții, sau **tulbure**, dominat de alge microscopice, care umbresc plantele de pe fund. Parametrul principal este cantitatea de nutrienți (fosfor) ``p`` care ajunge în lac, de exemplu din îngrășămintele de pe cîmpurile din jur. Se observă că:
* un lac limpede rămîne limpede pe măsură ce ``p`` crește, pînă la o valoare ``p_2``, după care devine **brusc** tulbure;
* pentru a-l face din nou limpede nu este suficient să coborîm ``p`` puțin sub ``p_2``: trebuie să coborîm pînă la o valoare ``p_1`` mult mai mică (M. Scheffer și colaboratorii, 2001).

a) Desenați calitativ diagrama de bifurcație (turbiditatea apei în funcție de ``p``). Ce tip de bifurcații au loc în ``p_1`` și ``p_2``?

b) Explicați de ce un lac tulbure nu redevine limpede pentru ``p_1 < p < p_2``. Ce altă intervenție, în afară de reducerea lui ``p``, ar putea readuce lacul în starea limpede?

c) Comparați cu comutatorul lac din secțiunea 3.3.
"""

# ╔═╡ e175a556-755e-40ec-877a-c95cf6e8f809
hint(md"a) O curbă în formă de S, ca pentru operonul lac: ramura de jos (lac limpede) este stabilă pentru ``p < p_2``, ramura de sus (lac tulbure) este stabilă pentru ``p > p_1``, iar ramura din mijloc (punctată) este instabilă și joacă rolul de prag. În ``p_1`` și ``p_2`` au loc bifurcații **șa-nod**. b) Pentru ``p_1 < p < p_2`` sistemul este **bistabil**, iar lacul tulbure se află în bazinul de atracție al stării tulburi: rămîne tulbure (**histerezis**). O perturbare mare poate muta sistemul peste prag, în bazinul stării limpezi, de exemplu prin **biomanipulare**: pescuirea peștilor care mănîncă zooplanctonul, astfel încît zooplanctonul să consume mai multe alge. c) Este exact aceeași structură: două bifurcații șa-nod, bistabilitate între ele, salturi bruște și memorie; rolul lui ``k`` îl joacă ``p`` (cu sensul inversat).")

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─329b5e95-328d-41d9-8d69-c41a0fc4b4b7
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000002
# ╟─bfbc972c-4703-4eb1-9405-cddb964d5bc2
# ╟─ad37b4e7-848d-41d4-b35b-995ad51c7e4d
# ╟─1e41598f-baa1-4e36-89bd-3e52818e5ec9
# ╟─7ff60562-f93a-4471-bda3-0bec3962c417
# ╟─672a0d81-5c6c-4906-a828-8c847dd4b80b
# ╟─39bf943b-72f1-42b4-ab31-a1a8840af5d9
# ╟─f1287f65-d265-45d6-8eb6-e8f1cae80852
# ╟─ecd0a0e6-78ee-435a-9c71-fb9c679a4d7a
# ╟─6cf2aef4-e72e-422d-b9d9-188d863172be
# ╟─866002fb-4d84-45ac-843d-22a888db5be7
# ╟─11a7f979-4026-486e-92c8-834077357a1f
# ╟─af3eceb0-6556-49d5-b0c7-acbfec2480ff
# ╟─4d4ae524-a518-46d9-952c-d27415be1d79
# ╟─aea3cccd-b20e-4b7b-b5d8-a76b71c1386b
# ╟─2a217a54-02cf-4500-a29e-9375a75eb4d9
# ╟─a7fb6f10-1998-4827-a6d4-9af07fbe67d1
# ╟─4cb83737-82f7-46f2-99e2-1f49222ada4a
# ╟─53cd73af-6c92-4bb4-97ac-ef4a3798e1f3
# ╟─30a7f680-64d8-4667-8288-80b0cbfe4e8d
# ╟─551bf8e5-e871-41d7-b24a-5baf75052a78
# ╟─b7f5c448-b5c7-4d19-9d0d-8897e5708ab0
# ╟─48ab4099-ca5b-47c6-a65b-6031fcb25efa
# ╟─25fdc3a9-caa3-49dc-89fc-b63ed025f11d
# ╟─907565a1-48e0-4481-981e-63bcfbd3b2ac
# ╟─6801d4b4-9f0f-452d-8fd1-84047b8a1b1b
# ╟─94c33587-fa82-4d98-bb41-de1070af33ba
# ╟─b0710d3e-1894-495c-9f37-7cfa1feeb686
# ╟─0c280fcc-f6af-4789-abf6-c07ece80a1cd
# ╟─bbf70a10-828d-4cc2-88af-9818e447263d
# ╟─9f93e844-d1a9-4e6c-9080-829b63067c6d
# ╟─d993ab39-7b69-442e-8491-8f78b30daf20
# ╟─2ccbaf4a-ef8c-4688-add0-33369b2bd0df
# ╟─7f4f3d54-78a9-4259-bea9-4522bd74618a
# ╟─00613d90-c130-4b6c-996e-5aaf9ca125cd
# ╟─344d0337-37a3-4571-a7e3-36f83e480844
# ╟─c3a69d29-5d2c-4c46-b372-3d205726705d
# ╟─2893cf23-e3f4-441d-b48e-f8513088facc
# ╟─fe319fa8-862d-4e2c-a6c1-5d5450c569cb
# ╟─b0201582-10bc-41bb-9421-23bf239f5609
# ╟─8d9323f5-4026-4fc8-a24b-aa17b7cac35f
# ╠═e800ae38-a356-4488-9280-4f47e2663756
# ╟─f27069f8-5f21-4121-b936-ab61edcd416f
# ╟─898f5912-d2b9-47e3-9a14-074911103b8f
# ╠═d3b2bf23-66eb-4989-a3f9-233c899afbda
# ╟─6ae7fccb-569d-4ea6-9da6-9e757d6f83ad
# ╟─e175a556-755e-40ec-877a-c95cf6e8f809
