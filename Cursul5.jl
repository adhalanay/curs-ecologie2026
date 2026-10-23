### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> title = "Cursul V"
#> description = "Sisteme dinamice în mai multe dimensiuni"
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

# ╔═╡ a1b2c3d4-e5f6-7890-abcd-ef1234567890
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using PlutoTeachingTools
	using Roots
	using LinearAlgebra
	using Printf

	# Setăm o temă vizuală curată și modernă pentru grafice
	CairoMakie.activate!(type = "svg")
	set_theme!(theme_minimal())

	TableOfContents()
end

# ╔═╡ daaa48f8-1afa-4828-8c41-aa2134dcd652
WidthOverDocs()

# ╔═╡ c10c3d78-56e0-4409-95d4-279cd27a0a08
begin
	# un pas al metodei Runge–Kutta de ordinul 4 pentru u' = F(u): o medie ponderată
	# a patru pante, mult mai precisă decît metoda Euler la același pas
	function rk4(F, u, h)
		k1 = F(u)
		k2 = F(u .+ h / 2 .* k1)
		k3 = F(u .+ h / 2 .* k2)
		k4 = F(u .+ h .* k3)
		u .+ h / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
	end

	# traiectoria sistemului u' = F(u) din plan pe [0, T] (pentru T < 0 mergem înapoi
	# în timp); ne oprim dacă traiectoria iese din pătratul [-lim, lim]²
	function traiectorie(F, u0, T; h = 0.01, lim = Inf)
		pas = sign(T) * h
		u = float.(Tuple(u0))
		X, Y, t = [u[1]], [u[2]], [0.0]
		for k in 1:round(Int, abs(T) / h)
			u = rk4(F, u, pas)
			push!(X, u[1]); push!(Y, u[2]); push!(t, k * pas)
			maximum(abs, u) > lim && break
		end
		X, Y, t
	end

	# metoda Euler în plan: n pași de lungime Δt, pornind din u0
	function euler2d(F, u0, Δt, n)
		X, Y = [float(u0[1])], [float(u0[2])]
		for _ in 1:n
			fx, fy = F((X[end], Y[end]))
			push!(X, X[end] + Δt * fx)
			push!(Y, Y[end] + Δt * fy)
		end
		X, Y
	end

	# cîmpul de vectori, desenat prin linii de curent (mai închise unde mișcarea este mai rapidă)
	cimp!(ax, F, xr, yr; densitate = 1.0) = streamplot!(ax, p -> Point2f(F(p)...),
		xr[1] .. xr[2], yr[1] .. yr[2], colormap = [:gray85, :gray35],
		linewidth = 1, arrow_size = 9, density = densitate)

	# punct de echilibru: negru dacă este stabil, alb dacă este instabil, gri pentru centru
	echilibru!(ax, x, y, stabil) = scatter!(ax, [x], [y], markersize = 16,
		strokewidth = 2.5, strokecolor = :black,
		color = isnothing(stabil) ? :gray60 : stabil ? :black : :white)

	# rîși (R) și iepuri (I): modelul Lotka–Volterra cu toți parametrii egali cu 1
	lotka_volterra(u) = (u[1] * u[2] - u[1], u[2] - u[1] * u[2])
	# mărime care rămîne constantă de-a lungul soluțiilor exacte
	H_lv(R, I) = R - log(R) + I - log(I)

	# competiția căprioare (D) – elani (M), cu coeficienții de competiție α și β
	competitie(α, β) = u -> (u[1] * (3 - u[1] - α * u[2]), u[2] * (2 - u[2] - β * u[1]))
	jacobian_comp(α, β, D, M) = [(3 - 2D - α * M) (-α * D); (-β * M) (2 - 2M - β * D)]

	# echilibrele modelului de competiție, ca tripletele (D, M, stabil)
	function echilibre_comp(α, β)
		E = [(0.0, 0.0), (3.0, 0.0), (0.0, 2.0)]
		d = 1 - α * β
		if abs(d) > 1e-9
			Ds, Ms = (3 - 2α) / d, (2 - 3β) / d
			Ds > 0 && Ms > 0 && push!(E, (Ds, Ms))
		end
		[(D, M, all(real.(eigvals(jacobian_comp(α, β, D, M))) .< 0)) for (D, M) in E]
	end

	# null-clinele modelului de competiție: D = 0, D + αM = 3 și M = 0, βD + M = 2
	function nullcline_comp!(ax, α, β, xmax)
		ds = range(0, xmax, length = 200)
		vlines!(ax, 0, color = :darkorange, linewidth = 3, label = "D-null-cline (D' = 0)")
		lines!(ax, ds, (3 .- ds) ./ α, color = :darkorange, linewidth = 3)
		hlines!(ax, 0, color = :teal, linewidth = 3, label = "M-null-cline (M' = 0)")
		lines!(ax, ds, 2 .- β .* ds, color = :teal, linewidth = 3)
	end
	nothing
end

# ╔═╡ b2c3d4e5-f6a7-8901-bcde-f12345678901
md"""
# 🌌 Cursul 5 — Sisteme dinamice în mai multe dimensiuni
"""

# ╔═╡ f0b5a93d-48f1-4539-99cd-7b78d49b8a6c
md"""
## 🎯 Obiectivele lecției
1. Să descriem un sistem cu două variabile de stare prin **cîmpul său de vectori** și să îi simulăm traiectoriile cu **metoda Euler în plan**.
2. Să găsim **punctele de echilibru** ale unui sistem plan și să recunoaștem tipurile lor: **nod**, **punct șa**, **spirală**, **centru**.
3. Să analizăm grafic un model de **competiție** cu ajutorul **null-clinelor** și să prezicem cine cîștigă competiția.
4. Să înțelegem **multistabilitatea**: **bazine de atracție**, **praguri** și **separatoare** (comutatorul genetic lac, efectele de prioritate).
5. Să trecem la **sisteme discrete în două dimensiuni**: populații structurate pe stadii de viață.
"""

# ╔═╡ 72d70b30-e0f6-4aa4-82aa-db9f5a09b43d
md"""
## 🔁 0. Ce știm deja
* **Cursul 2:** starea unui sistem este un punct din **spațiul stărilor**; o ecuație diferențială asociază fiecărei stări un vector de schimbare (**cîmpul de vectori**); traiectoriile se pot aproxima cu **metoda Euler**. Am scris deja modele cu două variabile: arcul ``(x, v)`` și sistemul Lotka–Volterra pentru rîși și iepuri.
* **Cursul 3:** pentru o singură variabilă, ``x' = f(x)``, echilibrele sînt soluțiile ecuației ``f(x^*) = 0``, iar stabilitatea se citește pe **linia de fază** sau din semnul lui ``f'(x^*)``.
* **Cursul 4:** în timp discret, ``x_{n+1} = f(x_n)``, echilibrele satisfac ``x^* = f(x^*)``, iar stabilitatea cere ``\lvert f'(x^*)\rvert < 1``.

**Ce se schimbă cînd avem două variabile?** Pe o dreaptă, o soluție a ecuației ``x' = f(x)`` se poate mișca doar într-un sens, monoton, și nu poate trece printr-un echilibru. În plan traiectoriile au **loc să ocolească**: se pot roti, pot oscila la nesfîrșit, se pot apropia de un echilibru pe o direcție și îndepărta pe alta. Tot ce vom face azi este să învățăm să „citim” aceste comportamente direct din ecuații, de cele mai multe ori **fără să le rezolvăm**.

> 💡 **Exemple din natură** cu mai multe variabile care interacționează:
> * 🐇 prădători și pradă (rîșii și iepurii de zăpadă din Cursul 1);
> * 🦌 specii care concurează pentru aceleași resurse;
> * 🧬 gene care se activează reciproc;
> * 🐻 populații împărțite pe vîrste (juvenili și adulți).
"""

# ╔═╡ c3d4e5f6-a7b8-9012-cdef-123456789012
md"""
## 🧭 1. Cîmpul de vectori și metoda Euler în plan

Un sistem cu două variabile de stare ``X`` și ``Y`` are forma
```math
\begin{aligned}
X' &= f(X, Y) \\
Y' &= g(X, Y) .
\end{aligned}
```
Viteza de schimbare a fiecărei variabile depinde de **ambele** variabile: aici apare interacțiunea. Starea sistemului este punctul ``(X, Y)`` din plan, iar ecuațiile asociază fiecărui punct vectorul ``\big(f(X,Y),\, g(X,Y)\big)``. O **traiectorie** este o curbă care, în fiecare punct al ei, este tangentă la vectorul din acel punct.

**Metoda Euler în plan.** Exact ca în Cursul 2, urmăm săgeata pe un interval scurt de timp ``\Delta t``, apoi recalculăm săgeata în noul punct:
```math
\begin{aligned}
X_{n+1} &= X_n + \Delta t\, f(X_n, Y_n) \\
Y_{n+1} &= Y_n + \Delta t\, g(X_n, Y_n) .
\end{aligned}
```
Atenție: ambele coordonate se actualizează **simultan**, folosind valorile vechi ``(X_n, Y_n)``.

**Exemplu: rîși și iepuri.** Reluăm modelul Lotka–Volterra din Cursul 2, cu toți parametrii egali cu ``1``. Notăm cu ``R`` rîșii și cu ``I`` iepurii (în unități convenabile):
```math
\begin{aligned}
R' &= R\,I - R \\
I' &= I - R\,I .
\end{aligned}
```
* iepurii se înmulțesc (termenul ``+I``) și sînt mîncați la întîlnirile cu rîșii (``-RI``);
* rîșii se hrănesc cu iepuri (``+RI``) și mor în lipsa hranei (``-R``).

De exemplu, în starea ``(R, I) = (2,\ 0.5)`` avem ``R' = 1 - 2 = -1`` și ``I' = 0.5 - 1 = -0.5``: săgeata arată în jos și spre stînga, adică **ambele** populații scad.

### 🎛️ Experimentați
Alegeți starea inițială și pasul ``\Delta t``. Curba mov este soluția calculată foarte precis (cu metoda Runge–Kutta, o variantă îmbunătățită a metodei Euler), iar punctele portocalii sînt pașii metodei Euler.
"""

# ╔═╡ c661a3e2-1f7a-4875-88ea-7fb972bb428b
md"""
Rîși la început ``R_0``: $(@bind R0_lv PlutoUI.Slider(0.2:0.1:3.0, show_value=true, default=2.0))

Iepuri la început ``I_0``: $(@bind I0_lv PlutoUI.Slider(0.2:0.1:3.0, show_value=true, default=0.5))

Pasul metodei Euler ``\Delta t``: $(@bind Δt_lv PlutoUI.Slider(0.01:0.01:0.5, show_value=true, default=0.1))
"""

# ╔═╡ 56324280-d312-4cf4-8e3c-de172f219a8c
let
	F = lotka_volterra
	u0 = (R0_lv, I0_lv)
	T = 15.0
	X, Y, t = traiectorie(F, u0, T)
	n = round(Int, T / Δt_lv)
	XE, YE = euler2d(F, u0, Δt_lv, n)
	tE = Δt_lv .* (0:n)
	HE = XE[end] > 0 && YE[end] > 0 ? @sprintf("%.3f", H_lv(XE[end], YE[end])) : "nedefinit"

	fig = Figure(size = (1200, 500))
	ax1 = Axis(fig[1, 1], title = "Spațiul stărilor: rîși și iepuri",
		subtitle = @sprintf("H = R − ln R + I − ln I:  la început %.3f,  la sfîrșitul metodei Euler %s",
			H_lv(u0...), HE),
		xlabel = "Rîși (R)", ylabel = "Iepuri (I)", aspect = 1, limits = (0, 4.5, 0, 4.5))
	cimp!(ax1, F, (0, 4.5), (0, 4.5))
	lines!(ax1, X, Y, color = :purple, linewidth = 3, label = "soluția precisă")
	scatterlines!(ax1, XE, YE, color = :darkorange, markersize = 5, linewidth = 1,
		label = "Euler, Δt = $(Δt_lv)")
	echilibru!(ax1, 1, 1, nothing)
	echilibru!(ax1, 0, 0, false)
	scatter!(ax1, [u0[1]], [u0[2]], color = :crimson, marker = :star5, markersize = 18,
		label = "starea inițială")
	axislegend(ax1, position = :rt, backgroundcolor = (:white, 0.85))

	ax2 = Axis(fig[1, 2], title = "Evoluția în timp (linii: soluția precisă, puncte: Euler)",
		xlabel = "Timp (t)", ylabel = "Populație", limits = (0, T, 0, 5))
	lines!(ax2, t, X, color = :crimson, linewidth = 3, label = "rîși R(t)")
	lines!(ax2, t, Y, color = :royalblue, linewidth = 3, label = "iepuri I(t)")
	scatter!(ax2, tE, XE, color = (:crimson, 0.5), markersize = 5)
	scatter!(ax2, tE, YE, color = (:royalblue, 0.5), markersize = 5)
	axislegend(ax2, position = :rt)
	fig
end

# ╔═╡ c7625076-300c-4814-bf6a-eeccb0ab6047
md"""
### 🤔 Ce am observat?
* Soluția precisă este o **curbă închisă**: după un timp populațiile revin exact la valorile inițiale, apoi totul se repetă. Avem **oscilații** fără nicio influență externă, ca în datele companiei Hudson's Bay din Cursul 1.
* **Maximul iepurilor precede maximul rîșilor.** Cînd iepurii sînt mulți, rîșii au hrană și se înmulțesc; cînd rîșii sînt mulți, iepurii sînt decimați; apoi rîșii flămînzi mor, iar iepurii își revin.
* **Metoda Euler spiralează spre exterior.** La fiecare pas mergem pe tangentă, iar tangenta la o curbă închisă „iese” din curbă. Erorile se adună de la un pas la altul: cu ``\Delta t`` mic spirala se strînge, cu ``\Delta t`` mare aproximarea devine inutilizabilă.

> 🧮 **De unde știm că traiectoriile exacte sînt închise?** Mărimea ``H(R, I) = R - \ln R + I - \ln I`` este **constantă** de-a lungul oricărei soluții exacte. Verificați: derivînd după ``t`` cu regula lanțului, ``\frac{d}{dt}H = \left(1 - \frac{1}{R}\right)R' + \left(1 - \frac{1}{I}\right)I' = 0``. Traiectoriile sînt deci curbele de nivel ale lui ``H``, care sînt închise. Metoda Euler **nu** păstrează ``H`` (vedeți subtitlul graficului). Despre astfel de sisteme, numite **conservative**, vom vorbi în Cursul 10.

### ✅ Verificare rapidă
1. Calculați vectorul de schimbare în starea ``(R, I) = (0.5,\ 2)``. În ce direcție arată? Ce înseamnă biologic?
2. Faceți un pas Euler cu ``\Delta t = 0.1`` pornind din ``(2,\ 0.5)``.
"""

# ╔═╡ 07fa44dc-1187-4e5c-88bd-8fba4eb57442
hint(md"1. ``R' = 0.5 \cdot 2 - 0.5 = 0.5`` și ``I' = 2 - 0.5 \cdot 2 = 1``: săgeata arată în sus și spre dreapta. Rîșii sînt puțini și au multă hrană, deci cresc; iepurii sînt mulți, dar prădătorii puțini, deci cresc și ei. 2. ``(R_1, I_1) = (2, 0.5) + 0.1\,(-1, -0.5) = (1.9,\ 0.45)``.")

# ╔═╡ 1fead660-881c-4b32-90db-fd06d113ac2c
md"""
## ⚖️ 2. Puncte de echilibru în plan

Un **punct de echilibru** ``(X^*, Y^*)`` este o stare în care nimic nu se schimbă, deci **ambele** viteze sînt nule:
```math
f(X^*, Y^*) = 0 \qquad \text{și} \qquad g(X^*, Y^*) = 0 .
```
Pentru rîși și iepuri: ``R(I - 1) = 0`` și ``I(1 - R) = 0``. Din prima ecuație, ``R = 0`` sau ``I = 1``. Dacă ``R = 0``, a doua ecuație dă ``I = 0``; dacă ``I = 1``, a doua dă ``R = 1``. Avem deci două echilibre: ``(0, 0)``, fără niciun animal, și ``(1, 1)``, în care cele două specii coexistă în număr constant.

**Stabilitate.** Ca în Cursul 3, un echilibru este **stabil** (atractor) dacă traiectoriile care pornesc aproape de el se apropie de el și **instabil** dacă unele dintre ele se îndepărtează. În plan însă ne putem apropia (sau îndepărta) în multe feluri: direct sau rotindu-ne, din toate direcțiile sau doar din unele. Apare și o posibilitate nouă: traiectoriile din jurul echilibrului ``(1, 1)`` nici nu se apropie, nici nu se îndepărtează, ci se **rotesc** în jurul lui la nesfîrșit.

### 2.1. Galeria echilibrelor
Cele mai simple exemple sînt sistemele **liniare** ``X' = aX + bY``, ``Y' = cX + dY``, care au echilibrul în origine. Ele prezintă șase comportamente tipice:

| Tip | Ce fac traiectoriile din apropiere | Stabil? |
| :--- | :--- | :---: |
| **nod stabil** | se apropie direct, fără să se rotească | da |
| **nod instabil** | se îndepărtează direct | nu |
| **punct șa** | se apropie pe o direcție și se îndepărtează pe alta | nu |
| **spirală stabilă** | se apropie rotindu-se | da |
| **spirală instabilă** | se îndepărtează rotindu-se | nu |
| **centru** | se rotesc pe curbe închise | neutru |

În figura de mai jos, echilibrele stabile sînt negre, cele instabile albe, iar centrul gri; curba roșie este o traiectorie particulară.
"""

# ╔═╡ e23780cd-5bbe-4343-b87b-ebd49295b640
let
	tipuri = [
		("Nod stabil", "x' = −x,  y' = −2y", u -> (-u[1], -2u[2]), true, (1.8, 1.5)),
		("Nod instabil", "x' = x,  y' = 2y", u -> (u[1], 2u[2]), false, (0.1, 0.05)),
		("Punct șa", "x' = x,  y' = −y", u -> (u[1], -u[2]), false, (0.02, 1.8)),
		("Spirală stabilă", "x' = y,  y' = −x − 0.5y", u -> (u[2], -u[1] - 0.5u[2]), true, (1.8, 0.0)),
		("Spirală instabilă", "x' = y,  y' = −x + 0.5y", u -> (u[2], -u[1] + 0.5u[2]), false, (0.1, 0.0)),
		("Centru", "x' = y,  y' = −x", u -> (u[2], -u[1]), nothing, (1.5, 0.0)),
	]
	fig = Figure(size = (1200, 820))
	for (i, (nume, ecuatii, F, stabil, u0)) in enumerate(tipuri)
		ax = Axis(fig[(i - 1) ÷ 3 + 1, (i - 1) % 3 + 1], title = nume, subtitle = ecuatii,
			xlabel = "x", ylabel = "y", aspect = 1, limits = (-2, 2, -2, 2))
		cimp!(ax, F, (-2, 2), (-2, 2), densitate = 0.8)
		X, Y, _ = traiectorie(F, u0, 20, lim = 2.5)
		lines!(ax, X, Y, color = :crimson, linewidth = 3)
		scatter!(ax, [u0[1]], [u0[2]], color = :crimson, markersize = 10)
		echilibru!(ax, 0, 0, stabil)
	end
	fig
end

# ╔═╡ 9b6c3fd8-d1a1-483d-bf2d-45a7c2f72132
md"""
> 🐎 **Punctul șa nu are corespondent în 1D.** Aproape toate traiectoriile trec pe lîngă el și pleacă, dar există exact două traiectorii care intră în el (în exemplul de mai sus, cele două jumătăți ale axei ``Oy``). Ele formează o curbă care împarte planul în două: traiectoriile de o parte pleacă spre dreapta, cele de cealaltă parte spre stînga. Vom vedea în secțiunea 4 că această curbă joacă rolul unui **prag**.

Azi recunoaștem tipul unui echilibru din desen. În Cursul 8 vom vedea cum se citește tipul direct din coeficienții ``a, b, c, d`` (cu ajutorul **valorilor proprii**), iar în Cursul 9 cum se aplică această idee sistemelor neliniare.

### 2.2. Arcul cu frecare
Adăugăm modelului arcului din Cursul 2 o forță de frecare proporțională cu viteza:
```math
\begin{aligned}
x' &= v \\
v' &= -x - c\,v ,
\end{aligned}
```
unde ``c`` este coeficientul de frecare. Echilibrul este ``(0, 0)``: arcul în repaus. Modificați ``c`` și urmăriți cum se schimbă **tipul** echilibrului, deși ecuațiile își păstrează forma.

Coeficientul de frecare ``c``: $(@bind c_arc PlutoUI.Slider(-1.0:0.05:3.0, show_value=true, default=0.5))
"""

# ╔═╡ 7cb25802-be2e-40da-b4ef-08422d5b4976
let
	F = u -> (u[2], -u[1] - c_arc * u[2])
	X, Y, t = traiectorie(F, (2.0, 0.0), 40, lim = 3)
	tip = abs(c_arc) < 1e-9 ? "centru (oscilații care nu se sting)" :
		c_arc < 0 ? "spirală instabilă (oscilații tot mai mari)" :
		c_arc < 2 ? "spirală stabilă (oscilații amortizate)" :
		"nod stabil (revenire fără oscilații)"

	fig = Figure(size = (1200, 480))
	ax1 = Axis(fig[1, 1], title = "Arcul cu frecare, c = $(c_arc)", subtitle = "echilibrul (0, 0): " * tip,
		xlabel = "Poziția x", ylabel = "Viteza v", aspect = 1, limits = (-2.5, 2.5, -2.5, 2.5))
	cimp!(ax1, F, (-2.5, 2.5), (-2.5, 2.5))
	lines!(ax1, X, Y, color = :crimson, linewidth = 3)
	scatter!(ax1, [2.0], [0.0], color = :crimson, markersize = 12)
	echilibru!(ax1, 0, 0, abs(c_arc) < 1e-9 ? nothing : c_arc > 0)

	ax2 = Axis(fig[1, 2], title = "Poziția în timp", xlabel = "Timp (t)", ylabel = "x(t)",
		limits = (0, 40, -3, 3))
	hlines!(ax2, 0, color = :gray50, linestyle = :dash)
	lines!(ax2, t, X, color = :crimson, linewidth = 2.5)
	fig
end

# ╔═╡ d95276d8-527f-44e4-8035-670f3c041851
md"""
**Ce observăm:**
* ``c = 0`` (fără frecare): **centru**, arcul oscilează la nesfîrșit cu aceeași amplitudine.
* ``0 < c < 2`` (frecare mică): **spirală stabilă**, oscilațiile se sting treptat.
* ``c \ge 2`` (frecare mare): **nod stabil**, arcul revine în repaus fără să mai treacă de partea cealaltă. Așa trebuie să funcționeze amortizoarele unei mașini.
* ``c < 0`` („frecare negativă”, adică sistemul primește energie, ca un leagăn împins la fiecare trecere): **spirală instabilă**. Un mecanism de acest tip, care nu explodează însă la nesfîrșit, explică oscilațiile clarinetului din Cursul 6.

Valoarea ``c = 0`` separă comportamente calitativ diferite: o mică schimbare a parametrului transformă un echilibru stabil într-unul instabil. Astfel de valori se numesc **puncte de bifurcație** (Cursul 6).

### ✅ Verificare rapidă
Ce tip de echilibru este ``(0, 0)`` pentru rîși și iepuri? Priviți cîmpul de vectori din secțiunea 1 de-a lungul celor două axe.
"""

# ╔═╡ 2d098eae-6e9a-4ac7-a892-851afdea616a
hint(md"Un **punct șa**. Pe axa ``R = 0`` (fără rîși) avem ``I' = I``: iepurii cresc exponențial și se îndepărtează de origine. Pe axa ``I = 0`` (fără iepuri) avem ``R' = -R``: rîșii mor și starea se apropie de origine. Echilibrul ``(1, 1)`` este, după desen, un **centru**.")

# ╔═╡ a3b4c5d6-e7f8-9012-abcd-123456789012
md"""
## 🦌 3. Null-cline: competiția căprioare–elani

Căprioarele (``D``) și elanii (``M``) se hrănesc cu aceleași plante. Fiecare specie, singură, ar crește logistic (Cursul 3), dar prezența celeilalte îi reduce rata de creștere:
```math
\begin{aligned}
D' &= D\,(3 - D - \alpha M) \\
M' &= M\,(2 - M - \beta D) .
\end{aligned}
```
* Fără competitor, căprioarele ajung la capacitatea de suport ``3``, iar elanii la ``2``.
* ``\alpha`` măsoară cît de mult îi încurcă un elan pe căprioare, iar ``\beta`` cît de mult îi încurcă o căprioară pe elani. Elanul este mult mai mare și mănîncă mai mult, deci este natural ca ``\alpha > \beta``. Începem cu ``\alpha = 1`` și ``\beta = 0.5``.

Sistemul este **neliniar** și nu îl putem rezolva explicit. Putem însă afla aproape totul despre el cu o metodă grafică.

### 3.1. Metoda null-clinelor
* **D-null-clina** este mulțimea punctelor în care ``D' = 0``. Acolo săgețile sînt **verticale** (``D`` nu se schimbă, se schimbă doar ``M``).
* **M-null-clina** este mulțimea punctelor în care ``M' = 0``. Acolo săgețile sînt **orizontale**.
* **Echilibrele** sînt punctele în care **o D-null-clină intersectează o M-null-clină** (atunci ambele viteze sînt nule).

Pentru modelul nostru:
* ``D' = 0 \iff D = 0`` sau ``D + \alpha M = 3`` (două drepte);
* ``M' = 0 \iff M = 0`` sau ``\beta D + M = 2`` (alte două drepte).

Cu ``\alpha = 1`` și ``\beta = 0.5`` găsim patru echilibre:
* ``D = 0`` și ``M = 0``: ``(0, 0)``;
* ``D + M = 3`` și ``M = 0``: ``(3, 0)``, doar căprioare;
* ``D = 0`` și ``0.5D + M = 2``: ``(0, 2)``, doar elani;
* ``D + M = 3`` și ``0.5D + M = 2``: scăzînd ecuațiile, ``0.5D = 1``, deci ``(2, 1)``, **coexistență**.

Atenție: intersecția dreptei ``D = 0`` cu dreapta ``D + M = 3`` **nu** este un echilibru, deoarece ambele sînt D-null-cline.

### 3.2. Direcția săgeților
Null-clinele împart cadranul în regiuni. Ca pe linia de fază din Cursul 3, în fiecare regiune semnele lui ``D'`` și ``M'`` sînt constante, deci este suficient să le calculăm într-un **punct test**:

| Regiunea | Punct test | ``D'`` | ``M'`` | Direcția săgeților |
| :--- | :---: | :---: | :---: | :---: |
| sub ambele null-cline | ``(0.5,\ 0.5)`` | ``+`` | ``+`` | ↗ |
| deasupra ambelor null-cline | ``(3,\ 3)`` | ``-`` | ``-`` | ↙ |
| între ele, în dreapta lui ``(2, 1)`` | ``(2.8,\ 0.4)`` | ``-`` | ``+`` | ↖ |
| între ele, în stînga lui ``(2, 1)`` | ``(1,\ 1.9)`` | ``+`` | ``-`` | ↘ |

Toate săgețile arată, „rotindu-se”, spre ``(2, 1)``: echilibrul de coexistență este **stabil**. Pe axa ``M = 0`` săgețile duc spre ``(3, 0)``, dar imediat deasupra axei elanii cresc: ``(3, 0)`` este un **punct șa**. La fel ``(0, 2)``, iar ``(0, 0)`` este un **nod instabil**.

### 🎛️ Experimentați
Mutați starea inițială (steaua roșie), apoi modificați coeficienții de competiție. Încercați, de exemplu, ``(\alpha, \beta) = (1, 1)``, ``(2, 0.5)`` și ``(2, 1)``.
"""

# ╔═╡ 4922edb7-72ae-43f4-866a-e9b15faae600
md"""
Coeficientul ``\alpha`` (efectul elanilor asupra căprioarelor): $(@bind α_c PlutoUI.Slider(0.1:0.1:3.0, show_value=true, default=1.0))

Coeficientul ``\beta`` (efectul căprioarelor asupra elanilor): $(@bind β_c PlutoUI.Slider(0.0:0.05:1.5, show_value=true, default=0.5))

Căprioare la început ``D_0``: $(@bind D0_c PlutoUI.Slider(0.0:0.1:4.5, show_value=true, default=0.5))

Elani la început ``M_0``: $(@bind M0_c PlutoUI.Slider(0.0:0.1:4.5, show_value=true, default=3.5))
"""

# ╔═╡ b4c5d6e7-f8a9-0123-bcde-234567890123
let
	F = competitie(α_c, β_c)
	X, Y, t = traiectorie(F, (D0_c, M0_c), 30)
	rezultat = abs(α_c - 1.5) < 1e-9 ? "caz limită" :
		α_c < 1.5 && β_c < 2 / 3 ? "coexistență" :
		α_c < 1.5 ? "căprioarele elimină elanii" :
		β_c < 2 / 3 ? "elanii elimină căprioarele" :
		"bistabilitate (cîștigă specia avantajată la început)"

	fig = Figure(size = (1200, 540))
	ax1 = Axis(fig[1, 1], title = "Null-cline și cîmpul de vectori (α = $(α_c), β = $(β_c))",
		subtitle = "rezultatul competiției: " * rezultat,
		xlabel = "Căprioare (D)", ylabel = "Elani (M)", aspect = 1, limits = (-0.05, 4.5, -0.05, 4.5))
	cimp!(ax1, F, (0, 4.5), (0, 4.5))
	nullcline_comp!(ax1, α_c, β_c, 4.5)
	lines!(ax1, X, Y, color = :crimson, linewidth = 3)
	scatter!(ax1, [D0_c], [M0_c], color = :crimson, marker = :star5, markersize = 18)
	for (D, M, stabil) in echilibre_comp(α_c, β_c)
		echilibru!(ax1, D, M, stabil)
	end
	axislegend(ax1, position = :rt, backgroundcolor = (:white, 0.85))

	ax2 = Axis(fig[1, 2], title = "Evoluția în timp", xlabel = "Timp (t)", ylabel = "Populație",
		limits = (0, 30, 0, nothing))
	lines!(ax2, t, X, color = :darkorange, linewidth = 3, label = "căprioare D(t)")
	lines!(ax2, t, Y, color = :teal, linewidth = 3, label = "elani M(t)")
	axislegend(ax2, position = :rt)
	fig
end

# ╔═╡ 75136446-2e33-4c5c-b968-cd2c3a2cc47e
md"""
### 3.3. Cine cîștigă? Criteriul invaziei
O întrebare naturală în ecologie: **poate o specie rară să se instaleze** pe un teritoriu ocupat de cealaltă?

* Dacă pe teritoriu trăiesc doar căprioare, la capacitatea de suport ``D = 3``, și sosesc cîțiva elani (``M`` foarte mic), atunci ``M' \approx M\,(2 - 3\beta)``. Elanii **pot invada** dacă ``2 - 3\beta > 0``, adică ``\beta < 2/3``.
* Analog, căprioarele pot invada un teritoriu ocupat de elani (``M = 2``) dacă ``3 - 2\alpha > 0``, adică ``\alpha < 3/2``.

Rezultă patru scenarii:

| | ``\alpha < 3/2``: căprioarele pot invada | ``\alpha > 3/2``: căprioarele nu pot invada |
| :--- | :---: | :---: |
| ``\beta < 2/3``: elanii pot invada | 🤝 **coexistență**, echilibrul interior este stabil | **elanii cîștigă**: ``(0, 2)`` |
| ``\beta > 2/3``: elanii nu pot invada | **căprioarele cîștigă**: ``(3, 0)`` | ⚔️ **bistabilitate**: cîștigă specia avantajată la început |

Echilibrul de coexistență, cînd există, este
```math
D^* = \frac{3 - 2\alpha}{1 - \alpha\beta}, \qquad M^* = \frac{2 - 3\beta}{1 - \alpha\beta} .
```

> 🌿 **Principiul excluderii competitive** (G. F. Gause, 1934): două specii care folosesc exact aceleași resurse nu pot coexista la nesfîrșit. Coexistența cere ca fiecare specie să se limiteze **mai mult pe sine** decît pe cealaltă; altfel una dintre ele dispare. Gause a observat acest lucru în experimente cu două specii de parameci (*Paramecium aurelia* și *P. caudatum*) crescute în aceeași cultură.

### ✅ Verificare rapidă
Pentru ``\alpha = \beta = 1`` (cele două specii se încurcă la fel de mult), cine cîștigă? Cum arată null-clinele? Verificați cu slider-ele.
"""

# ╔═╡ 74f4f50b-a680-439d-9cbe-d0af495811fe
hint(md"Avem ``\alpha = 1 < 3/2`` și ``\beta = 1 > 2/3``, deci **căprioarele cîștigă**. Null-clinele ``D + M = 3`` și ``D + M = 2`` sînt **paralele**, deci nu există echilibru de coexistență (numitorul ``1 - \alpha\beta`` este zero). Între cele două drepte căprioarele cresc, iar elanii scad, așa că orice traiectorie cu ``D_0 > 0`` ajunge la ``(3, 0)``.")

# ╔═╡ d6e7f8a9-b0c1-2345-defa-456789012345
md"""
## 🔀 4. Multistabilitate și bazine de atracție

Uneori un sistem are **mai multe** echilibre stabile, iar starea finală depinde de unde pornim. Mulțimea stărilor inițiale din care traiectoriile ajung la un anumit atractor se numește **bazinul de atracție** al acestuia. Bazinele vecine sînt despărțite de **praguri**: în 1D un echilibru instabil, în 2D o curbă numită **separatoare**.

### 4.1. Un comutator genetic: operonul lac
Bacteria *E. coli* preferă glucoza. Dacă în mediu există doar lactoză, ea pornește un grup de gene, **operonul lac**, care produce enzimele necesare pentru a o folosi. Una dintre ele, **permeaza**, pompează lactoza în celulă, iar lactoza din celulă activează producția de permează: un **feedback pozitiv** (Cursul 1). Notînd cu ``x`` concentrația de permează (adimensionalizată), un model simplu este
```math
x' = \underbrace{\frac{a + x^2}{1 + x^2}}_{\text{producție}} - \underbrace{k\,x}_{\text{degradare}}, \qquad a = 0.01,\ k = 0.4 .
```
* ``a`` este o producție de bază foarte mică: gena „scapă” puțin chiar cînd este oprită;
* termenul ``x^2/(1+x^2)`` este o **funcție Hill** (Cursul 1): producția crește brusc după un prag, apoi se saturează;
* ``kx`` reprezintă degradarea permeazei și diluarea ei prin diviziunea celulei.

Este un sistem cu o singură variabilă, dar cu **trei** echilibre: punctele în care producția este egală cu degradarea, adică intersecțiile celor două grafice.
"""

# ╔═╡ e7f8a9b0-c1d2-3456-efab-567890123456
let
	a, k = 0.01, 0.4
	productie(x) = (a + x^2) / (1 + x^2)
	F(x) = productie(x) - k * x
	ech = find_zeros(F, 0, 3)
	prag = ech[2]

	fig = Figure(size = (1200, 480))
	ax1 = Axis(fig[1, 1], title = "Producție și degradare",
		subtitle = @sprintf("echilibre: x* ≈ %.3f,  %.3f,  %.3f", ech...),
		xlabel = "Concentrația de permează (x)", ylabel = "Rata", limits = (0, 2.6, -0.05, 1.1))
	xs = range(0, 2.6, length = 400)
	lines!(ax1, xs, productie.(xs), color = :royalblue, linewidth = 3, label = "producție (a + x²)/(1 + x²)")
	lines!(ax1, xs, k .* xs, color = :crimson, linewidth = 3, label = "degradare k x")
	# linia de fază pe axa Ox: săgețile arată sensul în care se mișcă x
	repere = [0.0; ech; 2.6]
	for (p, q) in zip(repere[1:end-1], repere[2:end])
		s = sign(F((p + q) / 2))
		L = 0.5 * (q - p)
		arrows2d!(ax1, [(p + q) / 2 - s * L / 2], [0.0], [s * L], [0.0], color = :black,
			shaftwidth = 3, tipwidth = 14, tiplength = 12)
	end
	for (i, e) in enumerate(ech)
		echilibru!(ax1, e, productie(e), i != 2)
		echilibru!(ax1, e, 0.0, i != 2)
	end
	axislegend(ax1, position = :lt)

	ax2 = Axis(fig[1, 2], title = "Soluții pentru diferite valori inițiale",
		xlabel = "Timp (t)", ylabel = "x(t)", limits = (0, 25, 0, 3))
	hlines!(ax2, prag, color = :black, linestyle = :dash, label = @sprintf("prag x* ≈ %.3f", prag))
	for x0 in [0.1, 0.3, 0.42, 0.45, 0.48, 0.52, 1.0, 2.8]
		xs_t = [x0]
		for _ in 1:500
			push!(xs_t, rk4(F, xs_t[end], 0.05))
		end
		lines!(ax2, range(0, 25, length = 501), xs_t, linewidth = 2.5,
			color = x0 < prag ? :royalblue : :darkorange)
	end
	axislegend(ax2, position = :rt)
	fig
end

# ╔═╡ f8a9b0c1-d2e3-4567-fabc-678901234567
md"""
**Interpretare biologică:**
* ``x^* \approx 0.027`` este **stabil**: operonul este „**oprit**”, celula produce foarte puțină permează.
* ``x^* \approx 2.008`` este **stabil**: operonul este „**pornit**”, permeaza este produsă din plin.
* ``x^* \approx 0.465`` este **instabil**: el este **pragul** dintre cele două bazine de atracție, ``[0,\ 0.465)`` (albastru) și ``(0.465,\ \infty)`` (portocaliu).

Celula se comportă ca un **întrerupător cu memorie**: un impuls temporar (de exemplu, o creștere de scurtă durată a lactozei din mediu) care împinge ``x`` peste prag aprinde operonul, iar acesta **rămîne pornit** și după ce impulsul a trecut. Ce se întîmplă cînd parametrul ``k`` variază (două echilibre se ciocnesc și dispar) vom vedea în Cursul 6, la **bifurcația șa-nod**.

### 4.2. Bazine de atracție în plan: competiția puternică
Revenim la căprioare și elani, cu **competiție puternică**: ``\alpha = 2`` și ``\beta = 1`` (cazul bistabil din tabelul secțiunii 3.3; este modelul pe care îl vom liniariza în Cursul 9). Echilibrul de coexistență ``(1, 1)`` este acum un **punct șa**, iar ``(3, 0)`` și ``(0, 2)`` sînt ambele stabile.

Pentru fiecare stare inițială dintr-o grilă de ``121 \times 121`` puncte simulăm sistemul și colorăm punctul după echilibrul la care ajunge. Linia neagră punctată este **separatoarea**: cele două traiectorii care intră în punctul șa.

Coeficientul ``\alpha``: $(@bind α_b PlutoUI.Slider(1.6:0.1:3.0, show_value=true, default=2.0))

Coeficientul ``\beta``: $(@bind β_b PlutoUI.Slider(0.7:0.05:1.5, show_value=true, default=1.0))
"""

# ╔═╡ 788dd5e1-9da2-4081-89b6-5c139d276b82
let
	F = competitie(α_b, β_b)
	grila = range(0, 4, length = 121)
	# 1 dacă traiectoria ajunge la (3, 0), 0 dacă ajunge la (0, 2)
	bazin = [begin
			u = (D, M)
			for _ in 1:800
				u = rk4(F, u, 0.05)
			end
			u[1] > u[2] ? 1.0 : 0.0
		end for D in grila, M in grila]
	p = sum(bazin) / length(bazin)

	# separatoarea: pornim foarte aproape de punctul șa, pe direcția pe care traiectoriile
	# intră în el, și mergem înapoi în timp
	ech = echilibre_comp(α_b, β_b)
	Ds, Ms, _ = ech[4]
	val, vec = eigen(jacobian_comp(α_b, β_b, Ds, Ms))
	v = Tuple(real.(vec[:, argmin(real.(val))]))
	separatoare = [traiectorie(F, (Ds, Ms) .+ s * 1e-4 .* v, -30, lim = 6) for s in (1, -1)]

	fig = Figure(size = (1100, 650))
	ax = Axis(fig[1, 1], title = "Bazinele de atracție (α = $(α_b), β = $(β_b))",
		subtitle = @sprintf("bazinul căprioarelor ocupă %.0f%% din pătrat, al elanilor %.0f%%",
			100p, 100(1 - p)),
		xlabel = "Căprioare (D)", ylabel = "Elani (M)", aspect = 1, limits = (0, 4, 0, 4))
	heatmap!(ax, grila, grila, bazin, colormap = [(:teal, 0.3), (:darkorange, 0.3)],
		colorrange = (0, 1))
	cimp!(ax, F, (0, 4), (0, 4))
	nullcline_comp!(ax, α_b, β_b, 4)
	for (X, Y, _) in separatoare
		lines!(ax, X, Y, color = :black, linewidth = 3, linestyle = :dash)
	end
	for (D, M, stabil) in ech
		echilibru!(ax, D, M, stabil)
	end
	Legend(fig[1, 2],
		[PolyElement(color = (:darkorange, 0.3)), PolyElement(color = (:teal, 0.3)),
			LineElement(color = :black, linestyle = :dash, linewidth = 3)],
		["bazinul lui (3, 0): rămîn doar căprioarele", "bazinul lui (0, 2): rămîn doar elanii",
			"separatoarea"],
		orientation = :vertical, framevisible = false)
	fig
end

# ╔═╡ f3e849c7-86fa-4012-bb33-98c4c45cad3b
md"""
**Ce observăm:**
* Planul este împărțit de separatoare în **două bazine**. Starea finală nu depinde de parametri (aceiași pentru toate punctele), ci de **istoria** sistemului: de unde a pornit.
* În ecologie acest fenomen se numește **efect de prioritate**: specia care ocupă prima un teritoriu poate rezista invaziei celeilalte, chiar dacă, pornind din alte condiții, ar fi pierdut. În experimentele clasice ale lui Thomas Park cu două specii de gîndaci de făină (*Tribolium*), în anumite condiții de temperatură și umiditate specia cîștigătoare depindea de numărul inițial de indivizi din fiecare specie.
* Mărind ``\alpha`` (elanii devin competitori mai puternici), bazinul elanilor se lărgește.
* Cu cît un echilibru stabil este mai departe de separatoare, cu atît este mai **rezistent**: o perturbare mică (o iarnă grea, o epidemie) nu schimbă rezultatul, dar una mare poate muta sistemul în celălalt bazin, iar efectul este **ireversibil**. Astfel de praguri explică de ce unele ecosisteme (de exemplu, un lac cu apă limpede care devine tulbure din cauza poluării cu nutrienți) nu revin la starea inițială nici după ce cauza perturbării a dispărut.

### ✅ Verificare rapidă
Pentru ``\alpha = 2`` și ``\beta = 1``, o pădure are ``D = 1.5`` căprioare și ``M = 1`` elani (în unitățile modelului). O iarnă grea omoară jumătate din căprioare, iar elanii nu sînt afectați. Ce se întîmplă pe termen lung înainte și după iarnă?
"""

# ╔═╡ 823027a6-23f4-40ab-b647-d969fd6362ab
hint(md"Înainte de iarnă, ``(1.5,\ 1)`` este sub separatoare, în bazinul portocaliu: ar fi rămas doar căprioarele. După iarnă, ``(0.75,\ 1)`` este deasupra separatoarei, în bazinul elanilor: căprioarele vor dispărea. Un eveniment de scurtă durată a schimbat definitiv soarta comunității.")

# ╔═╡ a9b0c1d2-e3f4-5678-abcd-789012345678
md"""
## 🔢 5. Sisteme discrete în două dimensiuni

Ca în Cursul 4, multe populații se numără o dată pe an. Dacă împărțim populația în două clase, obținem un sistem discret în două dimensiuni:
```math
\begin{aligned}
X_{n+1} &= F(X_n, Y_n) \\
Y_{n+1} &= G(X_n, Y_n) .
\end{aligned}
```
Un **echilibru** este acum o stare care se reproduce identic de la un an la altul: ``F(X^*, Y^*) = X^*`` și ``G(X^*, Y^*) = Y^*``.

### 5.1. Un model structurat pe stadii
Împărțim femelele unei populații în **juvenili** ``J`` și **adulți** ``A``. În fiecare an:
* o fracțiune ``a`` dintre juvenili supraviețuiește și rămîne juvenilă, iar o fracțiune ``m`` supraviețuiește și se maturizează (restul de ``1 - a - m`` mor);
* fiecare adult produce în medie ``f`` pui (femele): **fecunditatea**;
* o fracțiune ``s`` dintre adulți supraviețuiește.

```math
\begin{aligned}
J_{n+1} &= a\,J_n + f\,A_n \\
A_{n+1} &= m\,J_n + s\,A_n
\end{aligned}
\qquad\Longleftrightarrow\qquad
\begin{pmatrix} J_{n+1} \\ A_{n+1} \end{pmatrix} =
\begin{pmatrix} a & f \\ m & s \end{pmatrix}
\begin{pmatrix} J_n \\ A_n \end{pmatrix} .
```
Sistemul este **liniar**: populația din anul următor se obține înmulțind vectorul ``(J_n, A_n)`` cu o matrice. Valorile implicite ale slider-elor (``a = 0.65``, ``m = 0.25``, ``f = 0.5``, ``s = 0.9``) descriu **ursul negru american**, pe care îl vom studia în detaliu în Cursul 7. Pornim cu ``100`` de juvenili și ``10`` adulți.

Juvenili care rămîn juvenili ``a``: $(@bind a_st PlutoUI.Slider(0.0:0.05:1.0, show_value=true, default=0.65))

Juvenili care se maturizează ``m``: $(@bind m_st PlutoUI.Slider(0.0:0.05:1.0, show_value=true, default=0.25))

Fecunditatea adulților ``f``: $(@bind f_st PlutoUI.Slider(0.0:0.05:2.0, show_value=true, default=0.5))

Supraviețuirea adulților ``s``: $(@bind s_st PlutoUI.Slider(0.0:0.01:1.0, show_value=true, default=0.9))
"""

# ╔═╡ c1d2e3f4-a5b6-7890-cdef-901234567890
let
	L = [a_st f_st; m_st s_st]
	ani = 30
	stari = [[100.0, 10.0]]
	for _ in 1:ani
		push!(stari, L * stari[end])
	end
	J, A = first.(stari), last.(stari)
	N = max.(J .+ A, 1e-6)

	# factorul de creștere pe termen lung și distribuția stabilă pe stadii (Cursul 7)
	valori = eigvals(L)
	i = argmax(abs.(valori))
	λ = abs(valori[i])
	v = abs.(real.(eigvecs(L)[:, i]))
	v = v ./ sum(v)

	avertisment = a_st + m_st > 1 ? "Atenție: a + m > 1, adică mai mulți juvenili decît existau!" : ""
	fig = Figure(size = (1300, 470))
	ax1 = Axis(fig[1, 1], title = "Spațiul stărilor (primii 12 ani)", subtitle = avertisment,
		xlabel = "Juvenili (J)", ylabel = "Adulți (A)")
	n_vizibil = 12
	Jmax, Amax = maximum(J[1:n_vizibil+1]), maximum(A[1:n_vizibil+1])
	r = min(Jmax / max(v[1], 1e-9), Amax / max(v[2], 1e-9))
	lines!(ax1, [0, r * v[1]], [0, r * v[2]], color = :purple, linestyle = :dash, linewidth = 2,
		label = "direcția distribuției stabile")
	scatterlines!(ax1, J[1:n_vizibil+1], A[1:n_vizibil+1], color = :darkorange, markersize = 9)
	text!(ax1, J[1:4], A[1:4], text = ["n = 0", "1", "2", "3"], offset = (6, 4), fontsize = 13)
	xlims!(ax1, 0, nothing); ylims!(ax1, 0, nothing)
	axislegend(ax1, position = :lt)

	ax2 = Axis(fig[1, 2], title = "Populația totală Nₙ = Jₙ + Aₙ (scară logaritmică)",
		subtitle = @sprintf("pe termen lung Nₙ₊₁ ≈ λ Nₙ, cu λ ≈ %.3f", λ),
		xlabel = "Anul (n)", ylabel = "Nₙ", yscale = log10)
	scatterlines!(ax2, 0:ani, N, color = :darkorange, markersize = 7)

	ax3 = Axis(fig[1, 3], title = "Structura populației",
		subtitle = @sprintf("proporția juvenililor tinde la %.2f", v[1]),
		xlabel = "Anul (n)", limits = (0, ani, 0, max(1.6, 1.1λ)))
	hlines!(ax3, [λ], color = :crimson, linestyle = :dash)
	hlines!(ax3, [v[1]], color = :royalblue, linestyle = :dash)
	scatterlines!(ax3, 1:ani, N[2:end] ./ N[1:end-1], color = :crimson, markersize = 7,
		label = "factorul anual Nₙ₊₁ / Nₙ")
	scatterlines!(ax3, 0:ani, J ./ N, color = :royalblue, markersize = 7,
		label = "proporția juvenililor Jₙ / Nₙ")
	axislegend(ax3, position = :rb)
	fig
end

# ╔═╡ e5942d67-a5cd-4b6b-9eb4-5d206ea9f431
md"""
**Ce observăm:**
* După cîțiva ani de **regim tranzitoriu**, populația totală crește **geometric**, ca în Cursul 4: ``N_{n+1} \approx \lambda N_n``. Pentru ursul negru ``\lambda = 1.15``, adică o creștere de ``15\%`` pe an. Pe scara logaritmică, creșterea geometrică apare ca o dreaptă.
* **Structura** populației se stabilizează: proporția juvenililor tinde la o valoare fixă (``50\%`` pentru urs), **indiferent de populația inițială**. Aceasta este **distribuția stabilă pe stadii**. În spațiul stărilor, punctele se aliniază pe o dreaptă prin origine.
* Factorul ``\lambda`` și direcția distribuției stabile sînt o **valoare proprie** și un **vector propriu** ai matricei; le vom calcula în Cursul 7.
* Singurul echilibru al sistemului (în afara cazului ``\lambda = 1``) este ``(0, 0)``: el este stabil (populația se stinge) dacă ``\lambda < 1`` și instabil dacă ``\lambda > 1``.

### ✅ Verificare rapidă
Braconajul reduce supraviețuirea adulților ``s``. Sub ce valoare a lui ``s`` populația de urși începe să scadă (ceilalți parametri rămînînd neschimbați)? Găsiți răspunsul cu slider-ul, apoi calculați-l: la limită, populația are un echilibru nenul ``(J^*, A^*)``.
"""

# ╔═╡ abd99d36-5deb-42d4-a713-9cb7d605b72d
hint(md"Un echilibru nenul satisface ``J^* = aJ^* + fA^*`` și ``A^* = mJ^* + sA^*``, adică ``(1 - a)J^* = fA^*`` și ``(1 - s)A^* = mJ^*``. Înmulțind cele două egalități și simplificînd cu ``J^* A^*``: ``(1 - a)(1 - s) = fm``. Pentru urs: ``0.35\,(1 - s) = 0.125``, deci ``s \approx 0.64``. Dacă mai puțin de ``64\%`` dintre adulți supraviețuiesc de la un an la altul, populația scade.")

# ╔═╡ e3f4a5b6-c7d8-9012-efab-123456789012
md"""
## 📝 6. Recapitulare și exerciții

### 🧠 Recapitulare
- 🧭 Un sistem ``X' = f(X, Y)``, ``Y' = g(X, Y)`` definește un **cîmp de vectori** în plan; traiectoriile sînt tangente la el. **Metoda Euler** actualizează **simultan** ambele variabile, iar erorile ei se acumulează (spirala rîși–iepuri).
- ⚖️ **Echilibre:** ``f(X^*, Y^*) = g(X^*, Y^*) = 0``. Tipuri: **nod** (stabil sau instabil), **punct șa**, **spirală** (stabilă sau instabilă), **centru**. Același echilibru își poate schimba tipul cînd variază un parametru (arcul cu frecare).
- 📏 **Null-cline:** pe ``X' = 0`` săgețile sînt verticale, pe ``Y' = 0`` orizontale; echilibrele sînt intersecțiile dintre o ``X``-null-clină și o ``Y``-null-clină; în fiecare regiune dintre null-cline direcția săgeților se află cu un punct test.
- ⚔️ **Competiție:** coexistență, excludere competitivă sau bistabilitate; decidem cu **criteriul invaziei**.
- 🔀 **Multistabilitate:** fiecare atractor are un **bazin de atracție**; bazinele sînt separate de un echilibru instabil (1D) sau de **separatoarea** unui punct șa (2D). Exemple: comutatorul lac, efectele de prioritate.
- 🔢 **Sisteme discrete 2D:** modelele pe stadii ``\mathbf{x}_{n+1} = L\,\mathbf{x}_n`` cresc pe termen lung geometric, cu un factor ``\lambda``, și au o **distribuție stabilă pe stadii**.

---
### ✍️ Exerciții
**Exercițiul 1 (metoda Euler în plan).** Pentru sistemul rîși–iepuri, porniți din ``(R_0, I_0) = (2,\ 0.5)``.

a) Calculați de mînă doi pași Euler cu ``\Delta t = 0.1``.

b) Calculați ``H = R - \ln R + I - \ln I`` în punctul inițial și după cei doi pași. Ce observați? Comparați cu subtitlul graficului din secțiunea 1.
"""

# ╔═╡ 9d4e2e82-bece-4ab8-b70c-eecfbaf065cc
hint(md"a) ``(R_1, I_1) = (1.9,\ 0.45)``; în acest punct vectorul este ``(1.9 \cdot 0.45 - 1.9,\ 0.45 - 1.9 \cdot 0.45) = (-1.045,\ -0.405)``, deci ``(R_2, I_2) = (1.7955,\ 0.4095)``. b) ``H(2, 0.5) = 2.5 - \ln 2 - \ln 0.5 = 2.5``, iar ``H(1.7955,\ 0.4095) \approx 2.5125``. Valoarea lui ``H`` crește ușor la fiecare pas: metoda Euler trece pe curbe de nivel tot mai „exterioare”, de aici spirala.")

# ╔═╡ 524c70f5-e5c8-49d7-a0a3-311e1e538c45
md"""
**Exercițiul 2 (null-cline).** Considerăm modelul de competiție
```math
\begin{aligned}
x' &= x\,(1 - x - 0.5\,y) \\
y' &= y\,(1 - y - 0.5\,x) .
\end{aligned}
```
a) Scrieți ecuațiile null-clinelor și găsiți toate cele patru puncte de echilibru.

b) Folosind criteriul invaziei, decideți dacă speciile coexistă.

c) Completați codul de mai jos (folosind funcțiile `cimp!`, `traiectorie` și `echilibru!` definite la începutul cursului) pentru a desena echilibrele și cîteva traiectorii.
"""

# ╔═╡ a1af1821-15ad-46c6-b3a2-9a6a56bf23ea
hint(md"a) ``x' = 0 \iff x = 0`` sau ``y = 2(1 - x)``; ``y' = 0 \iff y = 0`` sau ``y = 1 - 0.5x``. Echilibre: ``(0, 0)``, ``(1, 0)``, ``(0, 1)`` și, rezolvînd ``x + 0.5y = 1``, ``0.5x + y = 1``, punctul ``(2/3,\ 2/3)``. b) În ``(1, 0)``, specia ``y`` rară crește cu rata ``1 - 0.5 = 0.5 > 0``, deci poate invada; la fel ``x`` în ``(0, 1)``. Speciile coexistă, iar ``(2/3,\ 2/3)`` este stabil.")

# ╔═╡ 9b2a1875-5cfb-468b-88ef-d90924636bf6
# EXERCIȚIUL 2 c): completați și rulați
let
	F(u) = (u[1] * (1 - u[1] - 0.5u[2]), u[2] * (1 - u[2] - 0.5u[1]))
	fig = Figure(size = (600, 560))
	ax = Axis(fig[1, 1], title = "Exercițiul 2", xlabel = "x", ylabel = "y", aspect = 1,
		limits = (-0.02, 1.5, -0.02, 1.5))
	cimp!(ax, F, (0, 1.5), (0, 1.5))
	xs = range(0, 1.5, length = 100)
	lines!(ax, xs, 2 .* (1 .- xs), color = :darkorange, linewidth = 3, label = "x' = 0")
	lines!(ax, xs, 1 .- 0.5 .* xs, color = :teal, linewidth = 3, label = "y' = 0")
	# adăugați echilibrele, de exemplu: echilibru!(ax, 1, 0, false)
	# și cîteva traiectorii, de exemplu:
	X, Y, _ = traiectorie(F, (0.1, 1.4), 30)
	lines!(ax, X, Y, color = :crimson, linewidth = 2.5)
	axislegend(ax, position = :rt, backgroundcolor = (:white, 0.85))
	fig
end

# ╔═╡ 8c021753-4dc0-4545-9603-b96502a9dbee
md"""
**Exercițiul 3 (un model pe stadii pentru broaște).** O populație de broaște are două stadii, **mormoloci** ``T`` și **adulți** ``A``, iar modelul anual este
```math
\begin{pmatrix} T_{n+1} \\ A_{n+1} \end{pmatrix} =
\begin{pmatrix} 0 & 2 \\ 0.5 & 0.2 \end{pmatrix}
\begin{pmatrix} T_n \\ A_n \end{pmatrix} .
```
a) Ce semnificație biologică au numerele ``0``, ``2``, ``0.5`` și ``0.2`` din matrice?

b) Pornind cu ``100`` de mormoloci și niciun adult, calculați de mînă populația în primii trei ani. Ce observați despre structura populației?

c) Rulați codul de mai jos pentru a estima factorul de creștere pe termen lung ``\lambda``. Populația crește sau scade?

d) Folosind relația ``(1 - a)(1 - s) = fm`` din secțiunea 5, găsiți fracțiunea minimă de mormoloci care trebuie să devină adulți pentru ca populația să nu scadă.
"""

# ╔═╡ 87b40616-f8da-4b3d-b4ea-3adc6389b557
hint(md"a) ``0``: niciun mormoloc nu rămîne mormoloc (într-un an fie devine broască, fie moare); ``2``: fecunditatea, numărul mediu de mormoloci produși de un adult; ``0.5``: jumătate dintre mormoloci devin adulți; ``0.2``: ``20\%`` dintre adulți supraviețuiesc. b) ``(0,\ 50)``, ``(100,\ 10)``, ``(20,\ 52)``: structura oscilează puternic la început. c) ``\lambda \approx 1.105``: populația crește cu circa ``10.5\%`` pe an. d) ``1 \cdot 0.8 = 2m``, deci ``m = 0.4``: cel puțin ``40\%`` dintre mormoloci trebuie să ajungă adulți.")

# ╔═╡ 413db01b-faf3-4fc1-8896-672698fc6d4c
# EXERCIȚIUL 3 c): completați și rulați
let
	L = [0 2; 0.5 0.2]
	x = [100.0, 0.0]
	N = [sum(x)]
	for n in 1:40
		x = L * x
		push!(N, sum(x))
	end
	N[end] / N[end-1]   # estimarea factorului de creștere λ
end

# ╔═╡ 22cc6755-37fe-4602-a552-5b079d10bf7a
md"""
**Exercițiul 4 (discuție: reintroducerea unei specii).** 🌲 Într-o pădure modelul căprioare–elani are ``\alpha = 2`` și ``\beta = 1`` (secțiunea 4.2), iar acum trăiesc acolo doar căprioare, la capacitatea de suport ``D = 3``. Administrația vrea să reintroducă elanii.

a) Ce se întîmplă dacă sînt aduși cîțiva elani (de exemplu, ``M = 0.1``)? Explicați folosind criteriul invaziei.

b) Folosind figura bazinelor de atracție, estimați cîți elani ar trebui aduși deodată pentru ca reintroducerea să reușească.

c) Propuneți o strategie mai ieftină, care combină o intervenție temporară asupra căprioarelor cu aducerea elanilor. De ce funcționează?
"""

# ╔═╡ 4b7529b4-5f31-48f1-afad-357931053c94
hint(md"a) Elanii dispar: lîngă ``(3, 0)`` avem ``M' \approx M(2 - 3\beta) = -M < 0``, iar ``(3, 0)`` este stabil. b) Punctul ``(3, M_0)`` trebuie să fie deasupra separatoarei, care la ``D = 3`` trece pe la ``M \approx 2.2``: ar trebui aduși mai mulți elani decît capacitatea de suport a elanilor (``2``)! c) Dacă numărul căprioarelor este redus temporar (de exemplu, prin vînătoare) la ``D = 1.5``, separatoarea coboară la ``M \approx 1.33``: sînt necesari mai puțini elani. Odată trecut pragul, sistemul ajunge singur în noul echilibru și intervenția poate înceta.")

# ╔═╡ Cell order:
# ╟─a1b2c3d4-e5f6-7890-abcd-ef1234567890
# ╟─daaa48f8-1afa-4828-8c41-aa2134dcd652
# ╟─c10c3d78-56e0-4409-95d4-279cd27a0a08
# ╟─b2c3d4e5-f6a7-8901-bcde-f12345678901
# ╟─f0b5a93d-48f1-4539-99cd-7b78d49b8a6c
# ╟─72d70b30-e0f6-4aa4-82aa-db9f5a09b43d
# ╟─c3d4e5f6-a7b8-9012-cdef-123456789012
# ╟─c661a3e2-1f7a-4875-88ea-7fb972bb428b
# ╟─56324280-d312-4cf4-8e3c-de172f219a8c
# ╟─c7625076-300c-4814-bf6a-eeccb0ab6047
# ╟─07fa44dc-1187-4e5c-88bd-8fba4eb57442
# ╟─1fead660-881c-4b32-90db-fd06d113ac2c
# ╟─e23780cd-5bbe-4343-b87b-ebd49295b640
# ╟─9b6c3fd8-d1a1-483d-bf2d-45a7c2f72132
# ╟─7cb25802-be2e-40da-b4ef-08422d5b4976
# ╟─d95276d8-527f-44e4-8035-670f3c041851
# ╟─2d098eae-6e9a-4ac7-a892-851afdea616a
# ╟─a3b4c5d6-e7f8-9012-abcd-123456789012
# ╟─4922edb7-72ae-43f4-866a-e9b15faae600
# ╟─b4c5d6e7-f8a9-0123-bcde-234567890123
# ╟─75136446-2e33-4c5c-b968-cd2c3a2cc47e
# ╟─74f4f50b-a680-439d-9cbe-d0af495811fe
# ╟─d6e7f8a9-b0c1-2345-defa-456789012345
# ╟─e7f8a9b0-c1d2-3456-efab-567890123456
# ╟─f8a9b0c1-d2e3-4567-fabc-678901234567
# ╟─788dd5e1-9da2-4081-89b6-5c139d276b82
# ╟─f3e849c7-86fa-4012-bb33-98c4c45cad3b
# ╟─823027a6-23f4-40ab-b647-d969fd6362ab
# ╟─a9b0c1d2-e3f4-5678-abcd-789012345678
# ╟─c1d2e3f4-a5b6-7890-cdef-901234567890
# ╟─e5942d67-a5cd-4b6b-9eb4-5d206ea9f431
# ╟─abd99d36-5deb-42d4-a713-9cb7d605b72d
# ╟─e3f4a5b6-c7d8-9012-efab-123456789012
# ╟─9d4e2e82-bece-4ab8-b70c-eecfbaf065cc
# ╟─524c70f5-e5c8-49d7-a0a3-311e1e538c45
# ╟─a1af1821-15ad-46c6-b3a2-9a6a56bf23ea
# ╠═9b2a1875-5cfb-468b-88ef-d90924636bf6
# ╟─8c021753-4dc0-4545-9603-b96502a9dbee
# ╟─87b40616-f8da-4b3d-b4ea-3adc6389b557
# ╠═413db01b-faf3-4fc1-8896-672698fc6d4c
# ╟─22cc6755-37fe-4602-a552-5b079d10bf7a
# ╟─4b7529b4-5f31-48f1-afad-357931053c94
