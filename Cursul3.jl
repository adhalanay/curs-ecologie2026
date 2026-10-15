### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> title = "Cursul III"
#> description = "Integrala și ecuații diferențiale"
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

# ╔═╡ 5597b08d-2218-4f9d-b46b-228c1161f894
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using PlutoTeachingTools
	using OrdinaryDiffEq
	using QuadGK
	using Printf
end

# ╔═╡ a17d6370-a346-4630-b5a3-0f07fbf5474a
TableOfContents()

# ╔═╡ 5fe3884c-8911-4bf6-b88d-a5cb6e8ee05c
HTML("""
<style>
    :root {
        --bg-color: #282a36;
        --text-color: #f8f8f2;
        --accent-color: #8be9fd;
        --code-bg: #44475a;
        --code-text: #ff79c6;
        --quote-bg: #343746;
        --quote-border: #50fa7b;
        --quote-text: #f1fa8c;
    }
    body {
        background-color: var(--bg-color) !important;
        color: var(--text-color) !important;
    }
    main {
        max-width: 1900px !important;
    }
    .pluto-output {
        background-color: transparent !important;
        color: var(--text-color) !important;
    }
    .pluto-output h1, .pluto-output h2, .pluto-output h3 {
        color: var(--accent-color) !important;
        border-bottom: 2px solid #44475a !important;
    }
    .pluto-output h1 { border-bottom: 3px solid var(--accent-color) !important; }
    .pluto-output blockquote {
        border-left: 5px solid var(--quote-border) !important;
        background: var(--quote-bg) !important;
        color: var(--quote-text) !important;
    }
    .pluto-output code {
        background: var(--code-bg) !important;
        color: var(--code-text) !important;
    }
    .pluto-output a { color: #bd93f9 !important; }
    .pluto-output img {
        border-radius: 12px;
        box-shadow: 0 8px 16px rgba(0,0,0,0.5);
        margin: 20px 0;
        transition: transform 0.2s;
        max-width: 100%;
    }
    .pluto-output img:hover { transform: scale(1.02); }
    .pluto-output ul li::marker { color: var(--accent-color); }
    .pluto-cell { background-color: transparent !important; }
</style>
""")

# ╔═╡ 60547c1c-e309-4606-adca-f00e7b9c95f3
begin
	culoare_fundal = "#282a36"
	stil_ax = (
		backgroundcolor = culoare_fundal,
		xgridcolor = :grey30, ygridcolor = :grey30,
		xtickcolor = :white, ytickcolor = :white,
		xticklabelcolor = :white, yticklabelcolor = :white,
		xlabelcolor = :white, ylabelcolor = :white,
		titlecolor = :white, subtitlecolor = :grey80,
		leftspinecolor = :grey60, rightspinecolor = :grey60,
		topspinecolor = :grey60, bottomspinecolor = :grey60,
	)
	stil_leg = (backgroundcolor = (:black, 0.3), framecolor = :grey60, labelcolor = :white)
	figura(; kw...) = Figure(; backgroundcolor = culoare_fundal, kw...)
	# paleta Dracula
	cyan, roz, verde, galben = "#8be9fd", "#ff79c6", "#50fa7b", "#f1fa8c"
	portocaliu, mov, rosu = "#ffb86c", "#bd93f9", "#ff5555"
	nothing
end

# ╔═╡ cd5e96c4-a46b-4ba0-8254-abe5a518a50d
begin
	# Graficul lui f și linia de fază pe axa Ox: săgeată verde unde x crește (f > 0),
	# roșie unde x scade (f < 0); echilibrele stabile sînt pline, cele instabile goale.
	function linie_de_faza!(ax, f, echilibre, xmin, xmax; culoare = cyan)
		xs = range(xmin, xmax, length = 400)
		hlines!(ax, 0, color = :grey60, linewidth = 1)
		lines!(ax, xs, f.(xs), color = culoare, linewidth = 3)
		repere = sort(unique([xmin; echilibre; xmax]))
		for (p, q) in zip(repere[1:end-1], repere[2:end])
			s = sign(f((p + q) / 2))
			s == 0 && continue
			L = 0.4 * (q - p)
			arrows2d!(ax, [(p + q) / 2 - s * L / 2], [0.0], [s * L], [0.0],
				color = s > 0 ? verde : rosu, shaftwidth = 4, tipwidth = 16, tiplength = 14)
		end
		δ = 1e-3 * (xmax - xmin)
		for e in echilibre
			stabil = f(e - δ) > 0 && f(e + δ) < 0
			scatter!(ax, [e], [0.0], markersize = 18, strokewidth = 2.5,
				strokecolor = :white, color = stabil ? :white : culoare_fundal)
		end
	end

	# soluția explicită a ecuației logistice x' = r x (1 - x/K)
	logistica(t, x0; r, K) = K / (1 + (K / x0 - 1) * exp(-r * t))

	# soluția numerică a ecuației x' = f(x) pe [0, T]
	function traiectorie(f, x0, T)
		sol = solve(ODEProblem((u, p, t) -> f(u), float(x0), (0.0, T)), Tsit5(),
			saveat = T / 400, reltol = 1e-8, abstol = 1e-10)
		sol.t, sol.u
	end
	nothing
end

# ╔═╡ ec53dd46-1054-4dda-b83a-2ca3070e6edd
md"""
# 🧮 Cursul 3 — Integrala și ecuații diferențiale

**Durata estimată:** 100 minute
"""

# ╔═╡ 72400236-1ff8-426f-92e6-0647df943bb4
md"""
## 🎯 Obiectivele lecției
1. Să înțelegem **integrala** ca operația inversă derivării și ca **arie sub grafic**.
2. Să aproximăm integrale prin **sume Riemann** și să recunoaștem în ele **metoda Euler** din cursul trecut.
3. Să rezolvăm ecuațiile de **creștere și descompunere exponențială** și să calculăm timpul de dublare și de înjumătățire.
4. Să găsim **punctele de echilibru** ale unei ecuații ``x' = f(x)`` și să le determinăm **stabilitatea** (linia de fază și criteriul derivatei).
5. Să aplicăm aceste instrumente la **modelul logistic**, la **pescuit** și la **efectul Allee**.

"""

# ╔═╡ 6b69614e-1a69-4aa3-a0f7-79bfb85b32c7
md"""
## 🔁 0. Ce știm deja 
În cursul trecut am învățat să trecem de la o funcție ``x(t)`` la **rata ei de schimbare**:
```math
x'(t) = \lim_{\Delta t \to 0} \frac{x(t+\Delta t) - x(t)}{\Delta t},
```
și am folosit **metoda Euler**, ``x_{n+1} = x_n + \Delta t \cdot f(x_n)``, pentru a desena traiectoriile unei ecuații ``x' = f(x)``.

Astăzi facem drumul invers: **cunoaștem rata de schimbare și vrem să aflăm starea**. Vom întîlni două situații:
- rata depinde doar de timp, ``x' = g(t)`` (vitezometrul unei mașini, debitul unui rîu) ⟶ **integrala**;
- rata depinde de starea însăși, ``x' = f(x)`` (o populație crește cu atît mai repede cu cît este mai mare) ⟶ **ecuație diferențială**. De cele mai multe ori nu o putem rezolva explicit, dar putem afla **ce se întîmplă pe termen lung** fără să o rezolvăm.
"""

# ╔═╡ fb034a9e-44d1-4b25-b25a-245a62423cf4
md"""
## ↩️ 1. De la derivată la integrală 
Derivata unei funcții se calculează destul de ușor, pornind de la cîteva formule și aplicînd reguli de calcul. Putem inversa procesul, adică să determinăm ``F`` dacă știm doar ``F' = f``? O astfel de funcție ``F`` se numește **primitivă** a lui ``f``, iar procedura se numește **integrare**.

Pentru funcții simple nu este greu: dacă ``F'(x) = 2x``, atunci ``F(x) = x^2 + c`` pentru orice constantă ``c``, deoarece derivata unei constante este zero. Primitiva este deci determinată **pînă la o constantă**, pe care o fixăm cunoscînd o singură valoare, de exemplu ``F(0)``.

| ``f(x)`` | ``x^n\ (n \neq -1)`` | ``1/x`` | ``e^{x}`` | ``e^{kx}`` | ``\cos x`` | ``\sin x`` |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| ``F(x)`` | ``\dfrac{x^{n+1}}{n+1}`` | ``\ln\lvert x\rvert`` | ``e^{x}`` | ``\dfrac{e^{kx}}{k}`` | ``\sin x`` | ``-\cos x`` |

> 💡 Fiecare coloană se verifică **derivînd** rîndul de jos: trebuie să obținem rîndul de sus.

În general însă primitiva este mult mai greu de găsit decît derivata, iar uneori nici nu se poate exprima prin funcții elementare (de exemplu pentru ``e^{-x^2}``, curba lui Gauss). Mai mult, de multe ori **nu avem o formulă**, ci doar măsurători: avem înregistrările vitezometrului unei mașini (viteza ``X'(t)``) și vrem să aflăm distanța parcursă ``X(t)``.
"""

# ╔═╡ 2ebd78a7-8a94-4077-812d-86e308566903
begin
	viteza_masurata(t) = t - cos(t) * sin(t)
	t_masurat = range(0, 10, length = 15)
	fig1 = figura(size = (1000, 400))
	ax1 = Axis(fig1[1, 1]; stil_ax..., xlabel = "t  [s]", ylabel = "X'(t)  [m/s]",
		title = "Vitezometrul: cunoaștem X'(t) doar în cîteva momente")
	lines!(ax1, t_masurat, viteza_masurata.(t_masurat), color = cyan, linewidth = 2)
	scatter!(ax1, t_masurat, viteza_masurata.(t_masurat), color = roz, markersize = 12,
		strokewidth = 2, strokecolor = :white)
	fig1
end

# ╔═╡ 78cd9df2-8b48-41cd-8d9f-4ee9979bcde3
md"""
### ✅ Verificare rapidă
Găsiți funcția ``F`` cu ``F'(x) = 3x^2 + 2`` și ``F(0) = 1``.
"""

# ╔═╡ a26b0c6f-f66b-427c-b06f-eb182529ec47
hint(md"``F(x) = x^3 + 2x + c``, iar ``F(0) = c = 1``, deci ``F(x) = x^3 + 2x + 1``. Verificare: ``F'(x) = 3x^2 + 2``.")

# ╔═╡ e0da1c43-41c2-475e-a898-503424d2d511
md"""
## 📏 2. Sume Riemann și teorema fundamentală 

### 2.1. Metoda Euler pentru ``X' = g(t)``
Ne amintim de definiția derivatei:
```math
X'(t) \approx \frac{X(t+\Delta t) - X(t)}{\Delta t}
\quad\Longrightarrow\quad
X(t+\Delta t) \approx X(t) + \Delta t \cdot X'(t).
```
Este exact **metoda Euler**. Repetăm construcția pas cu pas pentru ``t = 0, \Delta t, 2\Delta t, \dots``:
```math
\begin{aligned}
X(\Delta t)  &\approx X(0) + X'(0)\,\Delta t, \\
X(2\Delta t) &\approx X(\Delta t) + X'(\Delta t)\,\Delta t = X(0) + X'(0)\,\Delta t + X'(\Delta t)\,\Delta t, \\
X(3\Delta t) &\approx X(0) + X'(0)\,\Delta t + X'(\Delta t)\,\Delta t + X'(2\Delta t)\,\Delta t.
\end{aligned}
```
După ``n`` pași, la momentul ``t = n\Delta t``:
```math
X(t) \approx X(0) + \sum_{k=0}^{n-1} X'(k\Delta t)\,\Delta t .
```
Suma din dreapta se numește **sumă Riemann**. Spre deosebire de cazul ``x' = f(x)``, aici rata ``X'`` nu depinde de ``X``, deci fiecare termen se calculează direct, fără să cunoaștem pașii anteriori.

### 2.2. Interpretarea geometrică: aria de sub grafic
Fiecare termen ``X'(k\Delta t)\,\Delta t`` este aria unui dreptunghi cu baza ``\Delta t`` și înălțimea ``X'(k\Delta t)``: presupunem că viteza este **constantă** pe fiecare subinterval. Suma Riemann este deci aria unei „scări” de dreptunghiuri care aproximează **aria de sub grafic**.

**Exemplu.** Viteza unei mașini care frînează este ``V(t) = 2\sqrt{1000 - t^3}`` m/s. Mașina pornește cu ``V(0) \approx 63`` m/s (``\approx 228`` km/h) și se oprește după ``10`` s. Ce distanță parcurge? Primitiva lui ``V`` nu se exprimă prin funcții elementare, așa că folosim sume Riemann.
"""

# ╔═╡ 10a166b2-1c43-4a69-879b-190cf49323d6
begin
	V(t) = 2 * sqrt(1000 - t^3)

	# suma Riemann la stînga (metoda Euler) pentru ∫ₐᵇ g(t) dt cu n dreptunghiuri
	function suma_riemann(g, a, b, n)
		Δ = (b - a) / n
		sum(g(a + k * Δ) for k in 0:n-1) * Δ
	end

	dist_exacta = quadgk(V, 0, 10)[1]
	t_fin = range(0, 10, length = 300)
	X_exact = [quadgk(V, 0, t)[1] for t in t_fin]
	nothing
end

# ╔═╡ 4ec95b87-81c2-40c6-b3cf-2ba63c805ec1
md"""
Numărul de dreptunghiuri ``n`` = $(@bind N_riemann PlutoUI.Slider([1, 2, 4, 5, 10, 20, 50, 100], default = 5, show_value = true))
"""

# ╔═╡ 6af67b03-6a21-4113-9578-a69cc795db98
begin
	Δn = 10 / N_riemann
	tk = Δn .* (0:N_riemann-1)            # capetele din stînga ale subintervalelor
	S_n = sum(V.(tk)) * Δn
	fig2 = figura(size = (1100, 420))
	ax2a = Axis(fig2[1, 1]; stil_ax..., xlabel = "t  [s]", ylabel = "V(t)  [m/s]",
		title = "Aria de sub grafic ≈ suma ariilor dreptunghiurilor",
		subtitle = @sprintf("n = %d:   Sₙ = %.2f m,   exact = %.2f m", N_riemann, S_n, dist_exacta))
	barplot!(ax2a, tk .+ Δn / 2, V.(tk), width = Δn, gap = 0,
		color = (roz, 0.35), strokecolor = roz, strokewidth = 1.5)
	lines!(ax2a, t_fin, V.(t_fin), color = cyan, linewidth = 3)
	ax2b = Axis(fig2[1, 2]; stil_ax..., xlabel = "t  [s]", ylabel = "X(t)  [m]",
		title = "Distanța parcursă: Euler vs. exact",
		subtitle = @sprintf("eroare la t = 10:  %.2f m", abs(S_n - dist_exacta)))
	lines!(ax2b, t_fin, X_exact, color = cyan, linewidth = 3, label = "exact (quadgk)")
	scatterlines!(ax2b, [tk; 10], [0; cumsum(V.(tk)) .* Δn], color = roz,
		markersize = 8, linewidth = 2, label = "Euler / sume Riemann")
	axislegend(ax2b; position = :rb, stil_leg...)
	fig2
end

# ╔═╡ e82ffe61-be7b-4888-adaa-fa92dccf6f1d
md"""
Cu cît ``n`` crește (``\Delta t`` scade), cu atît dreptunghiurile urmăresc mai bine graficul. Trecînd la limită ``\Delta t \to 0``, suma Riemann devine **integrala definită** ``\int_0^t X'(s)\,ds``.

**Teorema fundamentală a analizei (Leibniz–Newton).** Dacă ``F' = f`` pe ``[a,b]``, atunci
```math
\int_a^b f(t)\,dt = F(b) - F(a), \qquad \text{în particular} \qquad X(t) = X(0) + \int_0^t X'(s)\,ds .
```
Integrala are deci două fețe: **aria de sub grafic** (pe care o aproximăm numeric prin sume Riemann) și **diferența valorilor unei primitive** (pe care o calculăm exact cînd cunoaștem o primitivă). Acolo unde ``f < 0``, aria se socotește cu semnul minus.

*Exemplu:* ``\displaystyle\int_0^2 x^2\,dx = \frac{2^3}{3} - \frac{0^3}{3} = \frac{8}{3}``.

Pentru mașina noastră nu avem o primitivă, dar putem vedea cît de repede converg sumele Riemann către valoarea exactă (calculată cu `quadgk`, o metodă numerică mult mai precisă):
"""

# ╔═╡ aa48499a-994d-4baa-9ea7-68327f17b91e
let
	linii = ["| ``n`` | ``\\Delta t`` | Sumă Riemann [m] | Eroare [m] |", "|:---:|:---:|:---:|:---:|"]
	for n in (10, 100, 1000, 10_000)
		S = suma_riemann(V, 0, 10, n)
		push!(linii, @sprintf("| %d | %g | %.4f | %.4f |", n, 10 / n, S, abs(S - dist_exacta)))
	end
	push!(linii, @sprintf("| exact | — | %.4f | — |", dist_exacta))
	Markdown.parse(join(linii, "\n"))
end

# ╔═╡ 82e1c110-a244-4d5c-92a0-a0cbc26721e6
md"""
Eroarea scade de aproximativ 10 ori cînd ``n`` crește de 10 ori: eroarea metodei Euler este **proporțională cu ``\Delta t``**.

### 2.3. Valoarea medie
Dacă împărțim distanța la durata mișcării obținem viteza medie. În general, **valoarea medie** a lui ``f`` pe ``[a,b]`` este
```math
\bar{f} = \frac{1}{b-a}\int_a^b f(t)\,dt ,
```
adică înălțimea dreptunghiului cu baza ``[a,b]`` care are aceeași arie ca regiunea de sub grafic. Pentru mașina noastră ``\bar V = \frac{1}{10}\int_0^{10} V(t)\,dt \approx`` $(round(dist_exacta / 10, digits = 1)) m/s.

> 🧬 În farmacologie, aria de sub curba concentrației unui medicament în sînge (**AUC**, *area under the curve*) măsoară expunerea totală a organismului la medicament. Ea se calculează exact ca mai sus, din măsurători făcute la intervale de timp.

### ✅ Verificare rapidă
Calculați ``\int_0^1 (3x^2 + 2)\,dx`` folosind teorema fundamentală. Care este valoarea medie a lui ``3x^2 + 2`` pe ``[0, 1]``?
"""

# ╔═╡ 9d53abcc-1c11-42a5-9f7c-51d660d232b5
hint(md"O primitivă este ``F(x) = x^3 + 2x``, deci integrala este ``F(1) - F(0) = 3``. Intervalul are lungimea ``1``, deci și valoarea medie este ``3``.")

# ╔═╡ a721c2c7-b8c2-4fd9-842f-bbc69d6fe2d6
md"""
## 📈 3. Creștere și descompunere exponențială 
Trecem acum la ecuații în care rata de schimbare depinde de **stare**. Cel mai simplu caz este cel din cursul trecut: dacă ``r = b + i - d - e`` este rata netă de creștere per capita (nașteri + imigrare − decese − emigrare), atunci
```math
x' = r\,x .
```
Nu putem integra direct ca la ``X' = g(t)``, pentru că membrul drept conține necunoscuta ``x(t)``. Trucul: împărțim la ``x`` și observăm că ``x'/x = (\ln x)'``. Deci ``(\ln x)' = r`` și, integrînd de la ``0`` la ``t`` cu teorema fundamentală,
```math
\ln x(t) - \ln x(0) = r\,t \quad\Longrightarrow\quad x(t) = x(0)\,e^{rt}.
```
Verificare: ``\big(x(0)e^{rt}\big)' = r\,x(0)e^{rt} = r\,x(t)``. ✔

- ``r > 0``: **creștere exponențială** (bacterii într-un mediu bogat, o specie invazivă la început);
- ``r < 0``: **descompunere exponențială** (dezintegrare radioactivă, eliminarea unui medicament din sînge).

**Timpul de dublare / de înjumătățire.** Căutăm ``T`` cu ``x(T) = 2x(0)``, adică ``e^{rT} = 2``, respectiv ``x(T)=\frac{1}{2}x(0)`` și ``e^{rT}=\frac{1}{2}``:
```math
T_2 = \frac{\ln 2}{r} \quad (r>0), \qquad\qquad T_{1/2} = \frac{\ln 2}{\lvert r\rvert} \quad (r<0).
```
Timpul de dublare **nu depinde de mărimea populației**: o cultură de *E. coli* (``T_2 \approx 20`` min) are nevoie de 20 de minute ca să treacă de la o mie la două mii de celule, dar și de la un milion la două milioane. Carbonul-14 are ``T_{1/2} \approx 5730`` ani, ceea ce permite datarea fosilelor.

Rata ``r`` = $(@bind r_exp PlutoUI.Slider(-1.0:0.05:1.0, default = 0.3, show_value = true))
"""

# ╔═╡ 762f11dd-a697-4a63-9e26-559b4e0148f7
begin
	t_exp = range(0, 10, length = 300)
	x_exp = exp.(r_exp .* t_exp)
	fig3 = figura(size = (1100, 420))
	ax3a = Axis(fig3[1, 1]; stil_ax..., xlabel = "t", ylabel = "x(t)",
		title = @sprintf("x(t) = x₀·e^(rt),  x₀ = 1,  r = %.2f", r_exp))
	ax3b = Axis(fig3[1, 2]; stil_ax..., xlabel = "t", ylabel = "x(t)", yscale = log10,
		ytickformat = vs -> [@sprintf("%.3g", v) for v in vs],
		title = "Scară semi-logaritmică: exponențiala devine dreaptă")
	culoare_exp = r_exp > 0 ? verde : (r_exp < 0 ? roz : galben)
	for ax in (ax3a, ax3b)
		lines!(ax, t_exp, x_exp, color = culoare_exp, linewidth = 3)
	end
	if r_exp != 0
		T = log(2) / abs(r_exp)
		nivel = r_exp > 0 ? 2.0 : 0.5
		eticheta = r_exp > 0 ? @sprintf("T₂ = ln2/r ≈ %.2f", T) : @sprintf("T½ = ln2/|r| ≈ %.2f", T)
		if T <= 10
			lines!(ax3a, [0, T, T], [nivel, nivel, 0], color = galben, linestyle = :dash)
			scatter!(ax3a, [T], [nivel], color = galben, markersize = 12)
			text!(ax3a, T, nivel, text = eticheta, color = galben,
				align = (:left, r_exp > 0 ? :top : :bottom), offset = (8, r_exp > 0 ? -6 : 6))
		end
	end
	ylims!(ax3a, 0, nothing)
	fig3
end

# ╔═╡ a71c9839-5b58-4c85-939e-de9c1053135e
md"""
### ✅ Verificare rapidă
O populație de iepuri se dublează la fiecare 3 ani. Care este rata ``r``? De cîte ori crește populația în 12 ani?
"""

# ╔═╡ 05920e3c-c24d-45b0-a65a-c5c6df52dd4c
hint(md"``r = \ln 2 / 3 \approx 0.23`` pe an. În 12 ani au loc 4 dublări, deci populația crește de ``2^4 = 16`` ori (sau: ``e^{12r} = e^{4\ln 2} = 16``).")

# ╔═╡ 9b1be66f-ea19-41c6-8235-e90fa1257fd8
md"""
## 🌱 4. Modelul logistic și punctele de echilibru 
Modelul exponențial este valabil doar pe termen scurt: el prezice o creștere nelimitată, ceea ce este nerealist biologic. Resursele fiind finite, rata de creștere per capita trebuie să scadă cînd populația crește. Cea mai simplă alegere duce la **ecuația logistică** (Cursul 2):
```math
x' = f(x) = r\,x\left(1 - \frac{x}{K}\right),
```
unde ``K`` este **capacitatea de suport** a mediului. Această ecuație se poate rezolva explicit, separînd variabilele ca mai sus:
```math
x(t) = \frac{K}{1 + \left(\dfrac{K}{x_0} - 1\right)e^{-rt}} .
```
Formula nu este însă cel mai bun mod de a înțelege ce se întîmplă. Să privim **graficul lui ``f``** alături de soluții (``r = 0.5``, ``K = 100``):
"""

# ╔═╡ 30c59120-bee5-4743-8606-38ba3d9c1d43
begin
	r_log, K_log = 0.5, 100.0
	f_log(x) = r_log * x * (1 - x / K_log)
	fig4 = figura(size = (1100, 440))
	ax4a = Axis(fig4[1, 1]; stil_ax..., xlabel = "x", ylabel = "x' = f(x)",
		title = "f(x) = r x (1 − x/K) și linia de fază")
	vlines!(ax4a, K_log / 2, color = galben, linestyle = :dash)
	text!(ax4a, K_log / 2, f_log(K_log / 2), text = "maxim în K/2", color = galben,
		align = (:center, :bottom), offset = (0, 6))
	linie_de_faza!(ax4a, f_log, [0.0, K_log], 0.0, 1.35K_log)
	ylims!(ax4a, nothing, 1.4f_log(K_log / 2))
	ax4b = Axis(fig4[1, 2]; stil_ax..., xlabel = "t", ylabel = "x(t)",
		title = "Soluții pentru diferite valori inițiale x₀")
	hlines!(ax4b, K_log / 2, color = galben, linestyle = :dash)
	hlines!(ax4b, K_log, color = :grey60, linestyle = :dash)
	t_log = range(0, 25, length = 300)
	for (x0, c) in zip([2, 10, 30, 70, 135], [rosu, portocaliu, galben, verde, mov])
		lines!(ax4b, t_log, logistica.(t_log, x0; r = r_log, K = K_log),
			color = c, linewidth = 3, label = "x₀ = $x0")
	end
	axislegend(ax4b; position = :rb, stil_leg...)
	fig4
end

# ╔═╡ 8a77a20a-92cb-46a3-aba5-00affcf583a4
md"""
**Puncte de echilibru.** Un punct ``x^*`` cu ``f(x^*) = 0`` se numește **punct de echilibru** (stare staționară): dacă pornim din ``x^*``, atunci ``x' = 0`` și populația rămîne constantă pentru totdeauna. Pentru ecuația logistică ``f(x) = 0 \iff x^* = 0`` sau ``x^* = K``.

Din graficul lui ``f`` citim tot comportamentul:
- unde ``f(x) > 0`` (graficul deasupra axei), ``x`` **crește** — săgeată verde spre dreapta;
- unde ``f(x) < 0`` (graficul sub axă), ``x`` **scade** — săgeată roșie spre stînga;
- populația crește cel mai repede acolo unde ``f`` este maximă, adică la ``x = K/2``: acesta este **punctul de inflexiune** al curbei în formă de S.

Toate soluțiile cu ``x_0 > 0`` tind către ``K``, fie de jos, fie de sus, dar **niciuna nu traversează** nivelul ``K``: două soluții ale aceleiași ecuații ``x' = f(x)`` nu se pot intersecta.
"""

# ╔═╡ f231db04-5759-4209-b061-4d34f7f060a7
md"""
## ⚖️ 5. Stabilitatea echilibrelor
Un echilibru este **stabil** (atractor) dacă soluțiile care pornesc aproape de el se apropie de el, și **instabil** (repulsor) dacă se îndepărtează. În natură observăm doar echilibrele stabile: o mică perturbație (o iarnă grea, o epidemie) scoate sistemul dintr-un echilibru instabil și acesta nu se mai întoarce. Putem afla stabilitatea **fără să rezolvăm ecuația**.

### 5.1. Metoda semnului (linia de fază)
Fiind continuă, ``f`` nu își poate schimba semnul decît trecînd prin zero, deci între două echilibre consecutive ``f`` are semn constant. Este suficient să calculăm ``f`` într-un **punct test** din fiecare interval:
- ``f > 0`` la stînga lui ``x^*`` și ``f < 0`` la dreapta ⟹ săgețile arată spre ``x^*`` ⟹ **stabil** (● plin);
- ``f < 0`` la stînga și ``f > 0`` la dreapta ⟹ săgețile pleacă din ``x^*`` ⟹ **instabil** (○ gol);
- același semn de ambele părți ⟹ **semistabil** (atrage dintr-o parte, respinge din cealaltă).

Desenul obținut pe axa ``Ox`` se numește **linia de fază**. Pentru ecuația logistică: ``f(K/2) > 0`` și ``f(2K) < 0``, deci ``K`` este stabil, iar ``0`` este instabil.

### 5.2. Criteriul derivatei (liniarizare)
Fie ``x^*`` un echilibru și ``u = x - x^*`` o mică abatere. Aproximăm ``f`` prin tangenta sa în ``x^*`` (Cursul 2):
```math
u' = x' = f(x^* + u) \approx \underbrace{f(x^*)}_{=\,0} + f'(x^*)\,u .
```
Abaterea satisface deci ecuația exponențială ``u' = f'(x^*)\,u``, cu soluția ``u(t) = u(0)\,e^{f'(x^*)\,t}`` (secțiunea 3)! Prin urmare:
- ``f'(x^*) < 0 \implies`` abaterea se stinge exponențial ⟹ ``x^*`` **stabil**;
- ``f'(x^*) > 0 \implies`` abaterea crește exponențial ⟹ ``x^*`` **instabil**;
- ``f'(x^*) = 0 \implies`` criteriul **nu decide** (``x' = -x^3`` și ``x' = x^3`` au ambele ``f'(0) = 0``, dar ``0`` este stabil pentru prima ecuație și instabil pentru a doua).

În plus, ``\lvert f'(x^*)\rvert`` ne spune **cît de repede** revine sistemul la echilibru: timpul caracteristic de revenire este ``1/\lvert f'(x^*)\rvert``.

*Exemplu (logistic):* ``f(x) = rx - rx^2/K``, deci ``f'(x) = r\left(1 - \dfrac{2x}{K}\right)``:
- ``f'(0) = r > 0`` ⟹ ``0`` este instabil;
- ``f'(K) = -r < 0`` ⟹ ``K`` este stabil, iar o populație perturbată revine la ``K`` ca ``e^{-rt}``.

### 5.3. Exemplu: pescuit cu efort constant
O populație de pești crește logistic cu ``r = 0.2`` pe an și ``K = 1000`` tone, iar flota pescuiește în fiecare an o fracțiune ``E = 0.1`` din stoc (efort constant):
```math
x' = 0.2\,x\left(1 - \frac{x}{1000}\right) - 0.1\,x .
```
Echilibrele: ``x\big(0.2(1 - x/1000) - 0.1\big) = 0 \implies x^* = 0`` sau ``x^* = 500``. Puncte test: ``f(100) = 8 > 0`` și ``f(1000) = -100 < 0``. Cu derivata: ``f'(x) = 0.1 - 0.0004\,x``, deci ``f'(0) = 0.1 > 0`` (instabil) și ``f'(500) = -0.1 < 0`` (stabil). Captura la echilibru este ``0.1 \cdot 500 = 50`` tone pe an.
"""

# ╔═╡ 37602046-94fa-4853-a4cb-ffe93f6841a2
begin
	f_pesc(x) = 0.2x * (1 - x / 1000) - 0.1x
	fig5 = figura(size = (1100, 440))
	ax5a = Axis(fig5[1, 1]; stil_ax..., xlabel = "x  [tone]", ylabel = "x'",
		title = "Linia de fază cu pescuit (gri: fără pescuit)")
	x_pesc = range(0, 1100, length = 300)
	lines!(ax5a, x_pesc, 0.2 .* x_pesc .* (1 .- x_pesc ./ 1000), color = :grey60,
		linestyle = :dash, linewidth = 2)
	linie_de_faza!(ax5a, f_pesc, [0.0, 500.0], 0.0, 1100.0)
	ax5b = Axis(fig5[1, 2]; stil_ax..., xlabel = "t  [ani]", ylabel = "x(t)  [tone]",
		title = "Stocul de pește tinde la 500 de tone")
	hlines!(ax5b, 500, color = :white, linestyle = :dash)
	hlines!(ax5b, 1000, color = :grey60, linestyle = :dash)
	t_pesc = range(0, 100, length = 300)
	for (x0, c) in zip([50, 300, 900], [rosu, galben, verde])
		lines!(ax5b, t_pesc, logistica.(t_pesc, x0; r = 0.1, K = 500.0),
			color = c, linewidth = 3, label = "x₀ = $x0")
	end
	axislegend(ax5b; position = :rt, stil_leg...)
	fig5
end

# ╔═╡ 2c81e4da-6b65-4831-a7bc-e7cae78a78db
md"""
> 💡 Observați că ``x' = 0.1\,x(1 - x/500)``: pescuitul proporțional produce tot un model logistic, cu rata de creștere și capacitatea de suport mai mici. Populația nu este pusă în pericol. Situația se schimbă radical dacă se pescuiește o cantitate **fixă** ``H`` pe an (Exercițiul 4).

### ✅ Verificare rapidă
Găsiți echilibrele ecuației ``x' = x^2 - 1`` și stabilitatea lor, cu ambele metode.
"""

# ╔═╡ 70a6b506-4529-4ce3-b79e-6067759a8ff7
hint(md"``x^* = \pm 1``. Derivata: ``f'(x) = 2x``, deci ``f'(-1) = -2 < 0`` (stabil) și ``f'(1) = 2 > 0`` (instabil). Puncte test: ``f(-2) = 3 > 0``, ``f(0) = -1 < 0``, ``f(2) = 3 > 0``.")

# ╔═╡ d7fac950-643d-424a-b783-f774ac30c532
md"""
## 🐋 6. Efectul Allee 
În modelul logistic rata de creștere per capita ``f(x)/x = r(1 - x/K)`` este maximă cînd populația este foarte mică. Pentru multe specii acest lucru este fals: sub un anumit prag indivizii **nu își mai găsesc parteneri** (balene, rinoceri), nu mai pot **vîna sau se apăra în grup** (cîinii sălbatici africani, suricatele), iar plantele nu mai sînt **polenizate** suficient. Acesta este **efectul Allee** (W. C. Allee, 1931), modelat prin
```math
x' = f(x) = r\,x\left(1 - \frac{x}{K}\right)\left(\frac{x}{a} - 1\right), \qquad 0 < a < K,
```
unde ``a`` este **pragul Allee**. Acum rata per capita ``r(1 - x/K)(x/a - 1)`` este **negativă** pentru ``x < a``.

**Echilibre:** ``x^* = 0``, ``x^* = a``, ``x^* = K``. **Stabilitate** (criteriul derivatei):
```math
f'(0) = -r < 0 \;\; \text{(stabil)}, \qquad
f'(a) = r\left(1 - \frac{a}{K}\right) > 0 \;\; \text{(instabil — prag critic)}, \qquad
f'(K) = r\left(1 - \frac{K}{a}\right) < 0 \;\; \text{(stabil)}.
```
Spre deosebire de modelul logistic, **extincția ``x = 0`` devine stabilă**. Sistemul este **bistabil**: soarta populației depinde de unde pornește (``r = 0.5``, ``a = 20``, ``K = 100``):
"""

# ╔═╡ e07d6e57-38cb-493d-9aba-7a957c7660e6
hint(md"Derivați cu regula produsului. Într-un echilibru, unul dintre factori este zero, deci rămîne doar termenul în care se derivează exact acel factor. De exemplu, în ``x = a``: ``f'(a) = r\,a\,(1 - a/K)\cdot\frac{1}{a} = r(1 - a/K)``.")

# ╔═╡ 8f034303-c135-4b30-8cf4-0f63c8f9526c
begin
	r_A, a_A, K_A = 0.5, 20.0, 100.0
	f_A(x) = r_A * x * (1 - x / K_A) * (x / a_A - 1)
	fig6 = figura(size = (1100, 440))
	ax6a = Axis(fig6[1, 1]; stil_ax..., xlabel = "x", ylabel = "x' = f(x)",
		title = "Efectul Allee: f(x) și linia de fază")
	linie_de_faza!(ax6a, f_A, [0.0, a_A, K_A], 0.0, 125.0, culoare = portocaliu)
	for (x, s) in ((a_A, "a"), (K_A, "K"))
		text!(ax6a, x, 0, text = s, color = :white, align = (:center, :top), offset = (0, -14))
	end
	ax6b = Axis(fig6[1, 2]; stil_ax..., xlabel = "t", ylabel = "x(t)",
		title = "Sub prag: extincție (roșu). Peste prag: K (verde)")
	hlines!(ax6b, [a_A, K_A], color = :grey60, linestyle = :dash)
	for x0 in [3, 8, 13, 17, 19, 21, 23, 27, 35, 50, 75, 125]
		ts, xs = traiectorie(f_A, x0, 25.0)
		lines!(ax6b, ts, xs, color = x0 < a_A ? rosu : verde, linewidth = 2.5)
	end
	fig6
end

# ╔═╡ 83a3c927-ff48-4453-a2dc-273e6757bf89
md"""
Pragul ``a`` separă două **bazine de atracție**: ``(0, a)`` duce la extincție, ``(a, \infty)`` la capacitatea de suport. Explorați modelul modificînd parametrii și populația inițială (încercați și ``x_0 = a`` exact):

* **Pragul Allee ``a``:** $(@bind a_allee PlutoUI.Slider(5.0:5.0:60.0, default = 20.0, show_value = true))
* **Capacitatea de suport ``K``:** $(@bind K_allee PlutoUI.Slider(70.0:10.0:150.0, default = 100.0, show_value = true))
* **Rata de creștere ``r``:** $(@bind r_allee PlutoUI.Slider(0.1:0.1:1.0, default = 0.5, show_value = true))
* **Populația inițială ``x_0``:** $(@bind x0_allee PlutoUI.Slider(0.0:0.5:150.0, default = 25.0, show_value = true))
"""

# ╔═╡ 927f7ac6-2fb5-4f92-8619-8f6cd49fe3e4
begin
	g_allee(x) = r_allee * x * (1 - x / K_allee) * (x / a_allee - 1)
	T_allee = 12 / r_allee
	fig7 = figura(size = (1100, 440))
	ax7a = Axis(fig7[1, 1]; stil_ax..., xlabel = "x", ylabel = "x' = f(x)",
		title = @sprintf("a = %g,  K = %g,  r = %g", a_allee, K_allee, r_allee))
	linie_de_faza!(ax7a, g_allee, [0.0, a_allee, K_allee], 0.0, 1.25K_allee, culoare = portocaliu)
	scatter!(ax7a, [x0_allee], [g_allee(x0_allee)], color = galben, markersize = 14,
		marker = :diamond, strokewidth = 1.5, strokecolor = :white)
	ax7b = Axis(fig7[1, 2]; stil_ax..., xlabel = "t", ylabel = "x(t)",
		title = @sprintf("Soluția pentru x₀ = %g", x0_allee))
	hlines!(ax7b, [a_allee, K_allee], color = :grey60, linestyle = :dash)
	for x0 in range(0, 1.25K_allee, length = 16)[2:end]
		ts, xs = traiectorie(g_allee, x0, T_allee)
		lines!(ax7b, ts, xs, color = (x0 < a_allee ? rosu : verde, 0.3), linewidth = 1.5)
	end
	ts7, xs7 = traiectorie(g_allee, x0_allee, T_allee)
	lines!(ax7b, ts7, xs7, color = galben, linewidth = 4)
	fig7
end

# ╔═╡ df52c183-787f-4177-9f0b-6be4f079db0d
md"""
**Ce observăm:**
- ``r`` schimbă doar **viteza** (scara de timp), nu și echilibrele sau stabilitatea lor;
- cu cît ``x_0`` este mai aproape de ``a``, cu atît populația „ezită” mai mult înainte de a se decide (``f'(a)`` este mic);
- o populație aflată deasupra pragului poate fi împinsă sub el de o perturbație (vînătoare, distrugerea habitatului) și apoi dispare **chiar dacă perturbația încetează**.

> 🌍 Consecință practică: pentru speciile cu efect Allee, programele de conservare sau de reintroducere trebuie să pornească de la populații **mai mari decît pragul** ``a``. În Cursul 6 vom vedea cum, modificînd un parametru, echilibrele se pot ciocni și pot dispărea (**bifurcații**).
"""

# ╔═╡ b932a091-5b08-40b0-a8e2-f833e1620b27
md"""
## 📝 7. Recapitulare și exerciții 

### 🧠 Recapitulare
- ↩️ **Integrarea** este operația inversă derivării; primitiva este determinată pînă la o constantă.
- 📏 **Sumele Riemann** (metoda Euler pentru ``X' = g(t)``) aproximează **aria de sub grafic**; eroarea este proporțională cu ``\Delta t``.
- 🎓 **Teorema fundamentală:** ``\int_a^b f(t)\,dt = F(b) - F(a)``; valoarea medie este ``\frac{1}{b-a}\int_a^b f(t)\,dt``.
- 📈 ``x' = rx \implies x(t) = x_0e^{rt}``; timpul de dublare sau de înjumătățire este ``\ln 2/\lvert r\rvert``.
- ⚖️ **Echilibre:** ``f(x^*) = 0``. **Stabilitate:** linia de fază (semnul lui ``f``) sau semnul lui ``f'(x^*)``.
- 🌱 **Logistic:** ``0`` instabil, ``K`` stabil. 🐋 **Allee:** ``0`` și ``K`` stabile, pragul ``a`` instabil ⟹ bistabilitate.

---
### ✍️ Exerciții
**Exercițiul 1 (integrare numerică).** Funcția ``g(x) = e^{-x^2}`` nu are o primitivă elementară.

a) Folosiți sume Riemann (funcția `suma_riemann` de mai sus) pentru a aproxima ``\int_0^2 e^{-x^2}\,dx`` cu ``n = 10, 100, 1000``.

b) Comparați cu valoarea obținută cu `quadgk`. Cum scade eroarea?
"""

# ╔═╡ 5164df8b-4dc4-41d3-a3bb-bcb6e6b666c6
# EXERCIȚIUL 1: completați codul
# g(x) = ...
# [suma_riemann(g, 0, 2, n) for n in (10, 100, 1000)]
# quadgk(g, 0, 2)

# ╔═╡ 7ff35f19-e747-490f-923a-6e87b0649ab5
md"""
**Exercițiul 2 (descompunere exponențială).** Un izotop radioactiv are timpul de înjumătățire de ``5`` ani.

a) Scrieți ecuația ``N' = -kN`` și găsiți constanta ``k``.

b) Dacă avem inițial ``100`` g de substanță, după cît timp rămîn doar ``10`` g?
"""

# ╔═╡ a8cd80ae-d5b0-4782-8bff-8d329cbdc0a3
hint(md"``e^{-5k} = 1/2`` dă ``k = \ln 2 / 5``. Apoi rezolvați ``100\,e^{-kt} = 10``, adică ``t = \ln 10 / k``.")

# ╔═╡ ed2f3b66-9206-47fe-b4b0-8f4343d910b5
md"""
**Exercițiul 3 (stabilitate).** Fie ecuația ``x' = x^3 - 4x``.

a) Găsiți toate punctele de echilibru.

b) Determinați stabilitatea lor cu criteriul derivatei.

c) Desenați linia de fază pe hîrtie și verificați-o cu funcția `linie_de_faza!`.
"""

# ╔═╡ ec0729a3-68c2-4b05-b097-b726efd53c3d
# EXERCIȚIUL 3: verificați cu Julia
# f_ex3(x) = x^3 - 4x
# fig_ex3 = figura(); ax_ex3 = Axis(fig_ex3[1, 1]; stil_ax...)
# linie_de_faza!(ax_ex3, f_ex3, [ ... ], -3, 3); fig_ex3

# ╔═╡ a8ff171a-6246-4f1d-85af-80550f58f2cd
md"""
**Exercițiul 4 (pescuit cu cotă fixă — provocare).** Acum se pescuiește o cantitate fixă ``H`` pe an:
```math
x' = 0.5\,x\left(1 - \frac{x}{100}\right) - H .
```
a) Găsiți punctele de echilibru în funcție de ``H``.

b) Ce se întîmplă cu populația pentru ``H = 10``? Dar pentru ``H = 15``?

c) Care este cota maximă ``H`` pentru care populația poate supraviețui? Comparați-o cu maximul lui ``0.5\,x(1 - x/100)`` din secțiunea 4. Dispariția bruscă a echilibrelor cînd ``H`` crește este un exemplu de **bifurcație** (Cursul 6).
"""

# ╔═╡ 766543b4-e275-42d8-9802-ad07eecadb58
hint(md"Echilibrele sînt rădăcinile ecuației de gradul doi ``0.005x^2 - 0.5x + H = 0``; ele există doar dacă discriminantul ``0.25 - 0.02H \ge 0``.")

# ╔═╡ Cell order:
# ╟─5597b08d-2218-4f9d-b46b-228c1161f894
# ╟─a17d6370-a346-4630-b5a3-0f07fbf5474a
# ╟─5fe3884c-8911-4bf6-b88d-a5cb6e8ee05c
# ╟─60547c1c-e309-4606-adca-f00e7b9c95f3
# ╟─cd5e96c4-a46b-4ba0-8254-abe5a518a50d
# ╟─ec53dd46-1054-4dda-b83a-2ca3070e6edd
# ╟─72400236-1ff8-426f-92e6-0647df943bb4
# ╟─6b69614e-1a69-4aa3-a0f7-79bfb85b32c7
# ╟─fb034a9e-44d1-4b25-b25a-245a62423cf4
# ╟─2ebd78a7-8a94-4077-812d-86e308566903
# ╟─78cd9df2-8b48-41cd-8d9f-4ee9979bcde3
# ╟─a26b0c6f-f66b-427c-b06f-eb182529ec47
# ╟─e0da1c43-41c2-475e-a898-503424d2d511
# ╟─10a166b2-1c43-4a69-879b-190cf49323d6
# ╟─4ec95b87-81c2-40c6-b3cf-2ba63c805ec1
# ╟─6af67b03-6a21-4113-9578-a69cc795db98
# ╟─e82ffe61-be7b-4888-adaa-fa92dccf6f1d
# ╟─aa48499a-994d-4baa-9ea7-68327f17b91e
# ╟─82e1c110-a244-4d5c-92a0-a0cbc26721e6
# ╟─9d53abcc-1c11-42a5-9f7c-51d660d232b5
# ╟─a721c2c7-b8c2-4fd9-842f-bbc69d6fe2d6
# ╟─762f11dd-a697-4a63-9e26-559b4e0148f7
# ╟─a71c9839-5b58-4c85-939e-de9c1053135e
# ╟─05920e3c-c24d-45b0-a65a-c5c6df52dd4c
# ╟─9b1be66f-ea19-41c6-8235-e90fa1257fd8
# ╟─30c59120-bee5-4743-8606-38ba3d9c1d43
# ╟─8a77a20a-92cb-46a3-aba5-00affcf583a4
# ╟─f231db04-5759-4209-b061-4d34f7f060a7
# ╟─37602046-94fa-4853-a4cb-ffe93f6841a2
# ╟─2c81e4da-6b65-4831-a7bc-e7cae78a78db
# ╟─70a6b506-4529-4ce3-b79e-6067759a8ff7
# ╟─d7fac950-643d-424a-b783-f774ac30c532
# ╟─e07d6e57-38cb-493d-9aba-7a957c7660e6
# ╟─8f034303-c135-4b30-8cf4-0f63c8f9526c
# ╟─83a3c927-ff48-4453-a2dc-273e6757bf89
# ╟─927f7ac6-2fb5-4f92-8619-8f6cd49fe3e4
# ╟─df52c183-787f-4177-9f0b-6be4f079db0d
# ╟─b932a091-5b08-40b0-a8e2-f833e1620b27
# ╠═5164df8b-4dc4-41d3-a3bb-bcb6e6b666c6
# ╟─7ff35f19-e747-490f-923a-6e87b0649ab5
# ╟─a8cd80ae-d5b0-4782-8bff-8d329cbdc0a3
# ╟─ed2f3b66-9206-47fe-b4b0-8f4343d910b5
# ╠═ec0729a3-68c2-4b05-b097-b726efd53c3d
# ╟─a8ff171a-6246-4f1d-85af-80550f58f2cd
# ╟─766543b4-e275-42d8-9802-ad07eecadb58
