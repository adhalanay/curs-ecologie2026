### A Pluto.jl notebook ###
# v1.0.3

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
	using DifferentialEquations
	using PlutoTeachingTools
	using Roots
	using DynamicalSystems
	using LinearAlgebra
	TableOfContents()
end

# ╔═╡ daaa48f8-1afa-4828-8c41-aa2134dcd652
WidthOverDocs()

# ╔═╡ b2c3d4e5-f6a7-8901-bcde-f12345678901
md"""
# 🌌 Cursul 5: Sisteme Dinamice în Dimensiuni Superioare și Modele Discrete

Până acum am studiat ecuații diferențiale și sisteme discrete scalare (1D). Astăzi vom face un pas înainte și vom explora cum se comportă sistemele când avem **două sau mai multe variabile** care interacționează, precum și cum modelăm sistemele în **timp discret** (pași de timp).

---

## 📋 Planul Cursului (Outline)
Pentru a ne orienta, iată parcursul nostru de astăzi:
1. **Echilibre în 2D**: De la sisteme liniare (surse, scurgeri, puncte șa) la spirale și centre.
2. **Sisteme Neliniare și Null-cline**: Cum analizăm grafic competiția între specii.
3. **Bazine de Atracție și Comutatoare Biologice**: Modelul operonului lac.
4. **Sisteme Discrete și Modele Leslie**: Dinamica populațiilor structurate pe vârste.
5. **Rezumat și Exerciții Practice**.

---
"""

# ╔═╡ c3d4e5f6-a7b8-9012-cdef-123456789012
md"""
## 1. Echilibre în 2D: Sisteme Liniare

În 1D, soluțiile fie cresc la infinit, fie se stabilizează. În 2D, geometria planului ne oferă o bogăție de comportamente noi. 
Un **punct de echilibru** este punctul $(x^*, y^*)$ unde sistemul "se oprește" (derivatele sunt zero).

Să explorăm cele mai importante tipuri de echilibre liniare folosind **câmpuri de vectori** (streamplots). Săgețile ne arată direcția în care evoluează sistemul din orice punct de pornire.
"""

# ╔═╡ d4e5f6a7-b8c9-0123-defa-234567890123
begin
	# Definim o grilă comună pentru vizualizări
	grid_range = LinRange(-2.5, 2.5, 25)
	
	# Funcția pentru desenarea câmpului de vectori
	function plot_vector_field(field_func, title_text; limits=(-2.5, 2.5, -2.5, 2.5))
		fig = Figure(size=(500, 400))
		ax = Axis(fig[1, 1], 
			title=title_text, 
			xlabel="x", ylabel="y", 
			limits=limits,
			aspect=DataAspect()
		)
		streamplot!(ax, field_func, grid_range, grid_range, 
			colormap=:magma, linewidth=1, arrow_size=15)
		return fig
	end
end

# ╔═╡ e5f6a7b8-c9d0-1234-efab-345678901234
md"""
### 🔴 Sursa (Nod Instabil)
Dacă ambele variabile cresc independent ($x'=x, y'=y$), orice punct se îndepărtează de origine.
"""

# ╔═╡ f6a7b8c9-d0e1-2345-fabc-456789012345
plot_vector_field((x,y) -> Point2f(x, y), "Sursă: x'=x, y'=y")

# ╔═╡ a7b8c9d0-e1f2-3456-abcd-567890123456
md"""
### 🔵 Scurgerea (Nod Stabil)
Dacă ambele variabile descreasc ($x'=-x, y'=-y$), toate traiectoriile sunt "aspirate" spre origine.
"""

# ╔═╡ b8c9d0e1-f2a3-4567-bcde-678901234567
plot_vector_field(u -> Point2f(-u[1], -u[2]), "Scurgere: x'=-x, y'=-y")

# ╔═╡ c9d0e1f2-a3b4-5678-cdef-789012345678
md"""
### 🟣 Punctul Șa (Saddle Point)
Ce se întâmplă dacă combinăm o direcție stabilă cu una instabilă? ($x'=x, y'=-y$). 
Sistemul este atras pe o axă, dar respins pe cealaltă. Este un echilibru instabil, dar cu o structură geometrică fascinantă.
"""

# ╔═╡ d0e1f2a3-b4c5-6789-defa-890123456789
plot_vector_field(u -> Point2f(u[1], -u[2]), "Punct Șa: x'=x, y'=-y")

# ╔═╡ e1f2a3b4-c5d6-7890-efab-901234567890
md"""
### 🌀 Spirale și Centre
Când variabilele interacționează (de ex. poziția și viteza unui arc), apar mișcări de rotație.
*   **Spirală Stabilă**: Arc cu frecare (energia se disipă).
*   **Centru**: Arc fără frecare (oscilație perpetuă).
"""

# ╔═╡ f2a3b4c5-d6e7-8901-fabc-012345678901
begin
	fig_spirals = Figure(size=(800, 400))
	
	# Spirala Stabilă (x'=y, y'=-x-y)
	ax1 = Axis(fig_spirals[1, 1], title="Spirală Stabilă (Arc cu frecare)", xlabel="x", ylabel="y", aspect=DataAspect())
	streamplot!(ax1, u -> Point2f(u[2], -u[1]-u[2]), grid_range, grid_range, colormap=:magma, linewidth=1, arrow_size=10)
	
	# Centru (x'=y, y'=-x)
	ax2 = Axis(fig_spirals[1, 2], title="Centru (Arc fără frecare)", xlabel="x", ylabel="y", aspect=DataAspect())
	streamplot!(ax2, u -> Point2f(u[2], -u[1]), grid_range, grid_range, colormap=:magma, linewidth=1, arrow_size=10)
	
	fig_spirals
end

# ╔═╡ a3b4c5d6-e7f8-9012-abcd-123456789012
md"""
---
## 2. Sisteme Neliniare și Metoda Null-Clinelor

În natură, sistemele sunt rar pur liniare. Să analizăm un model de **competiție** între două specii: Căprioare ($D$) și Elani ($M$).

```math 
D' = 3D - D^2 - DM 
```
```math
M' = 2M - M^2 - 0.5DM 
```

Cum găsim echilibrele și stabilitatea lor fără să rezolvăm ecuațiile analitic? Folosim **Null-clinele**!
*   **D-null-clina**: Curba unde $D' = 0$ (vectorii sunt verticali).
*   **M-null-clina**: Curba unde $M' = 0$ (vectorii sunt orizontali).
Intersecțiile lor sunt **punctele de echilibru**.
"""

# ╔═╡ b4c5d6e7-f8a9-0123-bcde-234567890123
begin
	fig_competition = Figure(size=(600, 500))
	ax_comp = Axis(fig_competition[1, 1], 
		title="Competiție Căprioare vs Elani", 
		xlabel="Populația D (Căprioare)", 
		ylabel="Populația M (Elani)", 
		limits=(-0.1, 5, -0.1, 5),
		aspect=DataAspect()
	)
	
	# Câmpul de vectori
	f_comp(u) = Point2f(3u[1]-u[1]^2-u[1]*u[2], 2u[2]-u[2]^2-0.5*u[1]*u[2])
	streamplot!(ax_comp, f_comp, LinRange(0, 5, 25), LinRange(0, 5, 25), colormap=:magma, linewidth=1, arrow_size=10)
	
	# Null-clinele
	# D'=0 => D(3-D-M)=0 => D=0 sau M=3-D
	# M'=0 => M(2-M-0.5D)=0 => M=0 sau M=2-0.5D
	x_vals = LinRange(0, 5, 100)
	lines!(ax_comp, x_vals, 3 .- x_vals, color=:blue, linewidth=3, label="D-null-clina (M=3-D)")
	lines!(ax_comp, x_vals, 2 .- 0.5.*x_vals, color=:red, linewidth=3, label="M-null-clina (M=2-0.5D)")
	
	# Punctele de echilibru
	equilibria = [(0,0), (3,0), (0,2), (2,1)]
	for eq in equilibria
		scatter!(ax_comp, [eq[1]], [eq[2]], markersize=15, color=:yellow, strokecolor=:black, strokewidth=2)
	end
	
	axislegend(ax_comp, position=:lt)
	fig_competition
end

# ╔═╡ c5d6e7f8-a9b0-1234-cdef-345678901234
md"""
🔍 **Ce observăm pe grafic?**
*   Null-clinele împart planul în regiuni. În fiecare regiune, vectorii au o direcție generală constantă.
*   Punctul $(2,1)$ este un **nod stabil** (coexistență). Toate săgețile din jurul lui converg spre el.
*   Punctele $(3,0)$ și $(0,2)$ sunt **puncte șa** (instabile).
*   Originea $(0,0)$ este o **sursă** (instabilă).
"""

# ╔═╡ d6e7f8a9-b0c1-2345-defa-456789012345
md"""
---
## 3. Bazine de Atracție și Comutatoare Biologice

Uneori, un sistem are **mai multe echilibre stabile**. Regiunea din spațiul fazelor care "trimite" traiectoriile către un anumit echilibru se numește **bazin de atracție**.

### Exemplu: Operonul Lac (Comutatorul Genetic)
Bacteria *E. coli* produce permează pentru a importa lactoză. Producția este costisitoare, deci are loc doar dacă lactoza depășește un prag.
Modelul matematic (după gruparea termenilor):

```math
x' = \frac{a+x^2}{1+x^2} - kx 
```
Unde $f(x) = \frac{a+x^2}{1+x^2}$ (producția sigmoidală) și $g(x) = kx$ (degradarea liniară).
"""

# ╔═╡ e7f8a9b0-c1d2-3456-efab-567890123456
begin
	a_param = 0.01
	k_param = 0.4
	f_lac(x) = (a_param + x^2) / (1 + x^2)
	g_lac(x) = k_param * x

	fig_lac = Figure(size=(600, 400))
	ax_lac = Axis(fig_lac[1, 1], 
		title="Metoda Grafică pentru Echilibre (Operonul Lac)", 
		xlabel="Concentrația lactoză (x)", 
		ylabel="Rata de schimbare",
		limits=(0, 2.5, 0, 1.2)
	)
	
	x_plot = LinRange(0, 2.5, 200)
	lines!(ax_lac, x_plot, f_lac.(x_plot), color=:blue, linewidth=3, label="Producție: f(x)")
	lines!(ax_lac, x_plot, g_lac.(x_plot), color=:red, linewidth=3, label="Degradare: g(x)")
	
	# Găsirea rădăcinilor (intersecțiilor)
	h_lac(x) = f_lac(x) - g_lac(x)
	sol1 = find_zero(h_lac, (0, 0.2))
	sol2 = find_zero(h_lac, (0.3, 0.6))
	sol3 = find_zero(h_lac, (1, 3))
	
	sols = [sol1, sol2, sol3]
	scatter!(ax_lac, sols, f_lac.(sols), markersize=12, color=:black, label="Echilibre")
	
	axislegend(ax_lac, position=:rb)
	fig_lac
end

# ╔═╡ f8a9b0c1-d2e3-4567-fabc-678901234567
md"""
🔍 **Interpretare Biologică:**
*   **Echilibrele stabile** (capetele): Celula este fie "oprită" (fără permează), fie "pornită" (permează maximă).
*   **Echilibrul instabil** (mijlocul): Acționează ca un **prag (întrerupător)**. Dacă lactoza trece de acest punct, sistemul "sare" ireversibil spre starea activă.
"""

# ╔═╡ a9b0c1d2-e3f4-5678-abcd-789012345678
md"""
---
## 4. Sisteme Discrete Liniare și Explorare Interactivă

În loc de timp continuu ($t$), să considerăm timp discret ($n = 0, 1, 2...$), specific populațiilor cu reproducere anuală.

```math 
x_{n+1} = a x_n + b y_n ```
```math
y_{n+1} = c x_n + d y_n 
```

👇 **Folosiți sliders de mai jos pentru a modifica parametrii matricei și observați cum se schimbă stabilitatea sistemului!**
"""

# ╔═╡ b0c1d2e3-f4a5-6789-bcde-890123456789
md"""
*   Parametrul **a** = $(@bind a NumberField(0:0.1:1.5, default=1.1))
*   Parametrul **b** = $(@bind b NumberField(0:0.1:1.5, default=0.0))
*   Parametrul **c** = $(@bind c NumberField(0:0.1:1.5, default=0.0))
*   Parametrul **d** = $(@bind d NumberField(0:0.1:1.5, default=0.9))
*   Numărul de pași = $(@bind timp PlutoUI.Slider(10:200, default=50, show_value=true))
"""

# ╔═╡ c1d2e3f4-a5b6-7890-cdef-901234567890
begin
	function gen_discrete(u, p, n)
		x, y = u
		a, b, c, d = p
		return SVector(a*x + b*y, c*x + d*y)
	end
	
	u0 = [30.0, 80.0]
	p0 = [a, b, c, d]
	sislin = DeterministicIteratedMap(gen_discrete, u0, p0)
	X, t = trajectory(sislin, timp)
	
	fig_discrete = Figure(size=(900, 400))
	
	# Traiectoria în spațiul fazelor
	ax_disc1 = Axis(fig_discrete[1, 1], title="Spațiul Fazelor (x vs y)", xlabel="x", ylabel="y")
	scatter!(ax_disc1, u0[1], u0[2], markersize=15, color=:red, label="Start")
	lines!(ax_disc1, X[:,1], X[:,2], color=:blue, linewidth=2)
	scatter!(ax_disc1, X[:,1], X[:,2], markersize=5, color=:blue)
	axislegend(ax_disc1)
	
	# Seriile temporale
	ax_disc2 = Axis(fig_discrete[1, 2], title="Seriile Temporale", xlabel="Timp (n)", ylabel="Valoare")
	scatter!(ax_disc2, t, X[:,1], color=:blue, label="x_n")
	scatter!(ax_disc2, t, X[:,2], color=:red, label="y_n")
	axislegend(ax_disc2)
	
	fig_discrete
end

# ╔═╡ d2e3f4a5-b6c7-8901-defa-012345678901
md"""
---
## 5. Modele Leslie (Dinamica Populațiilor Structurate)

Un caz special de sisteme discrete liniare este **Modelul Leslie**, folosit în biologie pentru a modela populații împărțite pe clase de vârstă (ex: juvenili și adulți).

```math 
M = \begin{pmatrix} 0.65 & 0.5 \\ 0.25 & 0.9 \end{pmatrix}
```

*   **Diagonala**: Procentul care rămâne în aceeași clasă (supraviețuire).
*   **În afara diagonalei**: Tranziții între clase (creștere) și natalitate.

Dacă valorile proprii ale matricei sunt subunitare, populația se stinge. Dacă sunt supraunitare, populația crește exponențial.
"""

# ╔═╡ e3f4a5b6-c7d8-9012-efab-123456789012
md"""
---
## 🎯 Rezumat (Takeaways)

1.  **Geometria 2D**: Sistemele liniare în 2D pot forma surse, scurgeri, puncte șa, spirale sau centre, în funcție de interacțiunea dintre variabile.
2.  **Puterea Null-clinelor**: Pentru sistemele neliniare, intersecțiile null-clinelor ($x'=0$ și $y'=0$) ne oferă echilibrele, iar regiunile dintre ele ne arată direcția de evoluție.
3.  **Comutatoare Biologice**: Sistemele neliniare pot avea mai multe bazine de atracție, permițând celulelor să funcționeze ca niște "întrerupătoare" cu prag (ex: operonul lac).
4.  **Timp Discret**: Modelele matriciale (Leslie) sunt esențiale pentru a prezice viitorul populațiilor structurate pe vârste, pe baza ratelor de supraviețuire și natalitate.

---
"""

# ╔═╡ f4a5b6c7-d8e9-0123-fabc-234567890123
md"""
## 📝 Exerciții Practice

Încercați să rezolvați următoarele probleme pentru a vă consolida cunoștințele. Puteți scrie cod Julia în celule noi pentru a verifica!

**Exercițiul 1: Clasificarea Echilibrelor**
Avem sistemul liniar:
$x' = 2x + y$
$y' = x + 2y$

Scrieți matricea sistemului. Care sunt valorile proprii? Ce tip de echilibru reprezintă originea (sursă, scurgere, șa, spirală)?

**Exercițiul 2: Găsirea Null-clinelor**
Pentru sistemul de competiție:
$x' = x(1 - x - y)$
$y' = y(0.5 - 0.5x - y)$

Care sunt ecuațiile pentru x-null-clina și y-null-clina? Găsiți toate cele 4 puncte de echilibru.

**Exercițiul 3: Model Leslie**
Avem o populație de broaște cu două stadii (mormoloci și adulți). Matricea este:
$M = \begin{pmatrix} 0 & 2 \\ 0.5 & 0.2 \end{pmatrix}$

Ce semnificație biologică are valoarea `2` din matrice? Ce semnificație are valoarea `0.5`? Dacă începem cu 100 de mormoloci și 0 adulți, cum va arăta populația după mulți ani?
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-e5f6-7890-abcd-ef1234567890
# ╠═daaa48f8-1afa-4828-8c41-aa2134dcd652
# ╟─b2c3d4e5-f6a7-8901-bcde-f12345678901
# ╟─c3d4e5f6-a7b8-9012-cdef-123456789012
# ╠═d4e5f6a7-b8c9-0123-defa-234567890123
# ╟─e5f6a7b8-c9d0-1234-efab-345678901234
# ╠═f6a7b8c9-d0e1-2345-fabc-456789012345
# ╟─a7b8c9d0-e1f2-3456-abcd-567890123456
# ╠═b8c9d0e1-f2a3-4567-bcde-678901234567
# ╟─c9d0e1f2-a3b4-5678-cdef-789012345678
# ╠═d0e1f2a3-b4c5-6789-defa-890123456789
# ╟─e1f2a3b4-c5d6-7890-efab-901234567890
# ╟─f2a3b4c5-d6e7-8901-fabc-012345678901
# ╟─a3b4c5d6-e7f8-9012-abcd-123456789012
# ╟─b4c5d6e7-f8a9-0123-bcde-234567890123
# ╟─c5d6e7f8-a9b0-1234-cdef-345678901234
# ╟─d6e7f8a9-b0c1-2345-defa-456789012345
# ╠═e7f8a9b0-c1d2-3456-efab-567890123456
# ╟─f8a9b0c1-d2e3-4567-fabc-678901234567
# ╟─a9b0c1d2-e3f4-5678-abcd-789012345678
# ╟─b0c1d2e3-f4a5-6789-bcde-890123456789
# ╟─c1d2e3f4-a5b6-7890-cdef-901234567890
# ╟─d2e3f4a5-b6c7-8901-defa-012345678901
# ╟─e3f4a5b6-c7d8-9012-efab-123456789012
# ╟─f4a5b6c7-d8e9-0123-fabc-234567890123
