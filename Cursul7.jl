### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> title = "Cursul VII"
#> description = "Algebră liniară și modele matriciale"
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
	using Makie, CairoMakie
	using PlutoUI
	using PlutoTeachingTools
	using LinearAlgebra
	using Printf

	CairoMakie.activate!(type = "png", px_per_unit = 2)
	set_theme!(theme_minimal())
	TableOfContents()
end

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ ec01674b-79a9-4a71-9a12-2e04384599da
begin
	# stările x₀, x₁, …, xₙ ale modelului x_{k+1} = M x_k
	function iteratii(M, x0, n)
		stari = [float.(x0)]
		for _ in 1:n
			push!(stari, M * stari[end])
		end
		stari
	end

	# valoarea proprie dominantă λ₁ (cea de modul maxim) și vectorul propriu corespunzător,
	# normalizat astfel încît suma componentelor să fie 1 (distribuția stabilă pe stadii)
	function dominanta(M)
		valori, vectori = eigen(M)
		i = argmax(abs.(valori))
		v = real.(vectori[:, i])
		real(valori[i]), v ./ sum(v)
	end

	# graful ciclului de viață al unui model pe stadii cu matricea M: bucle pentru indivizii
	# care rămîn în stadiu, săgeți spre stadiul următor și, cu roșu, săgețile de reproducere
	# (primul rînd al matricei), care se întorc la primul stadiu
	function graf_ciclu(M, etichete; id = "g")
		n = size(M, 1)
		num(x) = @sprintf("%g", x)
		dx, r, y = n <= 3 ? 170 : 115, 22, 85
		xs = [45 + dx * (k - 1) for k in 1:n]
		adinc(k) = 36 + 22 * (k - 1)        # adîncimea arcului de reproducere din stadiul k
		W, H = xs[end] + 45, y + r + 0.75 * adinc(n) + 20
		el = String[]
		for k in 1:n
			x = xs[k]
			if M[k, k] > 0
				push!(el, """<path d="M $(x-11) $(y-r+3) C $(x-38) $(y-r-50), $(x+38) $(y-r-50), $(x+11) $(y-r+3)" fill="none" stroke="currentColor" stroke-width="1.5" marker-end="url(#$id-n)"/>""")
				push!(el, """<text x="$x" y="$(y-r-42)" text-anchor="middle">$(num(M[k, k]))</text>""")
			end
			if k < n && M[k+1, k] > 0
				push!(el, """<line x1="$(x+r)" y1="$y" x2="$(xs[k+1]-r-2)" y2="$y" stroke="currentColor" stroke-width="1.5" marker-end="url(#$id-n)"/>""")
				push!(el, """<text x="$(x+dx/2)" y="$(y-7)" text-anchor="middle">$(num(M[k+1, k]))</text>""")
			end
			if k > 1 && M[1, k] > 0
				h = adinc(k)
				push!(el, """<path d="M $x $(y+r) C $x $(y+r+h), $(xs[1]) $(y+r+h), $(xs[1]) $(y+r+3)" fill="none" stroke="#d1495b" stroke-width="1.5" marker-end="url(#$id-r)"/>""")
				push!(el, """<text x="$(x+7)" y="$(y+r+18)" text-anchor="start" fill="#d1495b">$(num(M[1, k]))</text>""")
			end
			push!(el, """<circle cx="$x" cy="$y" r="$r" fill="none" stroke="currentColor" stroke-width="2"/>""")
			push!(el, """<text x="$x" y="$(y+5)" text-anchor="middle" font-weight="bold">$(etichete[k])</text>""")
		end
		varf(cid, culoare) = """<marker id="$cid" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="$culoare"/></marker>"""
		HTML("""<svg width="$W" height="$H" viewBox="0 0 $W $H" style="font: 14px sans-serif; fill: currentColor; max-width: 100%">
		<defs>$(varf("$id-n", "currentColor"))$(varf("$id-r", "#d1495b"))</defs>
		$(join(el, "\n"))
		</svg>""")
	end
end

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 🧮 Cursul 7 — Algebră liniară și modele matriciale

**Durata estimată:** 100 minute
"""

# ╔═╡ 186cde3c-94c6-46cf-86b7-bc61461674d2
md"""
## 🎯 Obiectivele lecției
1. Să descriem starea unei populații structurate printr-un **vector** și evoluția ei printr-o **aplicație liniară** și **matricea** acesteia.
2. Să construim **modele matriciale de populație** (Leslie și pe stadii) pornind de la **graful ciclului de viață** și să compunem ani diferiți prin **produsul matricelor**.
3. Să calculăm **valorile și vectorii proprii** ai unei matrice ``2 \times 2`` cu ajutorul **ecuației caracteristice**.
4. Să explicăm comportamentul pe termen lung prin **valoarea proprie dominantă** (rata de creștere) și **distribuția stabilă pe stadii**, iar pe termen scurt prin **regimul tranzitoriu**.
5. Să folosim un model matricial pentru o decizie de **conservare**: țestoasa marină *Caretta caretta*.
"""

# ╔═╡ 644a8297-b028-4515-b29f-0600e67ac8f9
md"""
## 🔁 0. Ce știm deja (5 min)
* **Cursul 4:** modelul de creștere geometrică ``x_{n+1} = \lambda\,x_n`` are soluția ``x_n = \lambda^n x_0``: populația crește dacă ``\lambda > 1`` și scade spre ``0`` dacă ``0 < \lambda < 1``.
* **Cursul 5:** am împărțit femelele ursului negru american (*Ursus americanus*) în **juvenili** ``J`` și **adulți** ``A``. În fiecare an ``65\%`` dintre juvenili rămîn juvenili, ``25\%`` se maturizează, fiecare adult are în medie ``0.5`` pui femele, iar ``90\%`` dintre adulți supraviețuiesc:
```math
\begin{pmatrix} J_{n+1} \\ A_{n+1} \end{pmatrix} =
\begin{pmatrix} 0.65 & 0.5 \\ 0.25 & 0.9 \end{pmatrix}
\begin{pmatrix} J_n \\ A_n \end{pmatrix} .
```
Simulările au arătat că, după cîțiva ani, populația crește cu un factor ``\lambda \approx 1.15`` pe an, iar juvenilii ajung să fie ``50\%`` din populație, **oricare ar fi starea inițială**.

**Două întrebări pentru azi:**
1. 🧮 Cum aflăm factorul ``\lambda`` și structura finală a populației **direct din matrice**, fără să simulăm?
2. 🐢 Pentru o specie amenințată, ce stadiu de viață trebuie protejat ca declinul să se oprească: ouăle, juvenilii sau adulții?
"""

# ╔═╡ 6ece3812-853a-4c42-a863-d9430a9bdcb3
PlutoUI.LocalResource("./Black_Bear.jpg", :width => 600)

# ╔═╡ 647b474a-ec17-465b-8edb-c9f7ae631949
md"""
## 🧭 1. Vectori, aplicații liniare și matrice (15 min)

### 1.1. Starea populației este un vector
Starea unei populații cu două stadii este un **vector** din ``\mathbb{R}^2``, pe care îl scriem pe coloană:
```math
\mathbf{x} = \begin{pmatrix} J \\ A \end{pmatrix} .
```
Vectorii se adună și se înmulțesc cu numere (**scalari**) componentă cu componentă: dacă reunim două grupuri de urși, adunăm vectorii; dacă populația se dublează, înmulțim vectorul cu ``2``. O populație cu ``m`` stadii este un vector din ``\mathbb{R}^m``.

Doi vectori au un rol special: ``\mathbf{e}_1 = (1, 0)`` (un juvenil) și ``\mathbf{e}_2 = (0, 1)`` (un adult). Orice stare este o **combinație liniară** a lor:
```math
\begin{pmatrix} J \\ A \end{pmatrix} = J \begin{pmatrix} 1 \\ 0 \end{pmatrix} + A \begin{pmatrix} 0 \\ 1 \end{pmatrix} = J\,\mathbf{e}_1 + A\,\mathbf{e}_2 .
```
Mulțimea ``\{\mathbf{e}_1, \mathbf{e}_2\}`` este **baza canonică** a lui ``\mathbb{R}^2``. În general, o **bază** este o mulțime de vectori din care orice vector se obține, **în mod unic**, ca o combinație liniară. În secțiunea 3 vom găsi o bază mult mai potrivită pentru urși.

### 1.2. Aplicații liniare și matricea lor
O funcție ``f : \mathbb{R}^m \to \mathbb{R}^m`` este **liniară** dacă
```math
f(\mathbf{x} + \mathbf{y}) = f(\mathbf{x}) + f(\mathbf{y}), \qquad f(c\,\mathbf{x}) = c\,f(\mathbf{x}) .
```
Biologic, aceasta înseamnă că fiecare individ contribuie la populația de anul viitor **independent de ceilalți**: două grupuri de urși evoluează ca și cum ar fi separate, iar o populație de două ori mai mare produce de două ori mai mulți urși. Ipoteza este rezonabilă cînd populația este mult sub capacitatea de suport; altfel intervine competiția (modelul logistic, Cursurile 3 și 4).

O aplicație liniară este **complet determinată** de imaginile vectorilor bazei:
```math
f\begin{pmatrix} J \\ A \end{pmatrix} = f(J\,\mathbf{e}_1 + A\,\mathbf{e}_2) = J\,f(\mathbf{e}_1) + A\,f(\mathbf{e}_2) .
```
Așezăm ``f(\mathbf{e}_1)`` și ``f(\mathbf{e}_2)`` drept coloane într-un tabel: **matricea** lui ``f``. Produsul dintre o matrice și un vector este combinația coloanelor, cu coeficienții dați de componentele vectorului; fiecare componentă a rezultatului este un rînd al matricei înmulțit cu vectorul („rînd × coloană”):
```math
\begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}
\begin{pmatrix} x_1 \\ x_2 \end{pmatrix}
= x_1 \begin{pmatrix} a_{11} \\ a_{21} \end{pmatrix} + x_2 \begin{pmatrix} a_{12} \\ a_{22} \end{pmatrix}
= \begin{pmatrix} a_{11}x_1 + a_{12}x_2 \\ a_{21}x_1 + a_{22}x_2 \end{pmatrix} .
```

> 💡 **Cum citim matricea unui model de populație:** coloana ``j`` spune ce devine, după un an, **un singur individ** din stadiul ``j``. Pentru urs, un juvenil devine (în medie) ``0.65`` juvenili și ``0.25`` adulți, iar un adult devine ``0.5`` juvenili (puii lui) și ``0.9`` adulți (el însuși, dacă supraviețuiește).

Notăm cu ``L`` matricea ursului. Pentru ``100`` de juvenili și ``10`` adulți:
```math
L \begin{pmatrix} 100 \\ 10 \end{pmatrix} = 100\begin{pmatrix} 0.65 \\ 0.25 \end{pmatrix} + 10 \begin{pmatrix} 0.5 \\ 0.9 \end{pmatrix} = \begin{pmatrix} 70 \\ 34 \end{pmatrix} .
```
În Julia, matricele se scriu pe rînduri (elementele unui rînd separate prin spații, rîndurile prin `;`), iar produsul se scrie cu `*`:
"""

# ╔═╡ b187da58-91ab-4f2c-8955-6634c7c47dbe
let
	L = [0.65 0.5; 0.25 0.9]   # matricea ursului negru
	x = [100, 10]              # 100 de juvenili și 10 adulți
	L * x, 100 * L[:, 1] + 10 * L[:, 2]   # același rezultat, calculat în două feluri
end

# ╔═╡ ad470386-c9ca-4461-81e1-7fd11639b7b1
md"""
### ✅ Verificare rapidă
Calculați de mînă ``L \begin{pmatrix} 0 \\ 10 \end{pmatrix}`` și ``L \begin{pmatrix} 20 \\ 20 \end{pmatrix}``. Interpretați primul rezultat. Ce are special al doilea?
"""

# ╔═╡ 55b17dcd-ccfb-4c30-9cf2-d1d0341d44e5
hint(md"``L\,(0, 10) = (5, 9)``: zece adulți au în medie ``5`` pui femele, iar ``9`` dintre ei supraviețuiesc. ``L\,(20, 20) = (13 + 10,\ 5 + 18) = (23, 23) = 1.15 \cdot (20, 20)``: vectorul nu și-a schimbat direcția, ci doar s-a înmulțit cu ``1.15``, factorul de creștere din Cursul 5. Astfel de vectori sînt subiectul secțiunii 3.")

# ╔═╡ c025dcf9-37bd-43c9-a77f-540d33389e08
md"""
## 🔄 2. Modele matriciale de populație (20 min)

### 2.1. Graful ciclului de viață
Un model matricial se construiește pornind de la **graful ciclului de viață**: cîte un nod pentru fiecare stadiu și cîte o săgeată pentru fiecare mod în care un individ dintr-un stadiu contribuie, după un an, la un stadiu (eventual același). Săgeata de la stadiul ``j`` la stadiul ``i`` poartă valoarea ``l_{ij}``, elementul din rîndul ``i`` și coloana ``j`` al matricei. Pentru ursul negru (săgeata de reproducere este roșie):
"""

# ╔═╡ 7e9c12a7-ffc9-4874-8682-b713e14ac3dc
graf_ciclu([0.65 0.5; 0.25 0.9], ["J", "A"]; id = "urs")

# ╔═╡ 8a49308f-1ff7-4ad1-87c4-e25f8bb4d541
md"""
* săgețile care **pleacă** dintr-un stadiu formează coloana lui: supraviețuire, maturizare, reproducere;
* săgețile care **intră** într-un stadiu formează rîndul lui: din ce se compune stadiul în anul următor.

Două tipuri de modele sînt folosite des:
* **Modelul Leslie** (P. H. Leslie, 1945): populația este împărțită în **clase de vîrstă** de cîte un an. Cine supraviețuiește trece obligatoriu în clasa următoare, deci matricea are elemente nenule doar pe primul rînd (**fecunditățile** ``f_i``) și sub diagonală (**supraviețuirile** ``s_i``):
```math
L = \begin{pmatrix} f_1 & f_2 & f_3 \\ s_1 & 0 & 0 \\ 0 & s_2 & 0 \end{pmatrix} .
```
* **Modelul pe stadii** (L. P. Lefkovitch, 1965): clasele sînt **stadii de dezvoltare** (ou, larvă, juvenil, adult; sămînță, rozetă, plantă înflorită), iar un individ poate rămîne mai mulți ani în același stadiu, deci apar elemente nenule și pe diagonală. Modelul ursului este de acest tip. Pentru multe specii (plante, pești, țestoase) stadiul spune mai mult decît vîrsta: o plantă mică și bătrînă produce mai puține semințe decît una mare și tînără.

### ✅ Verificare rapidă
O plantă are trei stadii: **semințe** ``S``, **rozete** ``R`` și **plante înflorite** ``F``. În fiecare an ``10\%`` dintre semințe rămîn în sol fără să germineze, iar ``20\%`` germinează și devin rozete; ``50\%`` dintre rozete rămîn rozete, iar ``30\%`` înfloresc; o plantă înflorită produce ``40`` de semințe și apoi moare. Desenați graful ciclului de viață și scrieți matricea modelului.
"""

# ╔═╡ e1a18af2-f6ee-4945-a071-ae2138e86b49
hint(md"Cu stările în ordinea ``(S, R, F)``: ``L = \begin{pmatrix} 0.1 & 0 & 40 \\ 0.2 & 0.5 & 0 \\ 0 & 0.3 & 0 \end{pmatrix}``. Graful are bucle pe ``S`` (``0.1``) și pe ``R`` (``0.5``), săgețile ``S \to R`` (``0.2``) și ``R \to F`` (``0.3``) și săgeata de reproducere ``F \to S`` (``40``). Coloana lui ``F`` este ``(40, 0, 0)``: plantele înflorite nu supraviețuiesc.")

# ╔═╡ e3590799-a03c-441c-8493-254909f373c3
md"""
### 2.2. Ani buni și ani proști: produsul matricelor
Matricea ``L`` descrie un **an bun**. Într-un **an prost** (iarnă grea, hrană puțină) supraviețuirea și reproducerea scad:
```math
L = \begin{pmatrix} 0.65 & 0.5 \\ 0.25 & 0.9 \end{pmatrix} \ \text{(an bun)}, \qquad
M = \begin{pmatrix} 0.5 & 0.4 \\ 0.1 & 0.8 \end{pmatrix} \ \text{(an prost)} .
```
Dacă unui an bun îi urmează un an prost, după doi ani ``\mathbf{x}_2 = M\,\mathbf{x}_1 = M\,(L\,\mathbf{x}_0) = (ML)\,\mathbf{x}_0``. Compunerea a două aplicații liniare este tot o aplicație liniară, iar matricea ei este **produsul** ``ML``: elementul din rîndul ``i`` și coloana ``j`` al lui ``ML`` este rîndul ``i`` al lui ``M`` înmulțit cu coloana ``j`` a lui ``L``:
```math
ML = \begin{pmatrix} m_{11} & m_{12} \\ m_{21} & m_{22} \end{pmatrix}
\begin{pmatrix} l_{11} & l_{12} \\ l_{21} & l_{22} \end{pmatrix} =
\begin{pmatrix} m_{11}l_{11} + m_{12}l_{21} & m_{11}l_{12} + m_{12}l_{22} \\ m_{21}l_{11} + m_{22}l_{21} & m_{21}l_{12} + m_{22}l_{22} \end{pmatrix} .
```
Altfel spus, coloana ``j`` a lui ``ML`` este ``M`` aplicată coloanei ``j`` a lui ``L``: ce devine un individ din stadiul ``j`` după un an bun, urmărit încă un an prost.

> ⚠️ Spre deosebire de numere, **înmulțirea matricelor nu este comutativă**: în general ``ML \neq LM``. Matricea primului an se scrie **la dreapta**.
"""

# ╔═╡ 80504029-1e38-4b0f-a791-262e1991260e
let
	L = [0.65 0.5; 0.25 0.9]   # an bun
	M = [0.5 0.4; 0.1 0.8]     # an prost
	round.(M * L, digits = 3), round.(L * M, digits = 3)   # bun, apoi prost  vs.  prost, apoi bun
end

# ╔═╡ b76aecfd-8aa2-4f8f-8458-e04601a28830
let
	L = [0.65 0.5; 0.25 0.9]
	M = [0.5 0.4; 0.1 0.8]
	x0 = [100.0, 50.0]
	ani = 20
	ordine1 = [isodd(k) ? L : M for k in 1:ani]   # bun, prost, bun, ...
	ordine2 = [isodd(k) ? M : L for k in 1:ani]   # prost, bun, prost, ...
	traiectorie(ordine) = accumulate((x, P) -> P * x, ordine; init = x0)
	N1 = [sum(x0); sum.(traiectorie(ordine1))]
	N2 = [sum(x0); sum.(traiectorie(ordine2))]

	fig = Figure(size = (900, 430))
	ax = Axis(fig[1, 1], title = "Aceiași ani, în altă ordine",
		subtitle = @sprintf("după %d ani: %.1f urși (bun, prost, ...) față de %.1f urși (prost, bun, ...)",
			ani, N1[end], N2[end]),
		xlabel = "Anul (n)", ylabel = "Populația totală Jₙ + Aₙ")
	scatterlines!(ax, 0:ani, N1, color = :seagreen, linewidth = 2, label = "bun, prost, bun, ...")
	scatterlines!(ax, 0:ani, N2, color = :darkorange, linewidth = 2, linestyle = :dash,
		marker = :rect, label = "prost, bun, prost, ...")
	axislegend(ax, position = :lt)
	fig
end

# ╔═╡ 64bd52b4-2bd1-4633-85d6-a4869c80f87b
md"""
### 🤔 Ce am observat?
* Ambele succesiuni au același număr de ani buni și de ani proști, dar după un număr par de ani populațiile diferă: starea la un moment dat depinde de **ordinea** anilor, pentru că ``(ML)^k \mathbf{x}_0 \neq (LM)^k \mathbf{x}_0``.
* Totuși, la fiecare doi ani ambele populații se înmulțesc, pe termen lung, cu **același** factor, aproximativ ``1.035``: matricele ``ML`` și ``LM`` au aceleași valori proprii (Exercițiul 3). În secțiunea 3 vom vedea de ce tocmai valorile proprii decid comportamentul pe termen lung.
"""

# ╔═╡ d3e5263a-3fe4-4e6b-94c2-ab0de2d53b2b
md"""
## 🎯 3. Valori și vectori proprii (35 min)

### 3.1. Direcții care nu se rotesc
În Verificarea rapidă din secțiunea 1 am găsit ``L\,(20, 20) = 1.15\,(20, 20)``. În general, o matrice **rotește** și **alungește** vectorii. Căutăm vectorii nenuli ``\mathbf{v}`` pe care matricea doar îi alungește sau îi scurtează, fără să le schimbe direcția:
```math
L\,\mathbf{v} = \lambda\,\mathbf{v} .
```
Un astfel de vector ``\mathbf{v}`` se numește **vector propriu**, iar numărul ``\lambda`` **valoare proprie** a matricei ``L``. Dacă structura populației este dată de un vector propriu, ea nu se mai schimbă de la un an la altul, iar fiecare stadiu se înmulțește cu același factor ``\lambda``: exact comportamentul pe termen lung observat în Cursul 5.

### 🎛️ Experimentați
Vectorul albastru ``\mathbf{u}`` are lungimea ``1`` și face unghiul ``\theta`` cu axa juvenililor; vectorul roșu este imaginea lui, ``L\,\mathbf{u}``. Cînd ``\mathbf{u}`` parcurge cercul gri, ``L\,\mathbf{u}`` parcurge elipsa roșie. Căutați unghiurile pentru care cei doi vectori sînt pe aceeași dreaptă; în acest caz lungimea lui ``L\,\mathbf{u}`` este chiar valoarea proprie.

Unghiul ``\theta``: $(@bind θ_dir PlutoUI.Slider(0:1:180, show_value=true, default=0))

Arată direcțiile proprii: $(@bind arata_proprii PlutoUI.CheckBox(default=false))
"""

# ╔═╡ eb0af3c2-c015-496a-8e2d-112e8cf50f66
let
	L = [0.65 0.5; 0.25 0.9]
	u = [cosd(θ_dir), sind(θ_dir)]
	Lu = L * u
	unghi = acosd(clamp(u ⋅ Lu / norm(Lu), -1, 1))
	τ = range(0, 2π, length = 200)
	elipsa = [L * [cos(t), sin(t)] for t in τ]

	fig = Figure(size = (620, 600))
	ax = Axis(fig[1, 1], title = "Vectorul u și imaginea lui, L u",
		subtitle = @sprintf("unghiul dintre u și L u: %.1f°;  lungimea lui L u: %.3f", unghi, norm(Lu)),
		xlabel = "juvenili", ylabel = "adulți", aspect = DataAspect(), limits = (-1.5, 1.5, -1.5, 1.5))
	hlines!(ax, 0, color = :gray80)
	vlines!(ax, 0, color = :gray80)
	lines!(ax, cos.(τ), sin.(τ), color = :gray60, linewidth = 1.5)
	lines!(ax, first.(elipsa), last.(elipsa), color = (:crimson, 0.4), linewidth = 1.5)
	if arata_proprii
		for (λ, v) in ((1.15, [1, 1] / √2), (0.4, [-2, 1] / √5))
			lines!(ax, [-2v[1], 2v[1]], [-2v[2], 2v[2]], color = :purple, linestyle = :dash, linewidth = 2)
			text!(ax, 1.25v[1], 1.25v[2], text = "λ = $λ", color = :purple, offset = (8, 0), fontsize = 16)
		end
	end
	arrows2d!(ax, [0.0], [0.0], [u[1]], [u[2]], color = :royalblue,
		shaftwidth = 3, tipwidth = 14, tiplength = 12)
	arrows2d!(ax, [0.0], [0.0], [Lu[1]], [Lu[2]], color = :crimson,
		shaftwidth = 3, tipwidth = 14, tiplength = 12)
	text!(ax, u[1], u[2], text = "u", color = :royalblue, offset = (6, 6), fontsize = 18)
	text!(ax, Lu[1], Lu[2], text = "L u", color = :crimson, offset = (6, 6), fontsize = 18)
	fig
end

# ╔═╡ d5e51b29-2530-4e8d-a0cb-33ea3cb8586d
md"""
**Ce observăm:** există exact două direcții care nu se rotesc:
* ``\theta = 45^\circ``, direcția lui ``\mathbf{v}_1 = (1, 1)``: ``L\,\mathbf{v}_1 = 1.15\,\mathbf{v}_1``, vectorul se alungește;
* ``\theta \approx 153.4^\circ``, direcția lui ``\mathbf{v}_2 = (2, -1)``: ``L\,\mathbf{v}_2 = 0.4\,\mathbf{v}_2``, vectorul se scurtează.

Al doilea vector are o componentă negativă, deci nu reprezintă o populație reală; vom vedea însă că și el este necesar pentru a descrie evoluția populației.

### 3.2. Ecuația caracteristică
Cum găsim valorile proprii fără să le căutăm cu slider-ul? Pentru ``L = \begin{pmatrix} a & b \\ c & d \end{pmatrix}``, ecuația ``L\mathbf{v} = \lambda\mathbf{v}`` se scrie pe componente
```math
\begin{aligned}
(a - \lambda)\,v_1 + b\,v_2 &= 0 \\
c\,v_1 + (d - \lambda)\,v_2 &= 0 .
\end{aligned}
```
Sistemul are întotdeauna soluția ``v_1 = v_2 = 0``. El are și soluții **nenule** doar dacă cele două ecuații spun același lucru (coeficienții lor sînt proporționali), adică dacă **determinantul** sistemului este zero:
```math
(a - \lambda)(d - \lambda) - bc = 0
\quad\Longleftrightarrow\quad
\lambda^2 - (a + d)\,\lambda + (ad - bc) = 0 .
```
Aceasta este **ecuația caracteristică**, pe scurt ``\det(L - \lambda I) = 0``, unde ``I`` este matricea unitate. Coeficienții ei sînt **urma** ``\operatorname{tr} L = a + d`` (suma elementelor de pe diagonală) și **determinantul** ``\det L = ad - bc``. O matrice ``2 \times 2`` are deci cel mult două valori proprii; o matrice ``m \times m`` are cel mult ``m``, rădăcinile unui polinom de gradul ``m``.

**Pentru urs:** ``\operatorname{tr} L = 1.55`` și ``\det L = 0.65 \cdot 0.9 - 0.5 \cdot 0.25 = 0.46``, deci
```math
\lambda^2 - 1.55\,\lambda + 0.46 = 0, \qquad
\lambda_{1,2} = \frac{1.55 \pm \sqrt{1.55^2 - 4 \cdot 0.46}}{2} = \frac{1.55 \pm 0.75}{2}, \qquad
\lambda_1 = 1.15, \quad \lambda_2 = 0.4 .
```
**Vectorii proprii** se obțin din prima ecuație a sistemului (a doua spune același lucru):
* ``\lambda_1 = 1.15``: ``-0.5\,v_1 + 0.5\,v_2 = 0``, deci ``v_1 = v_2`` și ``\mathbf{v}_1 = (1, 1)``;
* ``\lambda_2 = 0.4``: ``0.25\,v_1 + 0.5\,v_2 = 0``, deci ``v_1 = -2v_2`` și ``\mathbf{v}_2 = (2, -1)``.

Orice multiplu nenul al unui vector propriu este tot un vector propriu. Pentru populații îl alegem de obicei cu suma componentelor ``1``: ``(0.5,\ 0.5)``, adică jumătate juvenili și jumătate adulți.

> ⚠️ Ecuația caracteristică poate avea și o **rădăcină dublă** (atunci matricea poate avea o singură direcție proprie) sau **rădăcini complexe** (atunci nicio direcție nu rămîne pe loc: matricea rotește toți vectorii). Aceste cazuri le vom studia în Cursul 8.
"""

# ╔═╡ a866b79f-d1e5-4d0b-8de1-50e1fb2b951d
md"""
### 3.3. Soluția modelului
Vectorii proprii ``\mathbf{v}_1 = (1, 1)`` și ``\mathbf{v}_2 = (2, -1)`` formează o **bază** a lui ``\mathbb{R}^2``: orice stare inițială se scrie în mod unic ca ``\mathbf{x}_0 = c_1\mathbf{v}_1 + c_2\mathbf{v}_2``. Matricea acționează pe fiecare termen separat (liniaritate!), înmulțindu-l cu valoarea proprie corespunzătoare:
```math
\mathbf{x}_1 = L\,\mathbf{x}_0 = c_1\lambda_1\mathbf{v}_1 + c_2\lambda_2\mathbf{v}_2, \qquad \dots, \qquad
\mathbf{x}_n = L^n\mathbf{x}_0 = c_1\,\lambda_1^n\,\mathbf{v}_1 + c_2\,\lambda_2^n\,\mathbf{v}_2 .
```
Problema în două dimensiuni s-a despărțit în două probleme de **creștere geometrică** din Cursul 4, cîte una pe fiecare direcție proprie.

**Exemplu:** pentru ``\mathbf{x}_0 = (100, 10)`` rezolvăm sistemul ``c_1 + 2c_2 = 100``, ``c_1 - c_2 = 10`` și găsim ``c_2 = 30``, ``c_1 = 40``. Deci
```math
J_n = 40 \cdot 1.15^n + 60 \cdot 0.4^n, \qquad A_n = 40 \cdot 1.15^n - 30 \cdot 0.4^n .
```
Pentru ``n = 1`` regăsim ``(70, 34)``, ca în secțiunea 1. Raportul dintre cei doi termeni se înmulțește în fiecare an cu ``0.4 / 1.15 \approx 0.35``, deci după ``5`` ani termenul cu ``0.4^n`` este, relativ, de aproximativ ``200`` de ori mai mic decît la început. Aceasta este durata **regimului tranzitoriu**; după el rămîne practic doar ``c_1 \cdot 1.15^n\,(1, 1)``.

### 🎛️ Experimentați
În figură, starea ``\mathbf{x}_n`` (punctul portocaliu) este suma dintre componenta pe direcția ``\mathbf{v}_1`` (săgeata mov, ``c_1\lambda_1^n\mathbf{v}_1``) și componenta pe direcția ``\mathbf{v}_2`` (săgeata verde, ``c_2\lambda_2^n\mathbf{v}_2``). Schimbați starea inițială și urmăriți cum dispare componenta verde.

Juvenili inițiali ``J_0``: $(@bind J0_s PlutoUI.Slider(0:10:200, show_value=true, default=100))

Adulți inițiali ``A_0``: $(@bind A0_s PlutoUI.Slider(0:10:200, show_value=true, default=10))

Anul ``n``: $(@bind n_an PlutoUI.Slider(0:1:12, show_value=true, default=0))
"""

# ╔═╡ c990f7b4-3b01-4189-93ac-78af3152bb4e
let
	L = [0.65 0.5; 0.25 0.9]
	v₁, v₂ = [1.0, 1.0], [2.0, -1.0]
	c₁, c₂ = (J0_s + 2A0_s) / 3, (J0_s - A0_s) / 3     # x₀ = c₁ v₁ + c₂ v₂
	stari = iteratii(L, [J0_s, A0_s], 12)
	x = stari[n_an+1]
	p₁ = c₁ * 1.15^n_an * v₁     # componenta pe direcția v₁
	p₂ = c₂ * 0.4^n_an * v₂      # componenta pe direcția v₂
	R = max(1.0, 1.1 * maximum(maximum, [stari[1:n_an+1]; [p₁]]))

	fig = Figure(size = (1200, 520))
	ax1 = Axis(fig[1, 1], title = "Starea în anul n = $n_an",
		subtitle = @sprintf("xₙ = %.1f · 1.15ⁿ · (1, 1) + %.1f · 0.4ⁿ · (2, −1)", c₁, c₂),
		xlabel = "Juvenili (J)", ylabel = "Adulți (A)", aspect = DataAspect(),
		limits = (-0.05R, R, -0.05R, R))
	lines!(ax1, [0, R], [0, R], color = :purple, linestyle = :dash, label = "direcția v₁ = (1, 1)")
	scatterlines!(ax1, first.(stari[1:n_an+1]), last.(stari[1:n_an+1]),
		color = (:darkorange, 0.5), markersize = 8)
	arrows2d!(ax1, [0.0], [0.0], [p₁[1]], [p₁[2]], color = :purple,
		shaftwidth = 3, tipwidth = 14, tiplength = 12)
	if norm(p₂) > 1e-3 * R
		arrows2d!(ax1, [p₁[1]], [p₁[2]], [p₂[1]], [p₂[2]], color = :seagreen,
			shaftwidth = 3, tipwidth = 14, tiplength = 12)
	end
	scatter!(ax1, [x[1]], [x[2]], color = :darkorange, markersize = 16, strokewidth = 1.5)
	axislegend(ax1, position = :lt)

	proportii = [s[1] / sum(s) for s in stari]
	ax2 = Axis(fig[1, 2], title = "Proporția juvenililor Jₙ / (Jₙ + Aₙ)",
		xlabel = "Anul (n)", limits = (-0.3, 12.3, 0, 1))
	hlines!(ax2, [0.5], color = :purple, linestyle = :dash)
	scatterlines!(ax2, 0:12, proportii, color = :darkorange, markersize = 8)
	scatter!(ax2, [n_an], [proportii[n_an+1]], color = :darkorange, markersize = 18, strokewidth = 1.5)
	fig
end

# ╔═╡ ac76b83d-2190-4a7d-b705-f92d51798d0d
md"""
**Ce observăm:** pentru orice populație reală ``c_1 = (J_0 + 2A_0)/3 > 0``, iar săgeata verde se micșorează rapid: starea se apropie de dreapta ``\mathbf{v}_1``, iar proporția juvenililor tinde la ``1/2``. Dacă pornim chiar de pe direcția ``\mathbf{v}_1`` (``J_0 = A_0``), atunci ``c_2 = 0`` și nu există regim tranzitoriu.

### 3.4. Valoarea proprie dominantă
Același raționament funcționează pentru o matrice cu ``m`` stadii care are ``m`` vectori proprii independenți:
```math
\mathbf{x}_n = c_1\lambda_1^n\mathbf{v}_1 + c_2\lambda_2^n\mathbf{v}_2 + \dots + c_m\lambda_m^n\mathbf{v}_m ,
\qquad \lvert\lambda_1\rvert > \lvert\lambda_2\rvert \ge \dots \ge \lvert\lambda_m\rvert .
```
Pe termen lung domină termenul cu valoarea proprie de **modul maxim**, ``\lambda_1``, numită **valoarea proprie dominantă**:
* ``\lambda_1`` este **rata de creștere asimptotică** a populației: creștere dacă ``\lambda_1 > 1``, declin dacă ``\lambda_1 < 1``;
* vectorul propriu ``\mathbf{v}_1`` dă **distribuția stabilă pe stadii**;
* raportul ``\lvert\lambda_2\rvert / \lambda_1`` arată cît de repede se stinge regimul tranzitoriu; dacă ``\lambda_2 < 0``, termenul ``\lambda_2^n`` își schimbă semnul în fiecare an, iar regimul tranzitoriu **oscilează**.

### 🎛️ Experimentați
Comparăm trei matrice, pornind de la ``10`` juvenili și ``50`` de adulți.

Scenariul: $(@bind scenariu PlutoUI.Select(["bun" => "An bun: L = (0.65  0.5; 0.25  0.9)", "prost" => "An prost: M = (0.5  0.4; 0.1  0.8)", "scurt" => "Adulți cu viață scurtă: S = (0.1  1.4; 0.4  0.2)"]))
"""

# ╔═╡ 0c6b2c43-3e55-4acc-8b17-c200d070bc0d
let
	matrici = Dict("bun" => [0.65 0.5; 0.25 0.9], "prost" => [0.5 0.4; 0.1 0.8],
		"scurt" => [0.1 1.4; 0.4 0.2])
	P = matrici[scenariu]
	λ = sort(eigvals(P), by = abs, rev = true)
	λ₁, v₁ = dominanta(P)
	ani = 25
	stari = iteratii(P, [10.0, 50.0], ani)
	J, A = first.(stari), last.(stari)

	fig = Figure(size = (1200, 430))
	ax1 = Axis(fig[1, 1], title = "Juvenili și adulți",
		subtitle = @sprintf("valorile proprii: λ₁ = %.2f, λ₂ = %.2f", λ[1], λ[2]),
		xlabel = "Anul (n)", ylabel = "Număr de urși")
	scatterlines!(ax1, 0:ani, J, color = :royalblue, linewidth = 2, label = "juvenili Jₙ")
	scatterlines!(ax1, 0:ani, A, color = :crimson, linewidth = 2, marker = :rect, label = "adulți Aₙ")
	axislegend(ax1, position = scenariu == "bun" ? :lt : :rt)
	ax2 = Axis(fig[1, 2], title = "Proporția juvenililor Jₙ / (Jₙ + Aₙ)",
		subtitle = @sprintf("distribuția stabilă: %.1f%% juvenili", 100 * v₁[1]),
		xlabel = "Anul (n)", limits = (0, ani, 0, 1))
	hlines!(ax2, [v₁[1]], color = :purple, linestyle = :dash)
	scatterlines!(ax2, 0:ani, J ./ (J .+ A), color = :purple, linewidth = 2)
	fig
end

# ╔═╡ abadcb6b-3501-4f24-9d75-c00506be86fe
md"""
**Ce observăm:**
* **an bun** (``\lambda_1 = 1.15``, ``\lambda_2 = 0.4``): creștere cu ``15\%`` pe an, cu jumătate juvenili;
* **an prost** (``\lambda_1 = 0.9``, ``\lambda_2 = 0.4``): declin cu ``10\%`` pe an. În primii doi ani numărul juvenililor **crește** (de la ``10`` la aproape ``29``): un efect tranzitoriu al stării inițiale cu prea mulți adulți (Verificarea rapidă de mai jos);
* **adulți cu viață scurtă** (``\lambda_1 = 0.9``, ``\lambda_2 = -0.6``): tot declin, dar în primii ani populația alternează între ani cu mulți juvenili și ani cu mulți adulți, pentru că ``(-0.6)^n`` își schimbă semnul.

> 🧭 **Teorema Perron–Frobenius.** De ce valoarea proprie dominantă a unui model de populație este un număr **real și pozitiv**, iar vectorul ei propriu are toate componentele **pozitive**, adică este o adevărată structură de populație? Dacă matricea are toate elementele ``\ge 0``, graful ciclului de viață permite să se ajungă din orice stadiu în orice alt stadiu, iar cel mai mare divizor comun al lungimilor ciclurilor grafului este ``1`` (de exemplu, dacă graful are o buclă), atunci există o valoare proprie reală ``\lambda_1 > 0``, mai mare în modul decît toate celelalte, cu un vector propriu cu toate componentele pozitive.

### 3.5. Calculul în Julia
Pentru matrice mai mari, valorile și vectorii proprii se calculează numeric, cu funcțiile `eigvals` și `eigvecs` din pachetul `LinearAlgebra`. Coeficienții ``c_1, c_2`` se obțin rezolvînd un sistem liniar: în Julia, `[1 2; 1 -1] \ [100, 10]`.
"""

# ╔═╡ dde6ae4f-e276-427a-9e00-817287f73923
let
	L = [0.65 0.5; 0.25 0.9]
	λ = eigvals(L)                 # valorile proprii, în ordine crescătoare
	V = eigvecs(L)                 # vectorii proprii, pe coloane (de lungime 1)
	v₁ = V[:, 2] / sum(V[:, 2])    # vectorul propriu al lui λ = 1.15, cu suma componentelor 1
	λ, v₁
end

# ╔═╡ 9eee84cd-35ac-489d-b00c-62a9185c623e
md"""
### ✅ Verificare rapidă
Pentru matricea anului prost ``M = \begin{pmatrix} 0.5 & 0.4 \\ 0.1 & 0.8 \end{pmatrix}``:

a) calculați valorile și vectorii proprii;

b) scrieți starea inițială ``(10, 50)`` ca o combinație a vectorilor proprii și explicați de ce numărul juvenililor crește la început.
"""

# ╔═╡ 65700614-5396-4e53-9841-fe64d2876b9c
hint(md"a) ``\operatorname{tr} M = 1.3``, ``\det M = 0.4 - 0.04 = 0.36``, deci ``\lambda^2 - 1.3\,\lambda + 0.36 = 0`` și ``\lambda_1 = 0.9``, ``\lambda_2 = 0.4``, cu vectorii proprii ``(1, 1)`` și ``(4, -1)``. b) ``(10, 50) = 42\,(1, 1) - 8\,(4, -1)``, deci ``J_n = 42 \cdot 0.9^n - 32 \cdot 0.4^n`` și ``A_n = 42 \cdot 0.9^n + 8 \cdot 0.4^n``. Termenul negativ ``-32 \cdot 0.4^n`` micșorează la început numărul juvenililor, dar se stinge mai repede decît ``42 \cdot 0.9^n``; de aceea ``J_n`` crește întîi (``J_1 = 25``, ``J_2 \approx 28.9``) și abia apoi scade. Biologic: o populație cu prea puțini juvenili față de structura stabilă își completează întîi clasa juvenililor.")

# ╔═╡ 1e8b77b3-525d-4b50-ab0f-31b4a875a204
md"""
## 🐢 4. Aplicație: cum salvăm țestoasa marină *Caretta caretta*? (15 min)
**Țestoasa marină cu cap mare** (*Caretta caretta*) trăiește peste ``50`` de ani și se maturizează abia după ``20`` de ani. Femelele revin pe plajele unde s-au născut ca să își depună ouăle. În anii 1980 populația de pe coasta de sud-est a SUA era în declin, iar eforturile de conservare se concentrau pe **protejarea cuiburilor** de pe plaje. D. Crouse, L. Crowder și H. Caswell (1987) au construit un model matricial cu ``7`` stadii pentru a afla unde ar avea efect cea mai mare intervenție:

| Stadiul | | Vîrsta (ani) |
| :--- | :--- | :--- |
| 1 | ouă și pui | sub ``1`` |
| 2 | juvenili mici | ``1``–``7`` |
| 3 | juvenili mari | ``8``–``15`` |
| 4 | subadulți | ``16``–``21`` |
| 5 | adulți la prima reproducere | ``22`` |
| 6 | adulți în al doilea an de reproducere | ``23`` |
| 7 | adulți maturi | ``24``–``54`` |

```math
T = \begin{pmatrix}
0 & 0 & 0 & 0 & 127 & 4 & 80 \\
0.6747 & 0.7370 & 0 & 0 & 0 & 0 & 0 \\
0 & 0.0486 & 0.6610 & 0 & 0 & 0 & 0 \\
0 & 0 & 0.0147 & 0.6907 & 0 & 0 & 0 \\
0 & 0 & 0 & 0.0518 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0.8091 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0.8091 & 0.8089
\end{pmatrix}
```
Primul rînd conține numărul mediu de ouă (din care vor ieși femele) depuse de o femelă într-un an; valorile diferă de la un stadiu la altul pentru că o femelă nu cuibărește în fiecare an. Graful ciclului de viață:
"""

# ╔═╡ 8c030d41-811e-4b2c-b1f8-e52ef8e9a05a
T_caretta = [0 0 0 0 127 4 80;
	0.6747 0.7370 0 0 0 0 0;
	0 0.0486 0.6610 0 0 0 0;
	0 0 0.0147 0.6907 0 0 0;
	0 0 0 0.0518 0 0 0;
	0 0 0 0 0.8091 0 0;
	0 0 0 0 0 0.8091 0.8089]

# ╔═╡ 2b18c49b-0f9c-4af3-a341-9ac488233baf
graf_ciclu(T_caretta, string.(1:7); id = "caretta")

# ╔═╡ cbadc430-762e-463c-9e6c-4503d1d7e932
dominanta(T_caretta)   # (λ₁, distribuția stabilă pe stadii)

# ╔═╡ 946cc39d-b4d7-48a3-9c45-dda0c304013a
md"""
**Ce spune modelul:** ``\lambda_1 \approx 0.945``. Fără intervenție, populația scade cu aproximativ ``5.5\%`` pe an, deci se înjumătățește în ``\ln 2 / \ln(1/0.945) \approx 12`` ani. În distribuția stabilă, ``88\%`` dintre țestoase sînt ouă, pui și juvenili mici, iar adulții (stadiile ``5``–``7``) sînt sub ``0.3\%``.

### 🎛️ Experimentați
Să presupunem că o măsură de protecție reduce **mortalitatea** anuală a unui stadiu cu un anumit procent (de exemplu, o scădere a mortalității de la ``30\%`` la ``15\%`` este o reducere cu ``50\%``), păstrînd aceeași proporție între indivizii care rămîn în stadiu și cei care trec mai departe. Cum se schimbă ``\lambda_1``? Ultima bară corespunde protejării simultane a stadiilor ``3``–``7``.

Reducerea mortalității (%): $(@bind reducere PlutoUI.Slider(0:5:100, show_value=true, default=50))
"""

# ╔═╡ e70aec2b-a471-4f9b-a030-57b3ada9efc0
let
	p = reducere / 100
	# reduce cu fracțiunea p mortalitatea stadiilor din lista stadii
	function protejeaza(T, stadii, p)
		T = copy(T)
		for k in stadii
			σ = T[k, k] + (k < 7 ? T[k+1, k] : 0.0)    # supraviețuirea anuală în stadiul k
			c = (σ + p * (1 - σ)) / σ
			T[k, k] *= c
			k < 7 && (T[k+1, k] *= c)
		end
		T
	end
	grupuri = [[[k] for k in 1:7]; [3:7]]
	λs = [dominanta(protejeaza(T_caretta, g, p))[1] for g in grupuri]
	λ₀, v₀ = dominanta(T_caretta)

	fig = Figure(size = (1250, 480))
	ax1 = Axis(fig[1, 1], title = "Rata de creștere λ₁ după protejarea unui stadiu",
		subtitle = "mortalitatea redusă cu $(reducere)%; linia gri: fără protecție (λ₁ ≈ 0.945)",
		xlabel = "Stadiul protejat", ylabel = "λ₁",
		xticks = (1:8, ["1\nouă, pui", "2\njuv. mici", "3\njuv. mari", "4\nsubadulți",
			"5\nprima repr.", "6\nal doilea an", "7\nadulți", "3–7"]),
		limits = (0.4, 8.6, 0.93, max(1.08, maximum(λs) + 0.03)))
	barplot!(ax1, 1:8, λs, fillto = 0.93, color = [λ > 1 ? :seagreen : :indianred for λ in λs],
		bar_labels = :y, label_formatter = x -> @sprintf("%.3f", x))
	hlines!(ax1, [1.0], color = :black, linestyle = :dash)
	hlines!(ax1, [λ₀], color = :gray50)

	ax2 = Axis(fig[1, 2], title = "Distribuția stabilă pe stadii (fără protecție)",
		xlabel = "Stadiul", ylabel = "% din populație", xticks = 1:7, limits = (0.4, 7.6, 0, 80))
	barplot!(ax2, 1:7, 100 .* v₀, color = :steelblue,
		bar_labels = :y, label_formatter = x -> @sprintf("%.2f%%", x))
	colsize!(fig.layout, 1, Relative(0.62))
	fig
end

# ╔═╡ 9eb95c7e-7816-4136-944d-53bb438b4467
md"""
### 🤔 Ce am observat?
* Chiar dacă **toate** ouăle și toți puii ar supraviețui (reducere cu ``100\%`` în stadiul ``1``), ``\lambda_1 \approx 0.965 < 1``: protejarea cuiburilor, singură, **nu poate opri** declinul.
* Cel mai mare efect îl are protejarea **juvenililor mari** (stadiul ``3``): o reducere a mortalității lor cu aproximativ ``57\%`` ajunge pentru ``\lambda_1 = 1``. Dacă reducem mortalitatea în toate stadiile ``3``–``7`` (juvenili mari, subadulți și adulți), ajunge o reducere de aproximativ ``26\%``.
* Stadiile cele mai **numeroase** (ouăle și juvenilii mici) nu sînt cele mai **importante** pentru creșterea populației. Un juvenil mare a supraviețuit deja anilor periculoși și are în față mulți ani de reproducere: pierderea lui costă mult mai mult decît pierderea unui ou.

Concluzia a avut consecințe practice. Juvenilii mari și adulții mor mai ales în **plasele de pescuit de creveți**; pe baza acestor rezultate, la sfîrșitul anilor 1980 și începutul anilor 1990 au devenit obligatorii în SUA dispozitivele **TED** (*Turtle Excluder Devices*), grilaje prin care țestoasele scapă din plase.

### ✅ Verificare rapidă
De ce reducerea mortalității în stadiile ``5`` și ``6`` aproape nu schimbă ``\lambda_1``, deși sînt stadii de adulți, în timp ce în stadiul ``7`` efectul este mare? (Indicație: cît timp petrece o țestoasă în fiecare dintre aceste stadii?)
"""

# ╔═╡ 7e47a6d6-e9d1-4650-97ea-f766bdb26025
hint(md"Stadiile ``5`` și ``6`` durează cîte un singur an, iar supraviețuirea în ele este deja mare (``0.81``): chiar și o supraviețuire de ``100\%`` adaugă puțin (``\lambda_1 \approx 0.953``). Stadiul ``7`` durează pînă la ``30`` de ani: o creștere a supraviețuirii anuale se acumulează pe toată durata vieții reproductive, deci multiplică numărul de ani în care o femelă depune ouă.")

# ╔═╡ 8f2c515a-37a9-4e4a-982f-9b57fa44d7d9
md"""
## 📝 5. Recapitulare și exerciții (10 min)

### 🧠 Recapitulare
- 🧭 Starea unei populații cu ``m`` stadii este un vector din ``\mathbb{R}^m``. O aplicație **liniară** este determinată de imaginile vectorilor bazei, așezate pe coloanele **matricei** ei; coloana ``j`` spune ce devine într-un an un individ din stadiul ``j``.
- 🔄 Modelele matriciale ``\mathbf{x}_{n+1} = L\,\mathbf{x}_n`` se construiesc din **graful ciclului de viață** (Leslie: clase de vîrstă; Lefkovitch: stadii). Succesiunea a doi ani corespunde **produsului** matricelor, care **nu este comutativ**.
- 🎯 **Valori și vectori proprii:** ``L\mathbf{v} = \lambda\mathbf{v}``. Valorile proprii sînt rădăcinile **ecuației caracteristice** ``\det(L - \lambda I) = 0``; pentru o matrice ``2 \times 2``: ``\lambda^2 - (\operatorname{tr} L)\,\lambda + \det L = 0``.
- 📈 ``\mathbf{x}_n = c_1\lambda_1^n\mathbf{v}_1 + c_2\lambda_2^n\mathbf{v}_2 + \dots``: pe termen lung domină **valoarea proprie dominantă** ``\lambda_1`` (rata de creștere), iar structura populației tinde la vectorul propriu ``\mathbf{v}_1`` (**distribuția stabilă pe stadii**). Ceilalți termeni descriu **regimul tranzitoriu**.
- 🐢 Modelele matriciale arată **unde** trebuie să acționăm: la țestoasa *Caretta caretta*, protejarea juvenililor mari și a adulților contează mult mai mult decît protejarea ouălor.

| Valorile proprii | Comportament | Exemple din curs |
| :--- | :--- | :--- |
| ``\lambda_1 > 1`` | creștere geometrică, cu o structură stabilă pe stadii | ursul, an bun (``\lambda_1 = 1.15``) |
| ``\lambda_1 < 1`` | declin geometric | ursul, an prost (``\lambda_1 = 0.9``); țestoasa (``\lambda_1 \approx 0.945``) |
| ``\lambda_2 < 0`` | regim tranzitoriu care oscilează | adulți cu viață scurtă (``\lambda_2 = -0.6``) |

---
### ✍️ Exerciții
**Exercițiul 1 (iepurii lui Fibonacci).** În cartea *Liber Abaci* (1202), Leonardo din Pisa (Fibonacci) urmărește o populație de iepuri: o pereche de iepuri tineri devine adultă după o lună, iar fiecare pereche adultă produce în fiecare lună o pereche de iepuri tineri; iepurii nu mor. Cu ``J_n`` perechi tinere și ``A_n`` perechi adulte în luna ``n``:
```math
J_{n+1} = A_n, \qquad A_{n+1} = J_n + A_n .
```
a) Scrieți matricea modelului și desenați graful ciclului de viață.

b) Calculați valorile proprii. Arătați că valoarea proprie dominantă este **numărul de aur** ``\varphi = \dfrac{1 + \sqrt{5}}{2} \approx 1.618``.

c) Pornind de la o pereche tînără, numărul total de perechi este ``1, 1, 2, 3, 5, 8, 13, \dots`` (șirul lui Fibonacci). Explicați de ce raportul a doi termeni consecutivi tinde la ``\varphi`` și găsiți distribuția stabilă.
"""

# ╔═╡ 4152ebbe-2912-48dd-a1e1-3bdecea37dd1
hint(md"a) ``\begin{pmatrix} 0 & 1 \\ 1 & 1 \end{pmatrix}``: o săgeată ``J \to A`` (``1``), o buclă pe ``A`` (``1``) și o săgeată de reproducere ``A \to J`` (``1``). b) ``\operatorname{tr} = 1``, ``\det = -1``, deci ``\lambda^2 - \lambda - 1 = 0`` și ``\lambda_{1,2} = \frac{1 \pm \sqrt{5}}{2}``: ``\lambda_1 = \varphi \approx 1.618``, ``\lambda_2 \approx -0.618``. c) Cum ``\lvert\lambda_2\rvert < \lambda_1``, termenul ``c_2\lambda_2^n\mathbf{v}_2`` se stinge (alternînd semnul), iar populația se înmulțește pe termen lung cu ``\varphi`` în fiecare lună. Vectorul propriu: ``-\varphi\,v_1 + v_2 = 0``, deci ``\mathbf{v}_1 = (1, \varphi)``; normalizat, ``(0.382,\ 0.618)``: ``38.2\%`` perechi tinere și ``61.8\%`` adulte.")

# ╔═╡ faf388a7-a8bb-4bd0-b938-4c43adc9763f
md"""
**Exercițiul 2 (pragul de supraviețuire, Cursul 5).** Pentru modelul pe stadii ``L = \begin{pmatrix} a & f \\ m & s \end{pmatrix}``:

a) Arătați că ``\lambda = 1`` este valoare proprie dacă și numai dacă ``(1 - a)(1 - s) = f\,m``. (Indicație: înlocuiți ``\lambda = 1`` în ecuația caracteristică.)

b) Pentru urs (``a = 0.65``, ``m = 0.25``, ``f = 0.5``), găsiți supraviețuirea adulților ``s`` pentru care ``\lambda_1 = 1``. Comparați cu Verificarea rapidă din secțiunea 5 a Cursului 5.

c) Braconajul reduce ``s`` de la ``0.9`` la ``0.6``. Calculați noua valoare proprie dominantă. Ce se întîmplă cu populația?
"""

# ╔═╡ 1502e673-a8cd-43c7-80f2-9c22556b87dd
hint(md"a) Ecuația caracteristică este ``\lambda^2 - (a + s)\lambda + (as - fm) = 0``; pentru ``\lambda = 1``: ``1 - a - s + as - fm = (1 - a)(1 - s) - fm = 0``. b) ``0.35\,(1 - s) = 0.125``, deci ``s = 1 - 0.125/0.35 \approx 0.643``, același prag ca în Cursul 5. c) ``\operatorname{tr} L = 1.25``, ``\det L = 0.39 - 0.125 = 0.265``, deci ``\lambda^2 - 1.25\,\lambda + 0.265 = 0`` și ``\lambda_1 = \frac{1.25 + \sqrt{0.5025}}{2} \approx 0.979``: populația scade cu aproximativ ``2\%`` pe an.")

# ╔═╡ 41290751-b1f0-4bbc-b10c-28b7e7192375
md"""
**Exercițiul 3 (ordinea anilor).**

a) Arătați că pentru orice două matrice ``2 \times 2`` avem ``\operatorname{tr}(ML) = \operatorname{tr}(LM)``. Folosind proprietatea ``\det(ML) = \det M \cdot \det L``, deduceți că ``ML`` și ``LM`` au aceleași valori proprii.

b) Calculați valorile proprii ale lui ``ML`` pentru matricele din secțiunea 2.2. Cu ce factor se înmulțește, în medie, populația într-un an, dacă anii buni și proști alternează? Crește sau scade populația?

c) Verificați că ``\lambda_1(ML) = 1.15 \cdot 0.9``. De ce se întîmplă acest lucru aici? Este adevărat pentru orice două matrice?
"""

# ╔═╡ aa7f4738-7827-4062-af06-8fcb13115cbd
hint(md"a) ``\operatorname{tr}(ML) = m_{11}l_{11} + m_{12}l_{21} + m_{21}l_{12} + m_{22}l_{22}``, o expresie care nu se schimbă dacă schimbăm rolurile lui ``M`` și ``L``; ``\det(ML) = \det M \det L = \det(LM)``. Ecuația caracteristică depinde doar de urmă și de determinant, deci este aceeași. b) ``\operatorname{tr}(ML) = 1.195``, ``\det(ML) = 0.36 \cdot 0.46 = 0.1656``, deci ``\lambda_1 = 1.035`` și ``\lambda_2 = 0.16``. În doi ani populația se înmulțește cu ``1.035``, deci în medie cu ``\sqrt{1.035} \approx 1.017`` pe an: o creștere lentă, de aproximativ ``1.7\%`` pe an. c) ``L`` și ``M`` au **același** vector propriu dominant ``(1, 1)``, deci ``ML\,(1, 1) = M\,(1.15\,(1, 1)) = 1.15 \cdot 0.9\,(1, 1)``. În general nu este adevărat: dacă vectorii proprii diferă, valoarea proprie dominantă a produsului nu este produsul valorilor proprii dominante.")

# ╔═╡ 10c32bb8-0dc9-4b50-ba9d-3a408b65f575
md"""
**Exercițiul 4 (un model Leslie).** O specie trăiește cel mult trei ani. Femelele de un an nu se reproduc, cele de doi ani au în medie ``2`` fiice, iar cele de trei ani ``4`` fiice. Jumătate dintre femelele de un an ajung la doi ani, iar un sfert dintre cele de doi ani ajung la trei ani.

a) Scrieți matricea Leslie și desenați graful ciclului de viață.

b) **Numărul net de reproducere** ``R_0 = f_1 + s_1 f_2 + s_1 s_2 f_3`` este numărul mediu de fiice pe care o femelă le are în toată viața ei. Calculați ``R_0``.

c) Rulați codul de mai jos și găsiți ``\lambda_1`` și distribuția stabilă pe vîrste.

d) Arătați că ecuația caracteristică se poate scrie ``1 = \dfrac{f_1}{\lambda} + \dfrac{s_1 f_2}{\lambda^2} + \dfrac{s_1 s_2 f_3}{\lambda^3}`` (**ecuația Euler–Lotka**) și deduceți că ``\lambda_1 > 1`` dacă și numai dacă ``R_0 > 1``.
"""

# ╔═╡ d97096cb-d48e-4486-9352-7ea595f0a37b
# EXERCIȚIUL 4 c): modificați matricea și rulați
let
	L = [0.0 2.0 4.0;
		0.5 0.0 0.0;
		0.0 0.25 0.0]
	λ₁, v₁ = dominanta(L)
	(λ₁ = λ₁, distributie_stabila = v₁)
end

# ╔═╡ 75069a38-129e-4993-94d0-93edc1090ce4
hint(md"a) ``L = \begin{pmatrix} 0 & 2 & 4 \\ 0.5 & 0 & 0 \\ 0 & 0.25 & 0 \end{pmatrix}``. b) ``R_0 = 0 + 0.5 \cdot 2 + 0.5 \cdot 0.25 \cdot 4 = 1.5``. c) ``\lambda_1 \approx 1.191``, cu distribuția stabilă ``(66\%,\ 28\%,\ 6\%)``. d) ``\det(L - \lambda I) = -\lambda^3 + f_1\lambda^2 + s_1 f_2\,\lambda + s_1 s_2 f_3``; împărțim ecuația ``\lambda^3 = f_1\lambda^2 + s_1 f_2\,\lambda + s_1 s_2 f_3`` prin ``\lambda^3``. Membrul drept este descrescător în ``\lambda > 0`` și este egal cu ``R_0`` pentru ``\lambda = 1``; deci soluția pozitivă ``\lambda_1`` este mai mare decît ``1`` exact cînd ``R_0 > 1``.")

# ╔═╡ 7c742114-39a2-4650-aede-bab345d9a859
md"""
**Exercițiul 5 (țestoasa: mai multe ouă?).** O altă idee de conservare ar fi creșterea **fecundității** (de exemplu, prin incubarea artificială a ouălor), adică înmulțirea primului rînd al matricei ``T`` cu un factor ``k``. Modificați ``k`` în codul de mai jos: cît de mare trebuie să fie ``k`` pentru ca populația să nu mai scadă? Este realist?
"""

# ╔═╡ 686914f1-06b9-440d-a8b5-77de679fd767
# EXERCIȚIUL 5: modificați factorul k și rulați
let
	k = 1.5
	T = copy(T_caretta)
	T[1, :] .*= k
	dominanta(T)[1]
end

# ╔═╡ c22194b2-90c7-4826-bee7-62d6db5197c8
hint(md"Pentru ``k = 1.5`` obținem ``\lambda_1 \approx 0.966``, iar ``\lambda_1 = 1`` abia pentru ``k \approx 2.73``: fiecare femelă ar trebui să producă de aproape trei ori mai mulți pui care ajung în mare, ceea ce nu este realist. Concluzia lui Crouse și a colaboratorilor ei rămîne valabilă: cheia este supraviețuirea juvenililor mari și a adulților.")

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─ec01674b-79a9-4a71-9a12-2e04384599da
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╟─186cde3c-94c6-46cf-86b7-bc61461674d2
# ╟─644a8297-b028-4515-b29f-0600e67ac8f9
# ╟─6ece3812-853a-4c42-a863-d9430a9bdcb3
# ╟─647b474a-ec17-465b-8edb-c9f7ae631949
# ╠═b187da58-91ab-4f2c-8955-6634c7c47dbe
# ╟─ad470386-c9ca-4461-81e1-7fd11639b7b1
# ╟─55b17dcd-ccfb-4c30-9cf2-d1d0341d44e5
# ╟─c025dcf9-37bd-43c9-a77f-540d33389e08
# ╟─7e9c12a7-ffc9-4874-8682-b713e14ac3dc
# ╟─8a49308f-1ff7-4ad1-87c4-e25f8bb4d541
# ╟─e1a18af2-f6ee-4945-a071-ae2138e86b49
# ╟─e3590799-a03c-441c-8493-254909f373c3
# ╠═80504029-1e38-4b0f-a791-262e1991260e
# ╟─b76aecfd-8aa2-4f8f-8458-e04601a28830
# ╟─64bd52b4-2bd1-4633-85d6-a4869c80f87b
# ╟─d3e5263a-3fe4-4e6b-94c2-ab0de2d53b2b
# ╟─eb0af3c2-c015-496a-8e2d-112e8cf50f66
# ╟─d5e51b29-2530-4e8d-a0cb-33ea3cb8586d
# ╟─a866b79f-d1e5-4d0b-8de1-50e1fb2b951d
# ╟─c990f7b4-3b01-4189-93ac-78af3152bb4e
# ╟─ac76b83d-2190-4a7d-b705-f92d51798d0d
# ╟─0c6b2c43-3e55-4acc-8b17-c200d070bc0d
# ╟─abadcb6b-3501-4f24-9d75-c00506be86fe
# ╠═dde6ae4f-e276-427a-9e00-817287f73923
# ╟─9eee84cd-35ac-489d-b00c-62a9185c623e
# ╟─65700614-5396-4e53-9841-fe64d2876b9c
# ╟─1e8b77b3-525d-4b50-ab0f-31b4a875a204
# ╟─8c030d41-811e-4b2c-b1f8-e52ef8e9a05a
# ╟─2b18c49b-0f9c-4af3-a341-9ac488233baf
# ╠═cbadc430-762e-463c-9e6c-4503d1d7e932
# ╟─946cc39d-b4d7-48a3-9c45-dda0c304013a
# ╟─e70aec2b-a471-4f9b-a030-57b3ada9efc0
# ╟─9eb95c7e-7816-4136-944d-53bb438b4467
# ╟─7e47a6d6-e9d1-4650-97ea-f766bdb26025
# ╟─8f2c515a-37a9-4e4a-982f-9b57fa44d7d9
# ╟─4152ebbe-2912-48dd-a1e1-3bdecea37dd1
# ╟─faf388a7-a8bb-4bd0-b938-4c43adc9763f
# ╟─1502e673-a8cd-43c7-80f2-9c22556b87dd
# ╟─41290751-b1f0-4bbc-b10c-28b7e7192375
# ╟─aa7f4738-7827-4062-af06-8fcb13115cbd
# ╟─10c32bb8-0dc9-4b50-ba9d-3a408b65f575
# ╠═d97096cb-d48e-4486-9352-7ea595f0a37b
# ╟─75069a38-129e-4993-94d0-93edc1090ce4
# ╟─7c742114-39a2-4650-aede-bab345d9a859
# ╠═686914f1-06b9-440d-a8b5-77de679fd767
# ╟─c22194b2-90c7-4826-bee7-62d6db5197c8
