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

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000000
begin
    using Pkg
    Pkg.activate(".")
    using Makie, CairoMakie
    using PlutoUI
    using PlutoTeachingTools
    using LinearAlgebra
    # using LoopVectorization
    TableOfContents()
end

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 🧮 Cursul 7 — Algebră Liniară și Modele Matriciale

## 🎯 Obiectivele lecției
1. Înțelegerea conceptelor de **aplicație liniară**, **bază** și **matrice**.
2. Modelarea dinamicii populațiilor structurate pe vârste folosind **modele Leslie**.
3. Analiza **comportamentului pe termen lung** (tranzient vs. asimptotic).
4. Introducerea **valorilor și vectorilor proprii** ca instrumente pentru determinarea riguroasă a dinamicii sistemelor liniare.
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000002
md"""
## 🗺️ 1. Aplicații Liniare și Sisteme Dinamice

Până acum am observat comportamentul sistemelor dinamice prin simulări. Pentru a studia *teoria calitativă* a acestora și a putea face calcule riguroase, avem nevoie de un instrument teoretic fundamental: **algebra liniară**.

Vom studia aplicații liniare $f:\mathbb{R}^n \to \mathbb{R}^n$ în două contexte paralele:
- **Sisteme discrete:** $\overrightarrow{X_{N+1}} = f\left(\overrightarrow{X_N}\right)$
- **Ecuații diferențiale (continue):** $\overrightarrow{X'} = f\left(\overrightarrow{X}\right)$

Elementele lui $\mathbb{R}^n$ vor fi notate cu $\overrightarrow{x}$ și scrise sub formă de coloană:
$$\overrightarrow{x} = \begin{pmatrix} x_1 \\ \vdots \\ x_n \end{pmatrix}$$

### 🧩 Combinații liniare și baze
În $\mathbb{R}^n$ avem o mulțime specială de vectori, **baza canonică**:
$$\overrightarrow{e_1} = (1,0,\dots,0), \quad \overrightarrow{e_2} = (0,1,\dots,0), \quad \dots, \quad \overrightarrow{e_n} = (0,0,\dots,1)$$
Orice alt vector $\overrightarrow{x}$ se scrie în mod unic ca o **combinație liniară** a acestora: $\overrightarrow{x} = x_1\overrightarrow{e_1} + \dots + x_n\overrightarrow{e_n}$.

O **bază** este orice sistem de vectori liniar independenți $\{b_1,\dots,b_n\}$ astfel încât orice vector din spațiu se poate scrie unic ca o combinație liniară a lor. Coordonatele lui $\overrightarrow{x}$ în raport cu o nouă bază (formată din coloanele unei matrice $\Lambda$) se obțin prin $\Lambda^{-1}\overrightarrow{x}$.

### 📐 Aplicații liniare
O aplicație $f:\mathbb{R}^n \to \mathbb{R}^m$ este **liniară** dacă păstrează adunarea și înmulțirea cu scalari:
$$f(\overrightarrow{x}+\overrightarrow{y}) = f(\overrightarrow{x}) + f(\overrightarrow{y})$$
$$f(a\cdot\overrightarrow{x}) = a\cdot f(\overrightarrow{x})$$

O proprietate remarcabilă a aplicațiilor liniare este că ele sunt **complet determinate** de acțiunea lor asupra vectorilor din baza canonică. De exemplu, pentru $f:\mathbb{R}^2 \to \mathbb{R}^2$, dacă știm unde se duc $\overrightarrow{e_1}$ și $\overrightarrow{e_2}$, știm unde se duce orice vector. Această informație este organizată într-un tabel numit **matricea** aplicației $f$:
$$f(\overrightarrow{x}) = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = A \cdot \overrightarrow{x}$$
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000003
md"""
## 🐻 2. Un Model Matricial de Populație (Ursul American)

Să aplicăm aceste concepte pentru a modela populația ursului american (*Ursus Americanus*). 
"""





# ╔═╡ 6ece3812-853a-4c42-a863-d9430a9bdcb3
PlutoUI.LocalResource("./Black_Bear.jpg", :width => 800)

# ╔═╡ 4865f5f6-92e9-4a74-94ec-9bd2aba29474
md"""
Femelele ating maturitatea sexuală la 3-4 ani și trăiesc 15-20 de ani. Vom împărți populația de femele în două clase de vârstă:
1. **Juvenili ($J$):** pui și animale nemature.
2. **Adulți ($A$):** animale mature.

**Regulile de tranziție (în fiecare an):**
- O femelă adultă naște în medie $0.5$ pui (femele) pe an.
- $10\%$ din juvenili mor, $25\%$ se maturizează (devin adulți), deci $65\%$ rămân juvenili.
- $10\%$ din adulți mor (speranța de viață ca adult e de ~10 ani), deci $90\%$ rămân adulți.

Modelul de evoluție este un sistem liniar discret:
$$\begin{pmatrix} J_{N+1} \\ A_{N+1} \end{pmatrix} = \begin{pmatrix} 0.65 & 0.5 \\ 0.25 & 0.9 \end{pmatrix} \begin{pmatrix} J_N \\ A_N \end{pmatrix} = A \cdot \begin{pmatrix} J_N \\ A_N \end{pmatrix}$$
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000004
let
    A = [0.65 0.5; 0.25 0.9]
    init = [100, 50]
    next_year = A * init
    println("Populația inițială: $(init[1]) juvenili, $(init[2]) adulți.")
    println("Populația în anul următor: $(round(next_year[1], digits=1)) juvenili, $(round(next_year[2], digits=1)) adulți.")
end

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000005
md"""
### 🔄 Compunerea funcțiilor și produsul matricelor

În natură, condițiile de mediu nu sunt constante. Să ne imaginăm că populația de urși trece printr-un "an bun" (guvernat de matricea $A$), urmat de un "an prost" (guvernat de matricea $B$). 

Cum calculăm populația după acești doi ani?
1. Starea după primul an (anul bun): $\overrightarrow{x}_1 = A \overrightarrow{x}_0$
2. Starea după al doilea an (anul prost): $\overrightarrow{x}_2 = B \overrightarrow{x}_1 = B (A \overrightarrow{x}_0)$

Datorită asociativității înmulțirii matricelor, putem grupa termenii:
$$\overrightarrow{x}_2 = (B \cdot A) \overrightarrow{x}_0$$
Matricea compusă $C = B \cdot A$ descrie astfel evoluția pe o perioadă de doi ani.

#### 🧮 Regula de înmulțire (Rând × Coloană)
Spre deosebire de numerele reale, **înmulțirea matricelor NU este comutativă**. Adică, în general, $B \cdot A \neq A \cdot B$. 
Pentru două matrice $2 \times 2$, produsul $C = B \cdot A$ se calculează astfel:

$$C = \begin{pmatrix} b_{11} & b_{12} \\ b_{21} & b_{22} \end{pmatrix} \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} = \begin{pmatrix} b_{11}a_{11} + b_{12}a_{21} & b_{11}a_{12} + b_{12}a_{22} \\ b_{21}a_{11} + b_{22}a_{21} & b_{21}a_{12} + b_{22}a_{22} \end{pmatrix}$$

> **💡 De ce contează ordinea în biologie?** 
> Dacă întâi se nasc puii și apoi o parte din populație moare ($A$ apoi $B$), rezultatul este diferit față de situația în care întâi mor o parte din adulți și abia apoi se nasc puii ($B$ apoi $A$). Acest lucru se întâmplă deoarece **natalitatea depinde de numărul de adulți rămași** în viață la momentul reproducerii.
"""

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000006
let
    A = [0.65 0.5; 0.25 0.9]
    B = [0.5 0.4; 0.1 0.8]
    init = [100, 50]
    C = B * A
    final_pop = C * init
    println("Matricea compusă C = B * A:\n", C)
    println("Populația finală: $(round.(final_pop, digits=1))")
end

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000007
md"""
## 📈 3. Comportarea pe Termen Lung

Modelele de acest tip se numesc **modele Leslie**. Pornind cu o condiție inițială $\overrightarrow{x}_0$, populația după $N$ ani este:
$$\overrightarrow{x}_N = A^N \overrightarrow{x}_0$$

Să observăm comportamentul pe termen lung pentru cei trei scenarii (an bun, an prost, și un al treilea model oscilant).
"""

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000008
md"""
### 🌲 Scenariul 1: Anul Bun (Creștere exponențială)
"""

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000009
let
    A = [0.65 0.5; 0.25 0.9]
    x = [50.0, 10.0]
    v_j, v_a = Float64[], Float64[]
    n_points = 20
    count = 1:n_points
    
    for n in count
        push!(v_j, x[1])
        push!(v_a, x[2])
        x = A * x
    end
    
    fig1 = Figure(size=(800, 450))
    ax = Axis(fig1[1,1], title="Populația de urși - An Bun", 
              xlabel="Anul (N)", ylabel="Număr de urși")
    colors = Makie.wong_colors()
    
    scatterlines!(ax, count, v_j, label="Juvenili", color=colors[1], marker=:circle, linewidth=2)
    scatterlines!(ax, count, v_a, label="Adulți", color=colors[2], marker=:rect, linewidth=2)
    axislegend(ax, position=:lt)
    fig1
end

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000010
md"""
**Interpretare:** După câteva iterații, ambele clase de vârstă cresc exponențial, iar **raportul** dintre juvenili și adulți devine constant.
"""

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000011
md"""
### 🥀 Scenariul 2: Anul Prost (Declin și Tranzient)
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000012
let
    A = [0.5 0.4; 0.1 0.8]
    x = [10.0, 50.0]
    v_j, v_a = Float64[], Float64[]
    n_points = 20
    count = 1:n_points
    
    for n in count
        push!(v_j, x[1])
        push!(v_a, x[2])
        x = A * x
    end
    
    fig2 = Figure(size=(800, 450))
    ax = Axis(fig2[1,1], title="Populația de urși - An Prost", 
              xlabel="Anul (N)", ylabel="Număr de urși")
    colors = Makie.wong_colors()
    
    scatterlines!(ax, count, v_j, label="Juvenili", color=colors[1], marker=:circle, linewidth=2)
    scatterlines!(ax, count, v_a, label="Adulți", color=colors[2], marker=:rect, linewidth=2)
    axislegend(ax, position=:lt)
    fig2
end

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000013
md"""
**Interpretare:** Populația scade spre zero. Observăm un **comportament tranzient**: la început, numărul de juvenili crește puțin înainte de a scădea. Acest lucru se datorează condiției inițiale (mulți adulți, puțini juvenili).
"""

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000014
md"""
### 🔄 Scenariul 3: Model Oscilant
"""

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000015
let
    A = [0.1 1.4; 0.4 0.2]
    x = [10.0, 50.0]
    v_j, v_a = Float64[], Float64[]
    n_points = 30
    count = 1:n_points
    
    for n in count
        push!(v_j, x[1])
        push!(v_a, x[2])
        x = A * x
    end
    
    fig3 = Figure(size=(800, 450))
    ax = Axis(fig3[1,1], title="Populația de urși - Model Oscilant", 
              xlabel="Anul (N)", ylabel="Număr de urși")
    colors = Makie.wong_colors()
    
    scatterlines!(ax, count, v_j, label="Juvenili", color=colors[1], marker=:circle, linewidth=2)
    scatterlines!(ax, count, v_a, label="Adulți", color=colors[2], marker=:rect, linewidth=2)
    axislegend(ax, position=:lt)
    fig3
end

# ╔═╡ d81d09a4-142a-4b1a-9af9-dfa9e3188fbb
md"""
### 📈 📉Scenariul 4: Ani buni și proști
"""

# ╔═╡ b76aecfd-8aa2-4f8f-8458-e04601a28830
let
    A = [0.65 0.5; 0.25 0.9]  # Matricea pentru An Bun
    B = [0.5 0.4; 0.1 0.8]   # Matricea pentru An Prost
    init = [100.0, 50.0]
    n_years = 10
    years = 0:n_years

    # Generăm secvențele de 10 ani
    # Secvența 1: Bun, Prost, Bun, Prost... (Se aplică A, apoi B)
    seq1 = [i % 2 == 1 ? A : B for i in 1:n_years]
    # Secvența 2: Prost, Bun, Prost, Bun... (Se aplică B, apoi A)
    seq2 = [i % 2 == 1 ? B : A for i in 1:n_years]

    # Funcție auxiliară pentru a calcula traiectoria populației
    function get_trajectory(init_pop, seq_matrices)
        traj = [copy(init_pop)]
        curr = copy(init_pop)
        for M in seq_matrices
            curr = M * curr
            push!(traj, copy(curr))
        end
        return traj
    end

    traj1 = get_trajectory(init, seq1)
    traj2 = get_trajectory(init, seq2)

    # Extragem datele pentru grafic
    J1 = [t[1] for t in traj1]; A1 = [t[2] for t in traj1]
    J2 = [t[1] for t in traj2]; A2 = [t[2] for t in traj2]

    fig = Figure(size=(1000, 500))
    ax = Axis(fig[1,1], 
              title="Impactul ordinii evenimentelor asupra populației (10 ani)",
              xlabel="Anul", 
              ylabel="Număr de urși",
              titlesize=20, xlabelsize=16, ylabelsize=16)
    
    colors = Makie.wong_colors()

    # Trasăm datele pentru Secvența 1 (Linii continue)
    scatterlines!(ax, years, J1, label="Juvenili (Bun → Prost)", color=colors[1], marker=:circle, linewidth=2)
    scatterlines!(ax, years, A1, label="Adulți (Bun → Prost)", color=colors[1], marker=:rect, linewidth=2)
    
    # Trasăm datele pentru Secvența 2 (Linii întrerupte)
    scatterlines!(ax, years, J2, label="Juvenili (Prost → Bun)", color=colors[2], marker=:circle, linestyle=:dash, linewidth=2)
    scatterlines!(ax, years, A2, label="Adulți (Prost → Bun)", color=colors[2], marker=:rect, linestyle=:dash, linewidth=2)

    axislegend(ax, position=:lt, fontsize=12)
    fig
end

# ╔═╡ 64bd52b4-2bd1-4633-85d6-a4869c80f87b
md"""
**Interpretare:** După cum se observă în grafic, deși ambele scenarii au exact același număr de ani buni și ani proști (5 din fiecare), **populația finală este diferită**. 

Dacă calculăm matricele compuse pe 2 ani, obținem:
$$C_1 = B \cdot A = \begin{pmatrix} 0.425 & 0.64 \\ 0.145 & 0.77 \end{pmatrix} \quad \text{vs.} \quad C_2 = A \cdot B = \begin{pmatrix} 0.375 & 0.66 \\ 0.13 & 0.76 \end{pmatrix}$$

Matematic, acest lucru confirmă că $B \cdot A \neq A \cdot B$. Biologic, acest lucru ne învață că **secvența evenimentelor de viață** (naștere, maturizare, deces) în interiorul unui an are un impact major asupra dinamicii pe termen lung a populației. Când construim modele Leslie, trebuie să fim foarte atenți la ordinea în care definim tranzițiile!
"""

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000016
md"""
## 🎯 4. Vectori și Valori Proprii

Până acum am determinat comportamentul pe termen lung *empiric*, iterând matricea de multe ori. Cum putem determina acest lucru *riguros*?

Răspunsul stă în căutarea unor vectori speciali care nu își schimbă direcția prin aplicarea matricei, ci doar sunt scalați. Căutăm vectori nenuli $\overrightarrow{v}$ și numere $\lambda$ astfel încât:
$$A \overrightarrow{v} = \lambda \overrightarrow{v}$$
Acești vectori $\overrightarrow{v}$ se numesc **vectori proprii**, iar numerele $\lambda$ se numesc **valori proprii**.

Dacă matricea $A$ are două valori proprii distincte $\lambda_1, \lambda_2$ cu vectorii proprii corespunzători $\overrightarrow{v_1}, \overrightarrow{v_2}$, atunci orice condiție inițială se poate scrie ca o combinație liniară:
$$\overrightarrow{x}_0 = c_1 \overrightarrow{v_1} + c_2 \overrightarrow{v_2}$$
Aplicând matricea de $N$ ori, obținem:
$$\overrightarrow{x}_N = A^N \overrightarrow{x}_0 = c_1 \lambda_1^N \overrightarrow{v_1} + c_2 \lambda_2^N \overrightarrow{v_2}$$

**Comportamentul pe termen lung** este dictat de **valoarea proprie dominantă** (cea cu modulul cel mai mare). Dacă $|\lambda_1| > |\lambda_2|$, termenul $\lambda_2^N$ devine neglijabil față de $\lambda_1^N$, iar populația va tinde să se alinieze cu direcția vectorului propriu $\overrightarrow{v_1}$.
"""

# ╔═╡ a1b2c3d4-0017-4a00-8000-000000000017
PlutoTeachingTools.danger(md"Nu orice matrice este diagonalizabilă (adică nu orice matrice are o bază de vectori proprii). Acest lucru se poate întâmpla dacă ecuația caracteristică nu are rădăcini reale distincte.")

# ╔═╡ a1b2c3d4-0018-4a00-8000-000000000018
md"""
### 🎛️ Explorare Interactivă: Convergența către Vectorul Propriu Dominant

Modifică numărul de iterații pentru a observa cum populația (orice vector inițial) se aliniază cu direcția vectorului propriu dominant.
"""

# ╔═╡ a1b2c3d4-0019-4a00-8000-000000000019
md"""
Numărul de iterații (ani): $N = {}$ $(@bind n_iter PlutoUI.Slider(1:50, default=10, show_value=true))
"""

# ╔═╡ a1b2c3d4-0020-4a00-8000-000000000020
let
    M = [0.65 0.5; 0.25 0.9] # Matricea pentru "Anul Bun"
    init = [50.0, 10.0]
    
    vals = eigvals(M)
    vecs = eigvecs(M)
    
    # Găsim valoarea proprie dominantă
    dom_idx = argmax(abs.(vals))
    dom_val = vals[dom_idx]
    dom_vec = vecs[:, dom_idx]
    dom_vec = dom_vec / sum(dom_vec) # Normalizăm pentru a obține raportul J/A
    
    xs, ys, ratios = Float64[], Float64[], Float64[]
    curr = copy(init)
    
    for i in 1:n_iter
        push!(xs, curr[1])
        push!(ys, curr[2])
        push!(ratios, curr[1] / (curr[1] + curr[2]))
        curr = M * curr
    end
    
    fig = Figure(size=(1000, 450))
    
    # Panoul 1: Spațiul stărilor
    ax1 = Axis(fig[1,1], title="Traiectoria în spațiul stărilor (J vs A)", 
               xlabel="Juvenili (J)", ylabel="Adulți (A)")
    lines!(ax1, xs, ys, color=:teal, linewidth=2)
    scatter!(ax1, xs, ys, color=:teal, markersize=8)
    scatter!(ax1, [init[1]], [init[2]], color=:red, markersize=12, label="Start")
    
    # Direcția vectorului propriu dominant
    max_val = max(maximum(xs), maximum(ys)) * 1.2
    lines!(ax1, [0, max_val*dom_vec[1]], [0, max_val*dom_vec[2]], 
           color=:orange, linestyle=:dash, linewidth=2, label="Direcția proprie dominantă")
    axislegend(ax1, position=:lt)
    
    # Panoul 2: Raportul pe iterații
    ax2 = Axis(fig[1,2], title="Raportul Juvenili / Total", 
               xlabel="Iterația (N)", ylabel="Raport J / (J+A)")
    lines!(ax2, 1:n_iter, ratios, color=:purple, linewidth=2)
    scatter!(ax2, 1:n_iter, ratios, color=:purple, markersize=8)
    hlines!(ax2, [dom_vec[1]], color=:orange, linestyle=:dash, 
            label="Raport asimptotic ($(round(dom_vec[1], digits=3)))")
    axislegend(ax2, position=:rt)
    
    fig
end

# ╔═╡ a1b2c3d4-0021-4a00-8000-000000000021
md"""
### 🧮 Calculul Valorilor și Vectorilor Proprii în Julia
Pentru matrice de ordin mai mare, ecuația caracteristică este $\det(\lambda I_n - M) = 0$. În Julia, putem folosi funcțiile din `LinearAlgebra`:
"""

# ╔═╡ a1b2c3d4-0022-4a00-8000-000000000022
let
    A = [1 2; 4 3]
    println("Matricea A:\n", A)
    println("\nValorile proprii (λ): ", eigvals(A))
    println("\nVectorii proprii (coloane):\n", eigvecs(A))
    println("\nInversa lui A:\n", inv(A))
end

# ╔═╡ a1b2c3d4-0023-4a00-8000-000000000023
md"""
## 📝 5. Concluzii și Sinteză

### Rezumat al conceptelor cheie

| Concept | Definiție / Semnificație | Exemplu din lecție |
| :--- | :--- | :--- |
| **📐 Aplicație Liniară** | O funcție care păstrează adunarea și înmulțirea cu scalari; reprezentată printr-o matrice. | Modelul Leslie de populație. |
| **🔄 Model Leslie** | Model matricial care descrie tranziția între clasele de vârstă (ex: Juvenili, Adulți). | Urșii americani ($J_{N+1}, A_{N+1}$). |
| **📈 Comportament Asimptotic** | Comportamentul pe termen lung al sistemului, independent de condițiile inițiale (după ce tranzienții dispar). | Creșterea exponențială sau declinul constant. |
| **🎯 Valoare Proprie ($\lambda$)** | Factorul de scalare al unui vector propriu. Determină viteza de creștere/descreștere pe o direcție. | $\lambda > 1 \implies$ creștere, $\lambda < 1 \implies$ declin. |
| **🧭 Vector Propriu ($\overrightarrow{v}$)** | O direcție care nu se rotește prin aplicarea matricei. Definește structura asimptotică a populației. | Raportul constant J/A pe termen lung. |

### Întrebări recapitulative
1. În modelul matricial al populației de urși, ce semnificație biologică au elementele de pe **diagonală** și cele **din afara diagonalei** ale matricei de tranziție?
2. De ce spunem că comportamentul pe termen lung al unui sistem liniar discret este dictat de **valoarea proprie dominantă**? Ce se întâmplă cu componentele corespunzătoare celorlalte valori proprii pe măsură ce numărul de iterații crește?
3. Cum ne ajută **vectorii proprii** să înțelegem "direcțiile naturale" de evoluție ale unui sistem dinamic și structura finală a populației?

---
**✨ Sfârșitul Cursului 7 ✨**
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╟─a1b2c3d4-0002-4a00-8000-000000000002
# ╟─a1b2c3d4-0003-4a00-8000-000000000003
# ╟─6ece3812-853a-4c42-a863-d9430a9bdcb3
# ╟─4865f5f6-92e9-4a74-94ec-9bd2aba29474
# ╠═a1b2c3d4-0004-4a00-8000-000000000004
# ╟─a1b2c3d4-0005-4a00-8000-000000000005
# ╟─a1b2c3d4-0006-4a00-8000-000000000006
# ╟─a1b2c3d4-0007-4a00-8000-000000000007
# ╟─a1b2c3d4-0008-4a00-8000-000000000008
# ╟─a1b2c3d4-0009-4a00-8000-000000000009
# ╟─a1b2c3d4-0010-4a00-8000-000000000010
# ╟─a1b2c3d4-0011-4a00-8000-000000000011
# ╟─a1b2c3d4-0012-4a00-8000-000000000012
# ╟─a1b2c3d4-0013-4a00-8000-000000000013
# ╟─a1b2c3d4-0014-4a00-8000-000000000014
# ╟─a1b2c3d4-0015-4a00-8000-000000000015
# ╟─d81d09a4-142a-4b1a-9af9-dfa9e3188fbb
# ╟─b76aecfd-8aa2-4f8f-8458-e04601a28830
# ╟─64bd52b4-2bd1-4633-85d6-a4869c80f87b
# ╟─a1b2c3d4-0016-4a00-8000-000000000016
# ╟─a1b2c3d4-0017-4a00-8000-000000000017
# ╟─a1b2c3d4-0018-4a00-8000-000000000018
# ╟─a1b2c3d4-0019-4a00-8000-000000000019
# ╟─a1b2c3d4-0020-4a00-8000-000000000020
# ╟─a1b2c3d4-0021-4a00-8000-000000000021
# ╠═a1b2c3d4-0022-4a00-8000-000000000022
# ╟─a1b2c3d4-0023-4a00-8000-000000000023
