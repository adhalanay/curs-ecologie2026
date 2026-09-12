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
    using DifferentialEquations
    using PlutoTeachingTools
    using DynamicalSystems
    using LinearAlgebra
    TableOfContents()
end

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 🧮 Cursul 8 — Sisteme Liniare și Analiza Stabilității

## 🎯 Obiectivele lecției
1. Înțelegerea intuitivă a **sistemelor decuplate** și a rolului vectorilor proprii.
2. Clasificarea **punctelor de echilibru** pentru sistemele liniare (discrete și continue).
3. Analiza comportamentului asimptotic: noduri, puncte șa, spirale și centre.
4. Explorarea interactivă a modului în care parametrii unei matrice dictează soarta sistemului dinamic.
"""

# ╔═╡ a1b2c3d4-0100-4a00-8000-000000000002
md"""
## 🧩 1. Fundamentul Teoretic: Sisteme Diagonale și Decuplarea

Pentru a înțelege *de ce* valorile proprii ($\lambda$) și vectorii proprii ($\vec{v}$) dictează comportamentul unui sistem complex, trebuie să începem cu cel mai simplu caz: **sistemele decuplate**, reprezentate prin matrice diagonale
```math
M=\begin{pmatrix}a & 0 \\ 0 & d\end{pmatrix}
```
cu $|a| \geq |d|$.

### 📐 1.1. Cazul Discret Decuplat
Să considerăm un sistem discret 2D unde variabila $x$ nu depinde de $y$ și invers:
$$\begin{pmatrix} x_{n+1} \\ y_{n+1} \end{pmatrix} = \begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} \begin{pmatrix} x_n \\ y_n \end{pmatrix}$$

Soluția este evidentă și complet independentă pe fiecare axă:
$$x_n = a^n \cdot x_0 \quad \text{și} \quad y_n = d^n \cdot y_0$$

**Reguli de evoluție:**
- Dacă $|a| > 1$, direcția $x$ **se amplifică** (creștere).
- Dacă $|a| < 1$, direcția $x$ **se atenuează** (descreștere spre 0).
- Dacă $a < 0$, sistemul **oscilează** pe acea axă (schimbă semnul la fiecare pas).
"""

# ╔═╡ a1b2c3d4-0101-4a00-8000-000000000003
let
    # Sistem diagonal: x crește (a=1.2), y scade (d=0.5)
    D = [1.2 0.0; 0.0 0.5]
    inits = [[1.0, 10.0], [5.0, 5.0], [10.0, 1.0], [-8.0, 6.0]]
    colors = Makie.wong_colors()
    
    fig = Figure(size=(1800, 550))
    ax = Axis(fig[1,1], title="Sistem Discret Diagonal (Decuplat)", 
              xlabel="x (λ₁ = 1.2 → Crește exponențial)", 
              ylabel="y (λ₂ = 0.5 → Scade spre 0)",
              aspect = DataAspect())
    
    for (i, init) in enumerate(inits)
        xs, ys = Float64[], Float64[]
        val = copy(init)
        for _ in 1:15
            push!(xs, val[1]); push!(ys, val[2])
            val = D * val
        end
        push!(xs, val[1]); push!(ys, val[2])
        lines!(ax, xs, ys, color=colors[i], linewidth=2)
        scatter!(ax, xs, ys, color=colors[i], markersize=8)
        scatter!(ax, [init[1]], [init[2]], markersize=12, color=colors[i], label="Start $i")
    end
    axislegend(ax, position=:lt, fontsize=10)
    fig
end

# ╔═╡ a1b2c3d4-0102-4a00-8000-000000000004
md"""
**Interpretare vizuală:** Deoarece sistemul este decuplat, traiectoriile nu se rotesc și nu se curbează diagonal. Ele se mișcă „în trepte” paralele cu axele, fiind atrase spre axa dominantă (în cazul nostru, axa $x$, unde valoarea proprie $1.2 > 1$).

### 🌊 1.2. Cazul Continuu Decuplat
Pentru ecuații diferențiale, un sistem diagonal arată astfel:
$$\begin{pmatrix} x' \\ y' \end{pmatrix} = \begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix}$$

Soluția implică funcții exponențiale continue:
$$x(t) = e^{at} \cdot x_0 \quad \text{și} \quad y(t) = e^{dt} \cdot y_0$$

Aici, semnul lui $a$ și $d$ dictează totul:
- Dacă $a > 0$, $x(t)$ crește exponențial.
- Dacă $a < 0$, $x(t)$ descrește spre 0.

### 🌉 1.3. Puntea către Sistemele Generale (Diagonalizarea)
Ce legătură are acest caz simplu cu o matrice oarecare $A$ (cu elemente pe toate pozițiile, unde $x$ depinde de $y$)? 

Teorema fundamentală a algebrei liniare ne spune că, dacă $A$ are valori și vectori proprii distincți, ea poate fi **diagonalizată**: $A = P D P^{-1}$.
- Matricea $D$ conține **valorile proprii** ($\lambda_1, \lambda_2$) pe diagonală.
- Matricea $P$ are ca **coloane vectorii proprii** ($\vec{v}_1, \vec{v}_2$).

**Concluzia majoră:** Dacă ne „schimbăm ochelarii” și privim sistemul din perspectiva unei noi baze de coordonate formată din vectorii proprii, matricea $A$ devine brusc o matrice diagonală $D$! 
Așadar, **valorile proprii** ne spun *cât de repede* crește/descrește sistemul, iar **vectorii proprii** ne indică *direcțiile geometrice* (axele) de-a lungul cărora are loc această evoluție decuplată.
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000005
md"""
---
## 🔄 2. Sisteme Liniare Discrete (Cazul General)

Revenind la un sistem discret general $\overrightarrow{X}_{n+1} = M \overrightarrow{X}_n$, comportamentul pe termen lung este dictat exclusiv de valorile proprii $\lambda_1, \lambda_2$ ale matricei $M$, exact ca în cazul diagonal, dar rotit în spațiu conform vectorilor proprii. 

Vom nota cu $\lambda$ valoarea proprie avînd modulul cel mai mare. Aceasta se numește valoarea **dominantă**.

**Regula de aur pentru sistemele discrete:**
- Dacă $|\lambda| > 1$, direcția corespunzătoare **se amplifică** (instabilitate).
- Dacă $|\lambda| < 1$, direcția corespunzătoare **se atenuează** (stabilitate).
- Dacă $\lambda < 0$, sistemul **oscilează** (alternează semnul) pe acea direcție.
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000006
md"""
### 🐻 2.1. Exemplu: Nod șa (Saddle Point) - Modelul Ursului (An Bun)

Matricea $M = \begin{pmatrix} 0.65 & 0.5 \\ 0.25 & 0.9 \end{pmatrix}$ are valorile proprii $\lambda_1 \approx 1.15$ și $\lambda_2 \approx 0.4$. Deoarece una este mai mare ca 1, iar cealaltă mai mică, originea este un **punct șa**.
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000007
let
    M = [0.65 0.5; 0.25 0.9]
    vals = eigvals(M)
    vecs = eigvecs(M)
    
    inits = [[10.0, 5.0], [10.0, -5.0], [5.0, 20.0], [-10.0, 10.0]]
    colors = Makie.wong_colors()
    
    fig = Figure(size=(1800, 450))
    ax = Axis(fig[1,1], title="Portret de fază: Punct șa (Saddle)", 
              xlabel="Juvenili (J)", ylabel="Adulți (A)")
    
    # Desenăm direcțiile proprii (eigenvectors)
    for i in 1:2
        v = vecs[:, i]
        lines!(ax, [-20v[1], 20v[1]], [-20v[2], 20v[2]], 
               color=:gray, linestyle=:dash, linewidth=1.5, label="Dir. proprie $(round(vals[i], digits=2))")
    end
    
    # Desenăm traiectoriile
    for (i, init) in enumerate(inits)
        xs, ys = Float64[], Float64[]
        val = copy(init)
        for _ in 1:15
            push!(xs, val[1]); push!(ys, val[2])
            val = M * val
        end
        push!(xs, val[1]); push!(ys, val[2])
        lines!(ax, xs, ys, color=colors[i], linewidth=2)
        scatter!(ax, xs, ys, color=colors[i], markersize=8)
        scatter!(ax, [init[1]], [init[2]], markersize=15, color=colors[i], label="Start")
    end
    
    axislegend(ax, position=:lt, fontsize=10)
    fig
end

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000008
md"""
**Interpretare:** Traiectoriile sunt atrase spre direcția proprie a valorii dominante ($\lambda_1 \approx 1.15$) și respinse de pe direcția valorii $\lambda_2 \approx 0.4$. Pe termen lung, populația va crește exponențial de-a lungul primei direcții.
"""

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000009
md"""
### 🥀 2.2. Nod Stabil (Stable Node) - Modelul Ursului (An Prost)

Matricea $M = \begin{pmatrix} 0.5 & 0.4 \\ 0.1 & 0.8 \end{pmatrix}$ are $\lambda_1 = 0.9$ și $\lambda_2 = 0.4$. Ambele sunt subunitare, deci originea este un **nod stabil** (toate traiectoriile converg la 0).
"""

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000010
let
    M = [0.5 0.4; 0.1 0.8]
    vals = eigvals(M)
    vecs = eigvecs(M)
    
    inits = [[10.0, 5.0], [10.0, -5.0], [5.0, 20.0], [-10.0, 10.0]]
    colors = Makie.wong_colors()
    
    fig = Figure(size=(1800, 450))
    ax = Axis(fig[1,1], title="Portret de fază: Nod Stabil", 
              xlabel="Juvenili (J)", ylabel="Adulți (A)")
    
    for i in 1:2
        v = vecs[:, i]
        lines!(ax, [-20v[1], 20v[1]], [-20v[2], 20v[2]], 
               color=:gray, linestyle=:dash, linewidth=1.5, label="Dir. proprie $(round(vals[i], digits=2))")
    end
    
    for (i, init) in enumerate(inits)
        xs, ys = Float64[], Float64[]
        val = copy(init)
        for _ in 1:20
            push!(xs, val[1]); push!(ys, val[2])
            val = M * val
        end
        push!(xs, val[1]); push!(ys, val[2])
        lines!(ax, xs, ys, color=colors[i], linewidth=2)
        scatter!(ax, xs, ys, color=colors[i], markersize=8)
    end
    scatter!(ax, [0], [0], markersize=20, color=:red, label="Echilibru (0,0)")
    axislegend(ax, position=:rt, fontsize=10)
    fig
end

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000011
md"""
### 🌀 2.3. Spirale și Oscilații (Valori proprii negative sau complexe)

Dacă o valoare proprie este **negativă** (ex: $\lambda = -0.6$), sistemul discret va oscila (sări de o parte și de alta a originii pe acea direcție). Dacă valorile proprii sunt **complexe**, traiectoria va descrie o **spirală**.
"""

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000012
let
    # Model oscilant cu valori proprii reale: -0.6 și 0.9
    M = [0.1 1.4; 0.4 0.2]
    vals = eigvals(M)
    
    inits = [[10.0, 5.0], [-10.0, 15.0]]
    colors = Makie.wong_colors()
    
    fig = Figure(size=(1800, 450))
    ax = Axis(fig[1,1], title="Comportament Oscilant (λ₁ = $(round(vals[1], digits=2)), λ₂ = $(round(vals[2], digits=2)))", 
              xlabel="Juvenili (J)", ylabel="Adulți (A)")
    
    for (i, init) in enumerate(inits)
        xs, ys = Float64[], Float64[]
        val = copy(init)
        for _ in 1:25
            push!(xs, val[1]); push!(ys, val[2])
            val = M * val
        end
        push!(xs, val[1]); push!(ys, val[2])
        lines!(ax, xs, ys, color=colors[i], linewidth=2)
        scatter!(ax, xs, ys, color=colors[i], markersize=8)
    end
    scatter!(ax, [0], [0], markersize=20, color=:red, label="Echilibru")
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000013
md"""
**Interpretare:** Deoarece $\lambda_1 = -0.6$ este negativă, componenta pe acea direcție proprie își schimbă semnul la fiecare pas, generând „zigzag-ul” (oscilația) traiectoriei înainte de a converge la echilibru (dictat de $\lambda_2 = 0.9 < 1$).
"""

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000014
md"""
---
## 🌊 3. Sisteme de Ecuații Diferențiale Liniare (Continue)

Pentru un sistem continuu $\dot{\overrightarrow{X}} = M \overrightarrow{X}$, soluția implică exponențiale de forma $e^{\lambda t}$. 

**Regula de aur pentru sistemele continue:**
- Dacă $Re(\lambda) > 0$, direcția **se amplifică** (instabilitate).
- Dacă $Re(\lambda) < 0$, direcția **se atenuează** (stabilitate).
- Dacă $\lambda$ este complex ($\alpha \pm i\beta$), sistemul oscilează (frecvența unghiulară este $\beta$), iar $\alpha$ dictează dacă spirala se strânge ($\alpha < 0$) sau se desfășoară ($\alpha > 0$).
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000015
md"""
### 🎛️ 3.1. Explorare Interactivă: Câmpul Vectorial 2D

Modifică parametrii $a, b, c, d$ ai matricei $M = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$ pentru a vedea cum se deformează câmpul vectorial și traiectoriile sistemului continuu $\dot{x} = ax+by, \dot{y} = cx+dy$.
"""

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000016
md"""
Parametrii matricei $M = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$:
"""

# ╔═╡ ddff5619-ce87-4ca0-bdf2-79a1633b68a1
md""" $a = {}$ $(@bind a_cont PlutoUI.Slider(-2.0:0.1:2.0, default=0.5, show_value=true)) """

# ╔═╡ d75b1821-4ce0-4d01-b1b4-e78a7971a996
md""" $b = {}$ $(@bind b_cont PlutoUI.Slider(-2.0:0.1:2.0, default=-1.0, show_value=true)) """

# ╔═╡ 76a8066e-9aa4-4086-bd41-346f1ca57212
md""" $c = {}$ $(@bind c_cont PlutoUI.Slider(-2.0:0.1:2.0, default=1.0, show_value=true)) """

# ╔═╡ f92728df-5870-4b73-a113-7f309e9fbc61
md""" $d = {}$ $(@bind d_cont PlutoUI.Slider(-2.0:0.1:2.0, default=0.5, show_value=true)) """

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000017
let
    M = [a_cont b_cont; c_cont d_cont]
    vals = eigvals(M)
    
    # Funcția pentru câmpul vectorial
    function field(x, y)
        return Point2f(a_cont*x + b_cont*y, c_cont*x + d_cont*y)
    end
    
    # Funcția pentru ODE (DynamicalSystems)
    function gen!(du, u, p, t)
        du[1] = a_cont*u[1] + b_cont*u[2]
        du[2] = c_cont*u[1] + d_cont*u[2]
        return nothing
    end
    
    fig = Figure(size=(1000, 500))
    
    # Panoul 1: Câmpul vectorial și o traiectorie
    ax1 = Axis(fig[1,1], title="Câmp vectorial și traiectorie", 
               xlabel="x", ylabel="y", aspect=DataAspect())
    
    xs = LinRange(-5, 5, 20)
    ys = LinRange(-5, 5, 20)
    streamplot!(ax1, field, xs, ys, colormap=:viridis, arrow_size=10, linewidth=1.5)
    
    # O traiectorie specifică
    initv = [3.0, 2.0]
    sisdiff = CoupledODEs(gen!, initv, ())
    Y, ti = trajectory(sisdiff, 10.0; Δt = 0.05)
    lines!(ax1, Y[:,1], Y[:,2], color=:red, linewidth=3, label="Traiectorie")
    scatter!(ax1, [initv[1]], [initv[2]], markersize=15, color=:red)
    scatter!(ax1, [0], [0], markersize=15, color=:black, marker=:xcross, label="Echilibru")
    
    # Panoul 2: Serii temporale
    ax2 = Axis(fig[1,2], title="Serii temporale", xlabel="Timp (t)", ylabel="Valoare")
    lines!(ax2, ti, Y[:,1], color=Makie.wong_colors()[1], linewidth=2, label="x(t)")
    lines!(ax2, ti, Y[:,2], color=Makie.wong_colors()[2], linewidth=2, label="y(t)")
    axislegend(ax2, position=:rt)
    
    # Text cu valorile proprii
    val_text = "λ₁ = $(round(vals[1], digits=2)), λ₂ = $(round(vals[2], digits=2))"
    Label(fig[0, :], text="Valori proprii: $val_text", tellwidth=false, fontsize=20, font=:bold)
    
    fig
end

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000018
md"""
**Experimentează cu glisoarele:**
- Încearcă să setezi $a=0, b=-1, c=1, d=0$. Ce obții? (Un **centru**, valori proprii $\pm i$).
- Setează $a=0.1, b=-1, c=1, d=0.1$. Ce se întâmplă cu spirala? (Devine **instabilă**, $\alpha > 0$).
"""

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000019
md"""
---
## 📝 4. Concluzii și Sinteză

### 📊 Tabel Rezumtiv: Geometria Echilibrelor Liniare

| Tipul Echilibrului | Condiție Valori Proprii (Discret $\lambda$) | Condiție Valori Proprii (Continu $\lambda$) | Comportament Geometric |
| :--- | :--- | :--- | :--- |
| **Nod Stabil** | $|\lambda_1| < 1, |\lambda_2| < 1$ (reale) | $Re(\lambda_1) < 0, Re(\lambda_2) < 0$ (reale) | Toate traiectoriile converg direct la origine. |
| **Nod Instabil** | $|\lambda_1| > 1, |\lambda_2| > 1$ (reale) | $Re(\lambda_1) > 0, Re(\lambda_2) > 0$ (reale) | Toate traiectoriile diverg de la origine. |
| **Punct șa (Saddle)** | Una $>1$, alta $<1$ în modul | Una $>0$, alta $<0$ (reale) | Atracție pe o direcție, respingere pe cealaltă. |
| **Spirală (Focar)** | Valori proprii complexe | Valori proprii complexe ($\alpha \pm i\beta$) | Traiectorii care se rotesc în jurul originii. |
| **Centru** (doar continuu) | N/A (specific sistemelor conservative) | $\lambda = \pm i\beta$ (pur imaginare) | Orbite închise (elipse) în jurul originii. |

### 💡 Întrebări de reflecție
1. **Discret vs. Continuu:** De ce o valoare proprie reală *negativă* (ex: $\lambda = -0.5$) generează oscilații în sistemul discret, dar nu și în cel continuu? (Indiciu: analizează forma soluției: $\lambda^n$ vs $e^{\lambda t}$).
2. **Centrul:** În sistemul continuu, un centru apare când valorile proprii sunt pur imaginare. De ce acest tip de echilibru este considerat „fragil” (structurally unstable) în comparație cu un nod sau o spirală?
3. **Modelul S-I:** În modelul epidemiologic S-I, matricea de tranziție are întotdeauna o valoare proprie $\lambda_1 = 1$. Ce semnificație biologică are această valoare proprie egală cu unitatea în contextul unei populații totale constante?

---
**✨ Sfârșitul Cursului 8 ✨**
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╟─a1b2c3d4-0100-4a00-8000-000000000002
# ╟─a1b2c3d4-0101-4a00-8000-000000000003
# ╟─a1b2c3d4-0102-4a00-8000-000000000004
# ╟─a1b2c3d4-0002-4a00-8000-000000000005
# ╟─a1b2c3d4-0003-4a00-8000-000000000006
# ╟─a1b2c3d4-0004-4a00-8000-000000000007
# ╟─a1b2c3d4-0005-4a00-8000-000000000008
# ╟─a1b2c3d4-0006-4a00-8000-000000000009
# ╟─a1b2c3d4-0007-4a00-8000-000000000010
# ╟─a1b2c3d4-0008-4a00-8000-000000000011
# ╟─a1b2c3d4-0009-4a00-8000-000000000012
# ╟─a1b2c3d4-0010-4a00-8000-000000000013
# ╟─a1b2c3d4-0011-4a00-8000-000000000014
# ╟─a1b2c3d4-0012-4a00-8000-000000000015
# ╟─a1b2c3d4-0013-4a00-8000-000000000016
# ╟─ddff5619-ce87-4ca0-bdf2-79a1633b68a1
# ╟─d75b1821-4ce0-4d01-b1b4-e78a7971a996
# ╟─76a8066e-9aa4-4086-bd41-346f1ca57212
# ╟─f92728df-5870-4b73-a113-7f309e9fbc61
# ╟─a1b2c3d4-0014-4a00-8000-000000000017
# ╟─a1b2c3d4-0015-4a00-8000-000000000018
# ╟─a1b2c3d4-0016-4a00-8000-000000000019
