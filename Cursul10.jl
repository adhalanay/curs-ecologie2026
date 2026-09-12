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
    using Makie, WGLMakie,CairoMakie
    using PlutoUI
    using DifferentialEquations
    using PlutoTeachingTools
    using DynamicalSystems
    using LinearAlgebra
    TableOfContents()
end

# ╔═╡ 211b3843-dbe3-4986-8c1b-1d1cd8d1e445
WGLMakie.activate!(; use_html_widgets = true)

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 🌌 Cursul 10 — Sisteme Conservative și Bifurcații Locale

**Durata estimată:** 100 minute

## 🎯 Obiectivele lecției
1. Înțelegerea **sistemelor conservative** și a conceptului de constantă de mișcare.
2. Analiza **bifurcației Hopf** și tranziția de la echilibru stabil la ciclu limită.
3. Clasificarea bifurcațiilor locale pentru **sisteme continue** (șa-nod, transcritică, trident).
4. Explorarea bifurcațiilor pentru **sisteme discrete** (fold, dublarea perioadei, Neimark-Sacker).
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000002
md"""
---
## ⚖️ 1. Sisteme Conservative și Centre (20 min)

Până acum, teorema Hartmann-Grobman ne-a permis să liniarizăm sistemele în jurul echilibrelor. Dar ce se întâmplă când valorile proprii sunt pur imaginare ($\lambda = \pm i\beta$)? Liniarizarea prezice un **centru**, dar termenii neliniari pot distruge această structură. 

Există totuși o clasă specială de sisteme unde centrele sunt garantate: **sistemele conservative**. Acestea posedă o funcție $H(x,y)$ (energia, sau o cantitate conservată) care rămâne constantă de-a lungul oricărei traiectorii ($\frac{dH}{dt} = 0$). 

> **Teoremă:** Dacă un sistem este conservativ și $(x^*, y^*)$ este un punct de extrem local (minim sau maxim) pentru $H$, atunci $(x^*, y^*)$ este un **centru** și toate orbitele din jurul lui sunt închise.
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000003
md"""
### 🦌 1.1. Modelul Lotka-Volterra (Prădător-Pradă)

Sistemul clasic Lotka-Volterra este conservativ. Cantitatea conservată este:
$$H(S,T) = c \ln S - dS - a T + b \ln T$$
Această funcție are un punct de extrem (maxim) în echilibrul netrivial, garantând că toate orbitele din jurul lui sunt închise (centre).
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000004
let
    h(x,y) = log(x) - x - y + log(y)
    f(x) = Point2f(x[1]*x[2] - x[1], -x[1]*x[2] + x[2])
    
    S = LinRange(0.1, 3.0, 100)
    T = LinRange(0.1, 3.0, 100)
    H = [h(s,t) for s ∈ S, t ∈ T]
    
    fig1 = Figure(size=(1000, 500))
    
    # Panoul 1: Suprafața 3D a funcției H
    ax1 = Axis3(fig1[1,1], title="Suprafața de energie H(S,T)", 
                xlabel="S (Pradă)", ylabel="T (Prădător)", zlabel="H")
    surface!(ax1, S, T, H, colormap=:viridis)
    contour3d!(ax1, S, T, H, levels=12, linewidth=1.5, color=:black)
    scatter!(ax1, [1.0], [1.0], [h(1.0,1.0)], markersize=15, color=:red, label="Echilibru (1,1)")
    axislegend(ax1, position=:lt)
    
    # Panoul 2: Portretul de fază
    ax2 = Axis(fig1[1,2], title="Portret de fază (Centre)", 
               xlabel="S (Pradă)", ylabel="T (Prădător)", aspect=DataAspect())
    streamplot!(ax2, f, S, T, colormap=:magma, linewidth=1.5, arrow_size=5)
    scatter!(ax2, [1.0], [1.0], markersize=12, color=:red, label="Echilibru (1,1)")
    axislegend(ax2, position=:lt)
    
    fig1
end

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000005
md"""
### 🕰️ 1.2. Pendulul Simplex și Orbite Homoclinice

Ecuația pendulului fără frecare este:
$$\begin{aligned} \theta' &= \omega \\ \omega' &= -\sin \theta \end{aligned}$$
Cantitatea conservată este energia totală: $E(\theta, \omega) = \frac{1}{2}\omega^2 - \cos\theta$.
- În $(0,0)$, energia este minimă $\implies$ **Centru**.
- În $(\pm\pi, 0)$, energia este maximă $\implies$ **Punct șa**.

Traiectoriile care pornesc dintr-un punct șa și se întorc în același punct șa (după un timp infinit) se numesc **orbite homoclinice**. Ele separă oscilațiile (librațiile) de rotațiile continue.
"""

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000006
let
    e(x,y) = 0.5 * y^2 - cos(x)
    g(x) = Point2f(x[2], -sin(x[1]))
    
    xs = LinRange(-2π, 2π, 200)
    ys = LinRange(-3.0, 3.0, 200)
    zs = [e(x,y) for x ∈ xs, y ∈ ys]
    
    fig2 = Figure(size=(1000, 500))
    
    ax21 = Axis3(fig2[1,1], title="Energia Pendulului E(θ,ω)", 
                 xlabel="θ (Unghi)", ylabel="ω (Viteză)", zlabel="E", elevation=1.2)
    surface!(ax21, xs, ys, zs, colormap=:viridis)
    contour3d!(ax21, xs, ys, zs, levels=12, linewidth=1.5, color=:black)
    scatter!(ax21, [0.0], [0.0], [e(0.0,0.0)], markersize=15, color=:red, label="Centru (0,0)")
    
    ax22 = Axis(fig2[1,2], title="Portret de fază (Orbite Homoclinice)", 
                xlabel="θ", ylabel="ω", aspect=DataAspect())
    streamplot!(ax22, g, xs, ys, colormap=:magma, linewidth=1.5, arrow_size=5)
    scatter!(ax22, [-2π, -π, 0, π, 2π], [0, 0, 0, 0, 0], markersize=10, color=:red, label="Puncte șa / Centre")
    axislegend(ax22, position=:lt)
    
    fig2
end

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000007
md"""
---
## 🌀 2. Bifurcația Poincaré-Andronov-Hopf (20 min)

Bifurcația Hopf este mecanismul principal prin care un echilibru stabil "naște" un ciclu limită. Aceasta apare când o pereche de valori proprii complexe conjugate traversează axa imaginară (partea reală trece de la negativ la pozitiv).

Să reluăm clarinetul lui Rayleigh, dar acum cu un parametru de control $c$ (frecarea):
$$\begin{aligned} x' &= v \\ v' &= -x - c(v^3 - v) \end{aligned}$$
Jacobianul în origine are valorile proprii $\lambda = \frac{c \pm \sqrt{c^2 - 4}}{2}$.
- Pentru $c > 0$: spirală stabilă.
- Pentru $c = 0$: centru (bifurcația).
- Pentru $c < 0$: spirală instabilă, dar apare un **ciclu limită stabil**.
"""

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000008
md"""
Parametrul de frecare: $c = {}$ $(@bind c_hopf PlutoUI.Slider(-1.0:0.05:1.0, default=0.0, show_value=true))
"""

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000009
let
    r(x) = Point2f(x[2], -x[1] - c_hopf * (x[2]^3 - x[2]))
    
    xs = LinRange(-2.5, 2.5, 100)
    ys = LinRange(-2.5, 2.5, 100)
    
    fig3 = Figure(size=(800, 500))
    ax3 = Axis(fig3[1,1], title="Bifurcația Hopf (c = $c_hopf)", 
               xlabel="x", ylabel="v", aspect=DataAspect())
    
    streamplot!(ax3, r, xs, ys, colormap=:magma, linewidth=1.5, arrow_size=5)
    scatter!(ax3, [0], [0], markersize=12, color=:red, label="Echilibru (0,0)")
    axislegend(ax3, position=:lt)
    
    fig3
end

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000010
md"""
**Interpretare:** Modificând glisorul, observați cum pentru $c > 0$ traiectoriile sunt absorbite de origine, iar pentru $c < 0$ originea devine instabilă, dar traiectoriile converg către o orbită închisă (ciclul limită).
"""

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000011
md"""
---
## 🛤️ 3. Bifurcații Locale în Sisteme Continue (20 min)

În dimensiuni superioare, bifurcațiile de codimensiune 1 (cele mai comune) pentru echilibrele sistemelor continue sunt:
1. **Șa-Nod (Saddle-Node):** Apare/dispare o pereche de echilibre (unul șa, unul nod).
2. **Transcritică:** Două echilibre se ciocnesc și fac schimb de stabilitate.
3. **Trident (Pitchfork):** Dintr-un echilibru apar alte două (simetrie spartă).

Acestea pot fi înțelese prin forme normale (ecuații simplificate care capturează esența locală).
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000012
md"""
Parametrul de bifurcație: $\mu = {}$ $(@bind μ_cont PlutoUI.Slider(-0.5:0.01:0.5, default=0.1, show_value=true))
"""

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000013
let
    # Formele normale pentru cele 3 bifurcații (cuplate cu o direcție stabilă y' = -y)
    s(x) = Point2f(μ_cont - x[1]^2, -x[2])       # Șa-Nod
    t1(x) = Point2f(μ_cont*x[1] - x[1]^2, -x[2]) # Transcritică
    t2(x) = Point2f(μ_cont*x[1] - x[1]^3, -x[2]) # Trident
    
    xs = LinRange(-2.0, 2.0, 100)
    ys = LinRange(-2.0, 2.0, 100)
    
    fig4 = Figure(size=(1200, 400))
    
    ax41 = Axis(fig4[1,1], title="Șa-Nod (μ = $μ_cont)", xlabel="x", ylabel="y")
    streamplot!(ax41, s, xs, ys, colormap=:viridis, linewidth=1.5, arrow_size=5)
    
    ax42 = Axis(fig4[1,2], title="Transcritică (μ = $μ_cont)", xlabel="x", ylabel="y")
    streamplot!(ax42, t1, xs, ys, colormap=:viridis, linewidth=1.5, arrow_size=5)
    
    ax43 = Axis(fig4[1,3], title="Trident (μ = $μ_cont)", xlabel="x", ylabel="y")
    streamplot!(ax43, t2, xs, ys, colormap=:viridis, linewidth=1.5, arrow_size=5)
    
    fig4
end

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000014
md"""
**Interpretare:** Urmăriți cum se modifică numărul și stabilitatea punctelor de echilibru (de-a lungul axei $x$, unde $y=0$) pe măsură ce $\mu$ trece prin zero.
"""

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000015
md"""
---
## 🔄 4. Bifurcații în Sisteme Discrete (30 min)

Pentru un sistem discret $\overrightarrow{x}_{n+1} = f(\overrightarrow{x}_n, \alpha)$, bifurcațiile apar când o valoare proprie a Jacobianului traversează cercul unitate în planul complex. Există trei tipuri principale:

1. **Fold (Tangentă):** $\lambda = 1$. (Analogul șa-nod/transcritică).
2. **Flip (Dublarea Perioadei):** $\lambda = -1$. Un punct fix stabil devine instabil și dă naștere unei orbite de perioadă 2.
3. **Neimark-Sacker (Torică):** $\lambda = e^{\pm i\theta}$ (complex cu modul 1). Analogul discret al bifurcației Hopf. Un punct fix stabil dă naștere unui **cerc invariant** (analogul ciclului limită).
"""

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000016
md"""
### 📉 4.1. Harta Hénon: Bifurcații Fold și Flip

Harta Hénon este un model clasic pentru sistemele discrete 2D:
$$\begin{aligned} x_{n+1} &= r - x_n^2 + 0.3 y_n \\ y_{n+1} &= 0.5 x_n \end{aligned}$$
"""

# ╔═╡ a1b2c3d4-0017-4a00-8000-000000000017
md"""
Parametrul $r$ (Fold / Flip): $r = {}$ $(@bind r_henon PlutoUI.Slider(-0.5:0.01:1.5, default=0.8, show_value=true))
"""

# ╔═╡ a1b2c3d4-0018-4a00-8000-000000000018
let
    function henon_rule(u, p, n)
        x, y = u
        r, b = p
        xn = r - x^2 + b * y
        yn = x
        return SVector(xn, yn)
    end
    
    p0 = [r_henon, 0.3]
    init_vals = [[0.0, 0.0], [0.5, 0.2], [-0.5, -0.2]]
    
    prbs = []
    total_time = 200
    for u in init_vals
        henon = DeterministicIteratedMap(henon_rule, u, p0)
        X, t = trajectory(henon, total_time)
        push!(prbs, (X, t))
    end
    
    fig5 = Figure(size=(1000, 450))
    Label(fig5[0,:], text="Harta Hénon (r = $r_henon)", fontsize=24, justification=:center, tellwidth=false)
    
    ax51 = Axis(fig5[1,1], title="Portret de fază (Discret)", xlabel="xₙ", ylabel="yₙ", aspect=DataAspect())
    ax52 = Axis(fig5[1,2], title="Serii temporale (xₙ)", xlabel="Iterația (n)", ylabel="xₙ")
    
    colors = Makie.wong_colors()
    for (i, (X, t)) in enumerate(prbs)
        scatter!(ax51, X[:,1], X[:,2], markersize=6, color=colors[i], label="Start $i")
        lines!(ax52, t, X[:,1], color=colors[i], linewidth=2)
    end
    axislegend(ax51, position=:lt)
    
    fig5
end

# ╔═╡ a1b2c3d4-0019-4a00-8000-000000000019
md"""
**Interpretare:** 
- Pentru $r \approx -0.2$, nu există puncte fixe (Bifurcație **Fold**).
- Pentru $r \approx 0.8$, punctul fix este stabil.
- Pentru $r > 1.0$, punctul fix devine instabil și apare o orbită de perioadă 2 (Bifurcație **Flip**).
"""

# ╔═╡ a1b2c3d4-0020-4a00-8000-000000000020
md"""
### 🪐 4.2. Bifurcația Neimark-Sacker (Cercul Invariant)

Pentru a observa o bifurcație Neimark-Sacker, folosim o hartă normală în coordonate polare (transformată în carteziene):
$$\begin{aligned} x_{n+1} &= a x_n - b y_n - x_n(x_n^2 + y_n^2) \\ y_{n+1} &= b x_n + a y_n - y_n(x_n^2 + y_n^2) \end{aligned}$$
Când $a^2 + b^2 = 1$ (și $b \neq 0$), valorile proprii sunt pe cercul unitate.
"""

# ╔═╡ a1b2c3d4-0021-4a00-8000-000000000021
md"""
Parametrul $a$ (controlul stabilității): $a = {}$ $(@bind r_ns PlutoUI.Slider(0.5:0.01:1.2, default=0.9, show_value=true))
"""

# ╔═╡ a1b2c3d4-0022-4a00-8000-000000000022
let
    function ns_rule(u, p, n)
        x, y = u
        a, b = p
        r2 = x^2 + y^2
        xn = a * x - b * y - x * r2
        yn = b * x + a * y - y * r2
        return SVector(xn, yn)
    end
    
    p0 = [r_ns, 0.5] # b=0.5 fixează unghiul de rotație
    init_vals = [[0.1, 0.1], [0.5, 0.0], [0.0, 0.8]]
    
    prbs = []
    total_time = 1000
    for u in init_vals
        ns_map = DeterministicIteratedMap(ns_rule, u, p0)
        X, t = trajectory(ns_map, total_time)
        push!(prbs, (X, t))
    end
    
    fig6 = Figure(size=(1000, 450))
    Label(fig6[0,:], text="Bifurcația Neimark-Sacker (a = $r_ns)", fontsize=24, justification=:center, tellwidth=false)
    
    ax61 = Axis(fig6[1,1], title="Portret de fază", xlabel="xₙ", ylabel="yₙ", aspect=DataAspect())
    ax62 = Axis(fig6[1,2], title="Serii temporale (xₙ)", xlabel="Iterația (n)", ylabel="xₙ")
    
    colors = Makie.wong_colors()
    for (i, (X, t)) in enumerate(prbs)
        scatter!(ax61, X[:,1], X[:,2], markersize=3, color=colors[i])
        lines!(ax62, t, X[:,1], color=colors[i], linewidth=1.5)
    end
    
    # Adăugăm cercul unitate ca referință
    circle_x = [cos(θ) for θ in LinRange(0, 2π, 100)]
    circle_y = [sin(θ) for θ in LinRange(0, 2π, 100)]
    lines!(ax61, circle_x, circle_y, linestyle=:dash, color=:black, linewidth=2, label="Cercul Invariant")
    axislegend(ax61, position=:lt)
    
    fig6
end

# ╔═╡ a1b2c3d4-0023-4a00-8000-000000000023
md"""
**Interpretare:** 
- Pentru $a < 1$ (ex: $0.9$), originea este un punct fix stabil (toate traiectoriile converg în centru).
- Pentru $a > 1$ (ex: $1.1$), originea devine instabilă, dar traiectoriile converg către un **cerc invariant** (analogul discret al ciclului limită). Aceasta este bifurcația Neimark-Sacker.
"""

# ╔═╡ a1b2c3d4-0024-4a00-8000-000000000024
md"""
---
## 📝 5. Concluzii și Sinteză (10 min)

### 📊 Tabel Rezumtiv: Bifurcații Locale

| Tipul Bifurcației | Condiție Valori Proprii (Continuu $\lambda$) | Condiție Valori Proprii (Discret $\lambda$) | Ce se întâmplă? |
| :--- | :--- | :--- | :--- |
| **Șa-Nod / Fold** | $\lambda = 0$ | $\lambda = 1$ | Apare/dispare o pereche de echilibre (sau puncte fixe). |
| **Transcritică** | $\lambda = 0$ | $\lambda = 1$ | Două echilibre se ciocnesc și schimbă stabilitatea. |
| **Trident** | $\lambda = 0$ | $\lambda = 1$ | Un echilibru se sparge în trei (simetrie). |
| **Hopf / Neimark-Sacker** | $\lambda = \pm i\beta$ | $\lambda = e^{\pm i\theta}$ | Echilibrul devine instabil, apare un ciclu limită (sau cerc invariant). |
| **Flip (Dublarea Perioadei)**| N/A (specific continuu) | $\lambda = -1$ | Punctul fix devine instabil, apare o orbită de perioadă 2. |

### 💡 Întrebări de reflecție
1. De ce teorema Hartmann-Grobman nu ne poate garanta existența unui centru pentru un sistem neliniar general, chiar dacă liniarizarea prezice valori proprii $\pm i\beta$? Cum ne ajută sistemele conservative să ocolim această problemă?
2. Care este diferența fundamentală din punct de vedere geometric între un ciclu limită (apărut prin bifurcația Hopf în sisteme continue) și un cerc invariant (apărut prin bifurcația Neimark-Sacker în sisteme discrete)?
3. În modelul pendulului, ce semnificație fizică au orbitele homoclinice și cum separă ele cele două regimuri dinamice majeri ale sistemului?

---
**✨ Sfârșitul Cursului 10 ✨**
"""

# ╔═╡ Cell order:
# ╠═a1b2c3d4-0000-4a00-8000-000000000000
# ╠═211b3843-dbe3-4986-8c1b-1d1cd8d1e445
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╟─a1b2c3d4-0002-4a00-8000-000000000002
# ╟─a1b2c3d4-0003-4a00-8000-000000000003
# ╠═a1b2c3d4-0004-4a00-8000-000000000004
# ╟─a1b2c3d4-0005-4a00-8000-000000000005
# ╠═a1b2c3d4-0006-4a00-8000-000000000006
# ╟─a1b2c3d4-0007-4a00-8000-000000000007
# ╟─a1b2c3d4-0008-4a00-8000-000000000008
# ╟─a1b2c3d4-0009-4a00-8000-000000000009
# ╟─a1b2c3d4-0010-4a00-8000-000000000010
# ╟─a1b2c3d4-0011-4a00-8000-000000000011
# ╟─a1b2c3d4-0012-4a00-8000-000000000012
# ╟─a1b2c3d4-0013-4a00-8000-000000000013
# ╟─a1b2c3d4-0014-4a00-8000-000000000014
# ╟─a1b2c3d4-0015-4a00-8000-000000000015
# ╟─a1b2c3d4-0016-4a00-8000-000000000016
# ╟─a1b2c3d4-0017-4a00-8000-000000000017
# ╟─a1b2c3d4-0018-4a00-8000-000000000018
# ╟─a1b2c3d4-0019-4a00-8000-000000000019
# ╟─a1b2c3d4-0020-4a00-8000-000000000020
# ╟─a1b2c3d4-0021-4a00-8000-000000000021
# ╟─a1b2c3d4-0022-4a00-8000-000000000022
# ╟─a1b2c3d4-0023-4a00-8000-000000000023
# ╟─a1b2c3d4-0024-4a00-8000-000000000024
