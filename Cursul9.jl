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
# 🧮 Cursul 9 — Liniarizarea Sistemelor Neliniare și Stabilitatea Locală

**Durata estimată:** 100 minute

## 🎯 Obiectivele lecției
1. Înțelegerea conceptului de **aproximare liniară** (plan tangent) pentru funcții de mai multe variabile.
2. Calculul și interpretarea **matricei Jacobiene** și a **Teoremei Hartmann-Grobman**.
3. Analiza cazurilor limită unde liniarizarea eșuează (valori proprii nule sau pur imaginare) și conceptul de **stabilitate structurală**.
4. Extinderea analizei de stabilitate locală către **sistemele discrete neliniare**.
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000002
md"""
## 📐 1. Aproximarea Liniară și Matricea Jacobiană (25 min)

Pentru a studia stabilitatea sistemelor neliniare, vom folosi aceeași intuiție pe care o avem din analiza matematică de o singură variabilă: **aproximarea unei funcții complexe cu una liniară** (mai ușor de studiat) în vecinătatea unui punct de interes.
"""



# ╔═╡ b574946d-dba3-43c2-9d90-011bc8e36f63
md"""
### 📏 1.1. Reamintire: Liniarizarea în 1D
Pentru o funcție scalară $y = f(x)$, aproximarea liniară în jurul unui punct $x_0$ este dată de ecuația tangentei:
$$f(x) \approx f(x_0) + f'(x_0)(x - x_0)$$

Dacă notăm variația (perturbația) cu $\Delta x = x - x_0$ și $\Delta f = f(x) - f(x_0)$, obținem relația fundamentală:
$$\Delta f \approx f'(x_0) \cdot \Delta x$$
**Semnificația:** Derivata $f'(x_0)$ acționează ca un **factor de scalare**. Ea ne spune cu ce se înmulțește o mică perturbație $\Delta x$ pentru a obține noua valoare a funcției.
"""

# ╔═╡ 8daab4b5-97af-41f2-b61c-32e6790bd036
md"""
### 🏔️ 1.2. Funcții scalare de 2 variabile: Planul Tangent
Când trecem la o funcție de două variabile $z = f(x, y)$, graficul nu mai este o curbă, ci o **suprafață**. Pentru a o aproxima liniar într-un punct $(x_0, y_0)$, nu mai este suficientă o singură "pantă". Avem nevoie de pantele pe ambele direcții: $x$ și $y$.

Acestea sunt date de **derivatele parțiale**:
- $\frac{\partial f}{\partial x}$: cât de repede se schimbă $f$ dacă ne mișcăm doar pe direcția $x$ (ținând $y$ constant).
- $\frac{\partial f}{\partial y}$: cât de repede se schimbă $f$ dacă ne mișcăm doar pe direcția $y$ (ținând $x$ constant).

Aproximarea liniară (planul tangent) se scrie astfel:
$$\Delta z \approx \frac{\partial f}{\partial x}\bigg|_{(x_0,y_0)} \Delta x + \frac{\partial f}{\partial y}\bigg|_{(x_0,y_0)} \Delta y$$
"""

# ╔═╡ 65f55b39-d3d6-49ae-9ce2-b541c90d064d
let
    # Funcția scalară: un paraboloid (suprafață)
    f(x, y) = 5 - x^2 - y^2
    
    # Punctul de contact
    x0, y0 = 1.0, 1.0
    z0 = f(x0, y0) # z0 = 3
    
    # Derivatele parțiale în punctul (1,1)
    df_dx = -2 * x0  # -2
    df_dy = -2 * y0  # -2
    
    # Ecuația planului tangent: z = z0 + df_dx*(x-x0) + df_dy*(y-y0)
    # z = 3 - 2(x-1) - 2(y-1) => z = 7 - 2x - 2y
    plan_tangent(x, y) = 7 - 2x - 2y

    xs = LinRange(-2.5, 2.5, 50)
    ys = LinRange(-2.5, 2.5, 50)
    
    zs_surf = [f(x, y) for x in xs, y in ys]
    zs_plan = [plan_tangent(x, y) for x in xs, y in ys]

    fig = Figure(size=(800, 600))
    ax = Axis3(fig[1, 1], 
               title="Aproximarea Liniară în 2D: Planul Tangent", 
               xlabel="x", ylabel="y", zlabel="z = f(x,y)",
               aspect=:data)
    
    # Desenăm suprafața originală (neliniară)
    surface!(ax, xs, ys, zs_surf, colormap=:viridis, alpha=0.7, label="Suprafața z = 5 - x² - y²")
    
    # Desenăm planul tangent (liniar)
    surface!(ax, xs, ys, zs_plan, color=:red, alpha=0.5, label="Planul Tangent (Aproximarea)")
    
    # Punctul de contact
    scatter!(ax, [x0], [y0], [z0], markersize=15, color=:black, label="Punctul de contact (1, 1, 3)")
    
    axislegend(ax, position=:lt)
    fig
end

# ╔═╡ 183bb2b5-fb19-4887-a535-409b441f2936
md"""
**Interpretare vizuală:** În vecinătatea imediată a punctului negru, suprafața roșie (planul) și suprafața colorată (funcția originală) sunt practic de nedistinguit. Cu cât ne îndepărtăm de punctul de contact, cu atât aproximarea liniară devine mai puțin precisă.
"""

# ╔═╡ a48a49bf-ac5a-4361-9ded-8c1d056e957c
md"""
### 🧭 1.3. Câmpuri vectoriale și Matricea Jacobiană
Un sistem dinamic 2D este definit de un câmp vectorial $\vec{F}(x,y) = \begin{pmatrix} f(x,y) \\ g(x,y) \end{pmatrix}$. 

Avem acum **două** funcții ($f$ și $g$), fiecare cu câte două variabile. Pentru a aproxima liniar acest câmp vectorial în jurul unui punct $(x_0, y_0)$, trebuie să construim planul tangent pentru *ambele* componente simultan:

$$\begin{aligned} \Delta f &\approx \frac{\partial f}{\partial x} \Delta x + \frac{\partial f}{\partial y} \Delta y \\ \Delta g &\approx \frac{\partial g}{\partial x} \Delta x + \frac{\partial g}{\partial y} \Delta y \end{aligned}$$

Pentru a scrie acest sistem compact, folosim **algebra liniară**. Extragem coeficienții (derivatele parțiale) într-o matrice, numită **Matricea Jacobiană** ($J$):

$$\begin{pmatrix} \Delta f \\ \Delta g \end{pmatrix} \approx \begin{pmatrix} \frac{\partial f}{\partial x} & \frac{\partial f}{\partial y} \\ \frac{\partial g}{\partial x} & \frac{\partial g}{\partial y} \end{pmatrix}_{(x_0,y_0)} \begin{pmatrix} \Delta x \\ \Delta y \end{pmatrix} \implies \Delta \vec{F} \approx J \cdot \Delta \vec{x}$$

**Concluzie intermediară:** Matricea Jacobiană este echivalentul multidimensional al derivatei. Ea "capturează" toate ratele de schimbare locale ale sistemului.
"""

# ╔═╡ c36059dc-75b1-4bbb-8ea9-0e38b6987d2a
md"""
### 🌉 1.4. De la Liniarizare la Sisteme Dinamice
Cum folosim acest instrument pentru ecuații diferențiale? Fie un sistem neliniar:
$$\vec{x}' = \vec{F}(\vec{x})$$

Presupunem că $\vec{x}^*$ este un **punct de echilibru**, deci $\vec{F}(\vec{x}^*) = \vec{0}$ (sistemul stă pe loc).
Vrem să studiem ce se întâmplă dacă sistemul este ușor perturbat de la echilibru. Definim perturbația $\vec{u}$:
$$\vec{u} = \vec{x} - \vec{x}^* \implies \vec{x} = \vec{x}^* + \vec{u}$$

Derivăm în raport cu timpul (deoarece $\vec{x}^*$ este constant, $\vec{x}' = \vec{u}'$) și aplicăm aproximarea liniară (Jacobianul):
$$\vec{u}' = \vec{F}(\vec{x}^* + \vec{u}) \approx \vec{F}(\vec{x}^*) + J(\vec{x}^*) \cdot \vec{u}$$

Deoarece $\vec{F}(\vec{x}^*) = \vec{0}$, ecuația devine un **sistem liniar cu coeficienți constanți**:
$$\vec{u}' = J(\vec{x}^*) \cdot \vec{u}$$

Am reușit să reducem studiul unui sistem neliniar complex la studiul unui sistem liniar simplu, guvernat de matricea $J$!
"""

# ╔═╡ 72e6323d-9ad8-4f0b-bb75-1f72fdc044e8
md"""
### 🛡️ 1.5. Teorema Hartmann-Grobman: Când este validă aproximarea?
Există o capcană: aproximarea liniară este valabilă *doar local* și *doar dacă echilibrul este hiperbolic*.
> **Teorema Hartmann-Grobman:** Fie $\vec{x}^*$ un punct de echilibru. Dacă **nici o valoare proprie** a matricei Jacobiene $J(\vec{x}^*)$ nu are partea reală egală cu zero ($Re(\lambda) \neq 0$), atunci comportamentul local (calitativ) al sistemului neliniar în jurul lui $\vec{x}^*$ este **identic** cu cel al sistemului liniarizat $\vec{u}' = J \vec{u}$.

**Ce înseamnă $Re(\lambda) \neq 0$?**
- Dacă $Re(\lambda) < 0$, traiectoriile sunt atrase spre echilibru (stabilitate).
- Dacă $Re(\lambda) > 0$, traiectoriile sunt respinse (instabilitate).
- Dacă $Re(\lambda) = 0$ (valori proprii nule sau pur imaginare, ex: $\lambda = \pm i$), teorema **nu se aplică**. Termenii neliniari (de ordin superior, pe care i-am ignorat la liniarizare) devin dominanți și pot schimba complet comportamentul (ex: un centru liniar poate deveni o spirală stabilă sau instabilă în sistemul neliniar).
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000003
md"""
---
## 🦌 2. Exemple de Analiză a Stabilității Locale (25 min)

Să aplicăm procedura standard pe două modele clasice:
1. Calculăm punctele de echilibru.
2. Calculăm matricea Jacobiană în acele puncte.
3. Extragem valorile proprii și clasificăm echilibrul (nod, spirală, șa).
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000004
md"""
### 🎷 2.1. Oscilatorul Rayleigh

Ecuația oscilatorului cu frecare neliniară este:
$$\begin{aligned} x' &= v \\ v' &= -x - (v^3 - v) \end{aligned}$$

Singurul punct de echilibru este $(0,0)$. Matricea Jacobiană evaluată în origine este:
$$J(0,0) = \begin{pmatrix} 0 & 1 \\ -1 & 1 \end{pmatrix}$$
Polinomul caracteristic este $\lambda^2 - \lambda + 1 = 0$, cu rădăcinile $\lambda = \frac{1 \pm i\sqrt{3}}{2}$. Deoarece partea reală este strict pozitivă ($\alpha = 0.5 > 0$), originea este o **spirală instabilă**.
"""

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000005
let
    f1(x) = Point2f(x[2], -x[1] - (x[2]^3 - x[2]))
    fig1 = Figure(size=(800, 500))
    ax1 = Axis(fig1[1,1], title="Oscilatorul Rayleigh (Spirală Instabilă)", 
               xlabel="x", ylabel="v")
    xs1 = LinRange(-2, 2, 20)
    ys1 = LinRange(-2, 2, 20)
    streamplot!(ax1, f1, xs1, ys1, colormap=:magma, linewidth=1.5, arrow_size=5)
    scatter!(ax1, [0], [0], markersize=15, color=:red, label="Echilibru (0,0)")
    axislegend(ax1, position=:rb)
    fig1
end

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000006
md"""
### 🦌 2.2. Competiția Căprioare-Elani

Modelul de competiție pentru resurse este:
$$\begin{aligned} D' &= 3D - 2MD - D^2 \\ M' &= 2M - DM - M^2 \end{aligned}$$

Punctul de echilibru netrivial (coexistență) este $(D^*, M^*) = (1, 1)$. Matricea Jacobiană în acest punct este:
$$J(1,1) = \begin{pmatrix} -1 & -2 \\ -1 & -1 \end{pmatrix}$$
Polinomul caracteristic este $\lambda^2 + 2\lambda - 1 = 0$, cu rădăcinile $\lambda = -1 \pm \sqrt{2}$. Deoarece $\sqrt{2} > 1$, avem o valoare proprie pozitivă și una negativă. Echilibrul este un **punct șa (instabil)**, ceea ce înseamnă că speciile nu pot coexista pe termen lung; una va elimina cealaltă în funcție de condițiile inițiale.
"""

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000007
let
    f2(x) = Point2f(3x[1] - 2x[1]x[2] - x[1]^2, 2x[2] - x[1]x[2] - x[2]^2)
    fig2 = Figure(size=(800, 500))
    ax2 = Axis(fig2[1,1], title="Competiția Căprioare-Elani (Punct șa)", 
               xlabel="Căprioare (D)", ylabel="Elani (M)")
    xs2 = LinRange(0, 3, 20)
    ys2 = LinRange(0, 3, 20)
    streamplot!(ax2, f2, xs2, ys2, colormap=:viridis, linewidth=1.5, arrow_size=5)
    scatter!(ax2, [1], [1], markersize=15, color=:red, label="Echilibru (1,1)")
    scatter!(ax2, [0], [0], markersize=15, color=:black, marker=:xcross, label="Origine")
    axislegend(ax2, position=:lt)
    fig2
end

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000008
md"""
---
## ⚠️ 3. Când Eșuează Teorema Hartmann-Grobman (25 min)

Teorema Hartmann-Grobman are o condiție critică: **partea reală a valorilor proprii trebuie să fie strict diferită de zero**. Ce se întâmplă când această condiție este încălcată?
"""

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000009
md"""
### 📏 3.1. Valoarea Proprie Zero și Stabilitatea Structurală

Să considerăm sistemul liniar (care este propria sa liniarizare):
$$\begin{aligned} x' &= x - 2y \\ y' &= 3x - 6y \end{aligned}$$
Matricea sistemului are valorile proprii $\lambda_1 = 0$ și $\lambda_2 = -5$. Deoarece avem o valoare proprie zero, sistemul nu are un echilibru izolat, ci o **linie întreagă de echilibre** (direcția $y = 0.5x$).
"""

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000010
begin
    f3(x) = Point2f(x[1] - 2x[2], 3x[1] - 6x[2])
    fig3 = Figure(size=(800, 500))
    ax3 = Axis(fig3[1,1], title="Valoare Proprie Zero (ϵ = 0)", 
               xlabel="x", ylabel="y")
    xs3 = LinRange(-2, 2, 20)
    ys3 = LinRange(-2, 2, 20)
    streamplot!(ax3, f3, xs3, ys3, colormap=:magma, linewidth=1.5, arrow_size=5)
    lines!(ax3, [-2, 2], [-1, 1], color=:blue, linestyle=:dash, linewidth=2, label="Linie de echilibre")
    scatter!(ax3, [0], [0], markersize=15, color=:red, label="Origine")
    axislegend(ax3, position=:lt)
    fig3
end

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000011
md"""
Acest tip de sistem este "fragil". Dacă perturbăm sistemul cu un mic parametru $\epsilon$, linia de echilibre se distruge. Modifică $\epsilon$ pentru a vedea cum se schimbă calitativ portretul de fază.
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000012
md""" Perturbarea sistemului: $\epsilon = {}$ $(@bind ϵ PlutoUI.Slider(-0.6:0.01:0.6, default=0.0, show_value=true)) """

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000013
let
    f4(x) = Point2f(x[1] - 2x[2], (3+ϵ)x[1] - 6x[2])
    fig4 = Figure(size=(800, 500))
    ax4 = Axis(fig4[1,1], title="Sistem Perturbat (ϵ = $ϵ)", 
               xlabel="x", ylabel="y")
    streamplot!(ax4, f4, xs3, ys3, colormap=:magma, linewidth=1.5, arrow_size=5)
    scatter!(ax4, [0], [0], markersize=15, color=:red, label="Echilibru (0,0)")
    axislegend(ax4, position=:lt)
    fig4
end

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000014
md"""
**Interpretare:** Pe măsură ce $\epsilon$ trece de la negativ la pozitiv, originea (care devine singurul echilibru) se transformă brusc dintr-un punct șa într-o spirală stabilă. Sistemele care **își păstrează comportamentul calitativ** la mici perturbații se numesc **structural stabile**. Modelele biologice realiste trebuie să fie structural stabile, altfel orice mică eroare de măsurare a parametrilor ar schimba complet predicțiile.
"""

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000015
md"""
### 🌀 3.2. Valori Proprii Pur Imaginare (Centru vs. Spirală)

Un alt caz unde Hartmann-Grobman eșuează este când valorile proprii sunt pur imaginare ($\lambda = \pm i\beta$). Liniarizarea prezice un **centru** (orbite închise). Totuși, termenii neliniari (de ordin superior) pot face ca traiectoriile să spiraleze lent spre interior (centru stabil) sau spre exterior (centru instabil). Liniarizarea **nu ne poate spune** care este situația reală!
"""

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000016
let
    f(x,y) = (x^2+y^2)/(1+x^2+y^2)
    f5(x) = Point2f(-x[2] + f(x[1],x[2])*x[1], x[1] + f(x[1],x[2])*x[2])
    fig5 = Figure(size=(800, 500))
    ax5 = Axis(fig5[1,1], title="Valori Proprii ±i (Termenii neliniari decid)", 
               xlabel="x", ylabel="y")
    xs5 = LinRange(-2, 2, 20)
    ys5 = LinRange(-2, 2, 20)
    streamplot!(ax5, f5, xs5, ys5, colormap=:magma, linewidth=1.5, arrow_size=5)
    scatter!(ax5, [0], [0], markersize=15, color=:red, label="Origine (Liniarizare: Centru)")
    axislegend(ax5, position=:lt)
    fig5
end

# ╔═╡ a1b2c3d4-0017-4a00-8000-000000000017
md"""
**Interpretare:** Deși Jacobianul în origine are $\lambda = \pm i$ (ceea ce ar sugera un centru), termenul neliniar $f(x,y)$ acționează ca o "frecare negativă" care împinge traiectoriile spre exterior, transformând originea într-o **spirală instabilă**.
"""

# ╔═╡ a1b2c3d4-0018-4a00-8000-000000000018
md"""
---
## 🔄 4. Sisteme Discrete Neliniare (20 min)

Pentru un sistem discret neliniar $\overrightarrow{x}_{n+1} = \overrightarrow{f}(\overrightarrow{x}_n)$, analiza stabilității locale în jurul unui punct fix $\overrightarrow{x}^*$ folosește aceeași unealtă: **matricea Jacobiană** $J$ evaluată în $\overrightarrow{x}^*$.

**Regula de stabilitate:**
- Dacă **toate** valorile proprii ale lui $J$ au modulul subunitar ($|\lambda| < 1$), echilibrul este **stabil**.
- Dacă **cel puțin o** valoare proprie are $|\lambda| > 1$, echilibrul este **instabil**.

*(Notă: În practică, se pot folosi criteriile lui Jury pentru a verifica această condiție direct din urma și determinantul lui $J$, fără a calcula explicit valorile proprii).*
"""

# ╔═╡ a1b2c3d4-0019-4a00-8000-000000000019
md"""
### 🐛 4.1. Modelul Gazdă-Parazitoid (Nicholson-Bailey)

Acest model descrie interacțiunea dintre o gazdă ($H$) și un parazitoid ($P$):
$$\begin{aligned} H_{n+1} &= k \cdot e^{-a P_n} \cdot H_n \\ P_{n+1} &= c \cdot (1 - e^{-a P_n}) \cdot H_n \end{aligned}$$

Punctul de echilibru netrivial (coexistență) există, iar matricea Jacobiană în acest punct are urma și determinantul:
$$\text{tr} J = 1 + \frac{\ln k}{k-1}, \quad \det J = \ln k + \frac{\ln k}{k-1}$$
Deoarece $k > 1$ (rata de reproducere a gazdei), se poate demonstra matematic că $\det J > 1$ întotdeauna. Conform condițiilor de stabilitate pentru sisteme discrete, acest lucru implică faptul că **cel puțin o valoare proprie are $|\lambda| > 1$**.
"""

# ╔═╡ a1b2c3d4-0020-4a00-8000-000000000020
let
    a = 0.005
    k = 1.05
    c = 3.0
    fn(x,y) = exp(-a*y)
    function nicholson(u,p,n)
        x, y = u
        k,c = p
        xn = k * fn(x,y) * x
        yn = c * (1 - fn(x,y)) * x
        return SVector(xn, yn)
    end
    
    p0 = [k,c]
    init_vals = [[50.0, 10.0], [80.0, 20.0], [30.0, 5.0]]
    prbs = []
    total_time = 50
    for u in init_vals
        nichols = DeterministicIteratedMap(nicholson, u, p0)
        X, t = trajectory(nichols, total_time)
        push!(prbs, (X,t))
    end

    fig6 = Figure(size=(1000, 500))
    ax6 = Axis(fig6[1,1], title="Portretul de fază (Gazdă-Parazitoid)", 
               xlabel="Gazde (H)", ylabel="Paraziți (P)")
    colors = Makie.wong_colors()
    for (i, (X,t)) in enumerate(prbs)
        lines!(ax6, X[:,1], X[:,2], color=colors[i], linewidth=2)
        scatter!(ax6, [X[1,1]], [X[1,2]], markersize=12, color=colors[i], label="Start $i")
    end
    axislegend(ax6, position=:lt)

    ax61 = Axis(fig6[1,2], title="Seriile temporale", 
                xlabel="Generația (n)", ylabel="Populație")
    for (i, (X,t)) in enumerate(prbs)
        lines!(ax61, t, X[:,1], color=colors[i], linewidth=2, label="Gazde $i")
        lines!(ax61, t, X[:,2], color=colors[i], linewidth=2, linestyle=:dash, label="Paraziți $i")
    end
    axislegend(ax61, position=:rt)
    fig6
end

# ╔═╡ a1b2c3d4-0021-4a00-8000-000000000021
md"""
**Interpretare:** Graficul confirmă teoria: echilibrul netrivial este instabil. Traiectoriile nu converg spre coexistență, ci descriu oscilații din ce în ce mai ample (spirale care se îndepărtează de echilibru), ducând în final la extincția uneia sau a ambelor populații în modelul determinist.
"""

# ╔═╡ a1b2c3d4-0022-4a00-8000-000000000022
md"""
---
## 📝 5. Concluzii și Sinteză (10 min)

### 📊 Tabel Rezumtiv: Liniarizare și Stabilitate

| Concept | Definiție / Condiție | Implicație Biologică / Fizică |
| :--- | :--- | :--- |
| **Matricea Jacobiană** | Aproximarea liniară a unui sistem neliniar în jurul unui echilibru. | Permite folosirea algebrei liniare pentru a prezice soarta sistemului. |
| **Teorema Hartmann-Grobman** | Dacă $Re(\lambda) \neq 0$, sistemul neliniar se comportă local ca cel liniar. | Putem trage concluzii riguroase despre stabilitate fără a rezolva ecuațiile neliniare. |
| **Stabilitate Structurală** | Sistemul își păstrează calitativ portretul de fază la mici perturbații. | Esențială pentru modele biologice robuste (ex: necesită $Re(\lambda) \neq 0$). |
| **Criteriul Discret** | Echilibru stabil dacă toate $|\lambda| < 1$ pentru Jacobian. | Prezice coexistența sau oscilații divergente în modelele cu generații discrete. |

### 💡 Întrebări recapitulative
1. De ce teorema Hartmann-Grobman nu ne poate garanta comportamentul sistemului neliniar dacă o valoare proprie a matricei Jacobiene are partea reală zero?
2. Ce înseamnă ca un model matematic să fie "structural stabil" și de ce este această proprietate o cerință minimă atunci când încercăm să descriem un ecosistem real?
3. În modelul gazdă-parazitoid (Nicholson-Bailey), am demonstrat analitic și vizual că echilibrul este instabil. Ce mecanisme biologice (neliniare sau stocastice) ar putea explica totuși coexistența pe termen lung observată în natură între astfel de specii?

---
**✨ Sfârșitul Cursului 9 ✨**
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╠═a1b2c3d4-0002-4a00-8000-000000000002
# ╟─b574946d-dba3-43c2-9d90-011bc8e36f63
# ╟─8daab4b5-97af-41f2-b61c-32e6790bd036
# ╟─65f55b39-d3d6-49ae-9ce2-b541c90d064d
# ╟─183bb2b5-fb19-4887-a535-409b441f2936
# ╟─a48a49bf-ac5a-4361-9ded-8c1d056e957c
# ╠═c36059dc-75b1-4bbb-8ea9-0e38b6987d2a
# ╟─72e6323d-9ad8-4f0b-bb75-1f72fdc044e8
# ╟─a1b2c3d4-0003-4a00-8000-000000000003
# ╟─a1b2c3d4-0004-4a00-8000-000000000004
# ╟─a1b2c3d4-0005-4a00-8000-000000000005
# ╟─a1b2c3d4-0006-4a00-8000-000000000006
# ╠═a1b2c3d4-0007-4a00-8000-000000000007
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
