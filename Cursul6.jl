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
    TableOfContents()
end

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 📈 Cursul 6 — Comportament Periodic și Bifurcații în Sistemele Dinamice

**Durata estimată:** 100 minute
"""

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000002
md"""
## 🎯 Obiectivele lecției
1. Înțelegerea conceptelor de **atractor** și **ciclu limită**.
2. Explorarea modelelor biologice care prezintă oscilații stabile (clarinet, hormoni, genetică).
3. Analiza **bifurcațiilor locale**: transcritică, șa-nod și Hopf.
4. Interpretarea tranzițiilor calitative ale sistemelor dinamice în funcție de parametri.
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000002
md"""
## 🌀 1. Comportament Periodic și Atractori (30 min)

Până acum am studiat sisteme care tind către puncte de echilibru. În natură, însă, multe fenomene prezintă **oscilații persistente**: bătăile inimii, ciclurile circadiene, oscilațiile glicemiei sau interacțiunea prădător-pradă. 
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000003
md"""
### 🌌 1.0. Ce este un Atractor?

În spațiul fazelor, un **atractor** este o mulțime de stări către care un sistem tinde să evolueze, indiferent de condițiile inițiale (atâta timp cât acestea sunt într-o vecinătate numită *bazin de atracție*). 

Există trei tipuri principale de atractori în sistemele dinamice continue:

| Tipul Atractorului | Descriere | Exemplu Biologic / Fizic |
| :--- | :--- | :--- |
| **📍 Punct Fix (Echilibru)** | Sistemul evoluează către o stare constantă și se oprește. | Populație care atinge capacitatea de suport a mediului. |
| **🔄 Ciclu Limită** | Sistemul oscilează periodic; traiectoria în spațiul fazelor este o curbă închisă izolată. | Bătăile inimii, oscilații glicolitice, clarinet. |
| **🌪️ Atractor Ciudat (Haotic)** | Traiectorii care nu se repetă niciodată, dar rămân mărginite. Extrem de sensibil la condițiile inițiale. | Dinamica haotică a unor ecosisteme, turbulențe. |

În acest curs, ne vom concentra în special pe **ciclurile limită** și pe modul în care ele apar sau dispar prin **bifurcații**.
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000003
md"""
### 🎷 1.1. Clarinetul lui Rayleigh

Un prim exemplu clasic de ciclu limită este modelul clarinetului lui Rayleigh. Acesta descrie o oscilație auto-întreținută, unde energia este injectată la amplitudini mici și disipată la amplitudini mari.

Ecuația este:
$$\begin{aligned} x' &= v \\ v' &= -x - (v^3 - v) \end{aligned}$$

Termenul neliniar $-(v^3 - v)$ acționează ca o **frecare negativă** pentru $|v| < 1$ (amplifică oscilațiile) și ca o **frecare pozitivă** pentru $|v| > 1$ (le atenuează).
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000004
let
    function f_rayleigh!(du, u, p, t)
        du[1] = u[2]
        du[2] = -u[1] - (u[2]^3 - u[2])
    end
    
    tspan = (0.0, 30.0)
    colors = Makie.wong_colors()
    inits = ([0.1, 0.1], [2.0, 2.0], [0.5, -1.5])
    
    fig = Figure(size=(1000, 450))
    
    # Phase portrait
    ax1 = Axis(fig[1, 1], title="Portretul de fază", xlabel="x", ylabel="v")
    xs = LinRange(-2.5, 2.5, 20)
    ys = LinRange(-2.5, 2.5, 20)
    f_vec(x) = Point2f(x[2], -x[1] - (x[2]^3 - x[2]))
    streamplot!(ax1, f_vec, xs, ys, colormap=:magma, linewidth=1, arrow_size=5)
    scatter!(ax1, [0], [0], markersize=15, color=:red, label="Echilibru (0,0)")
    axislegend(ax1, position=:rb)

    # Time series
    ax2 = Axis(fig[1, 2], title="Serii temporale", xlabel="Timp (t)", ylabel="x(t), v(t)")
    
    for (i, u0) in enumerate(inits)
        prob = ODEProblem(f_rayleigh!, u0, tspan, ())
        sol = solve(prob, Tsit5(), saveat=0.1)
        lines!(ax2, sol.t, [u[1] for u in sol.u], color=colors[i], linewidth=2, label="x₀=$(u0[1])")
    end
    axislegend(ax2, position=:rt)
    
    fig
end

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000005
md"""
**Observație:** Indiferent de condiția inițială (mică sau mare), traiectoriile converg către aceeași orbită periodică (ciclul limită). Aceasta este semnătura unui atractor de tip ciclu limită.
"""

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000006
md"""
### 🧬 1.2. Reglarea Hormonală (Feedback Negativ cu Întârziere)

Multe oscilații biologice apar din mecanisme de **feedback negativ** combinate cu **întârzieri** (delay). Să considerăm un hormon $G$ produs de o gonadă, controlat de hormonul hipofizar $P$, care la rândul lui este controlat de hormonul hipotalamic $H$.

$$\begin{aligned} H' &= \frac{1}{1+G^n} - k_1 H \\ P' &= H - k_2 P \\ G' &= P - k_3 G \end{aligned}$$

Când $G$ este mare, producția de $H$ este inhibată. Parametrul $n$ controlează "puterea" sau "cooperativitatea" feedback-ului.
"""

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000007
md"""
Cooperativitatea feedback-ului: $n = {}$ $(@bind n_horm PlutoUI.Slider(2:10, default=8, show_value=true))
"""

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000008
let
    function f_horm!(du, u, p, t)
        n_val, k1, k2, k3 = p
        du[1] = 1.0 / (u[3]^n_val + 1.0) - k1 * u[1]
        du[2] = u[1] - k2 * u[2]
        du[3] = u[2] - k3 * u[3]
    end
    
    p_horm = [n_horm, 0.2, 0.2, 0.2]
    u0_horm = [0.2, 0.2, 0.2]
    tspan = (0.0, 100.0)
    prob = ODEProblem(f_horm!, u0_horm, tspan, p_horm)
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    fig = Figure(size=(1000, 450))
    
    ax1 = Axis(fig[1, 1], title="Serii temporale (n = $n_horm)", xlabel="Timp", ylabel="Concentrație")
    colors = Makie.wong_colors()
    labels = ["H (Hipotalamus)", "P (Hipofiză)", "G (Gonadă)"]
    for i in 1:3
        lines!(ax1, sol.t, [u[i] for u in sol.u], color=colors[i], linewidth=2, label=labels[i])
    end
    axislegend(ax1, position=:rt)
    
    ax2 = Axis3(fig[1, 2], title="Spațiul fazelor 3D", xlabel="H", ylabel="P", zlabel="G")
    lines!(ax2, [u[1] for u in sol.u], [u[2] for u in sol.u], [u[3] for u in sol.u], linewidth=2, color=:teal)
    scatter!(ax2, [u0_horm[1]], [u0_horm[2]], [u0_horm[3]], markersize=15, color=:red, label="Start")
    
    fig
end

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000009
md"""
**Interpretare:** Pentru $n < 8$, sistemul converge la un echilibru stabil (punct fix). Pentru $n \ge 8$, feedback-ul este suficient de puternic și neliniar încât echilibrul devine instabil, apărând un **ciclu limită**. Întârzierea introdusă de cele 3 ecuații este esențială pentru apariția oscilațiilor.
"""

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000010
md"""
### 🧫 1.3. Expresii Genetice Oscilante (Hes1)

Factorul de transcripție `Hes1` are un comportament oscilant cu o perioadă de aproximativ 2 ore, esențial pentru dezvoltarea embrionară. Modelul matematic include ARN-ul mesager ($Y$), proteina ($X$) și un factor de degradare ($Z$).

$$\begin{aligned} X' &= -AXZ + BY - CX \\ Y' &= \frac{E}{1+X^2} - DY \\ Z' &= -AXZ + \frac{F}{1+X^2} - GZ \end{aligned}$$
"""

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000011
let
    function hes1!(du, u, p, t)
        A, B, C, D, E, F, G = p
        du[1] = -A * u[1] * u[3] + B * u[2] - C * u[1]
        du[2] = E / (1 + u[1]^2) - D * u[2]
        du[3] = -A * u[1] * u[3] + F / (1 + u[1]^2) - G * u[3]
    end
    
    p_hes = [0.022, 0.3, 0.031, 0.028, 0.5, 20.0, 0.3]
    u0_hes = [0.2, 0.2, 0.2]
    tspan = (0.0, 150.0)
    prob = ODEProblem(hes1!, u0_hes, tspan, p_hes)
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    fig = Figure(size=(1000, 450))
    ax1 = Axis(fig[1, 1], title="Modelul Hes1 - Serii temporale", xlabel="Timp", ylabel="Concentrație")
    colors = Makie.wong_colors()
    labels = ["Proteina (X)", "ARNm (Y)", "Degradare (Z)"]
    for i in 1:3
        lines!(ax1, sol.t, [u[i] for u in sol.u], color=colors[i], linewidth=2, label=labels[i])
    end
    axislegend(ax1, position=:rt)
    
    ax2 = Axis3(fig[1, 2], title="Atractorul 3D", xlabel="X", ylabel="Y", zlabel="Z")
    lines!(ax2, [u[1] for u in sol.u], [u[2] for u in sol.u], [u[3] for u in sol.u], linewidth=2, color=:purple)
    
    fig
end

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000012
md"""
---
## 🦋 2. Bifurcații Locale (50 min)

Am văzut că modificarea unui parametru poate schimba complet comportamentul unui sistem (de la echilibru stabil la oscilații). O **bifurcație** este o modificare calitativă a structurii spațiului fazelor (numărul sau stabilitatea echilibrelor) atunci când un parametru trece printr-o valoare critică.
"""

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000013
md"""
### 🔄 2.1. Bifurcația Transcritică (Efectul Allee)

Să reconsiderăm modelul populației cu efect Allee:
$$x' = 0.1 x \left(1-\frac{x}{k}\right)\left(\frac{x}{a}-1\right)$$
Echilibrele sunt $x^* \in \{0, a, k\}$. Când parametrul $a$ (pragul Allee) trece prin $k$ (capacitatea de suport), două echilibre se ciocnesc și **își schimbă stabilitatea**. Aceasta este o bifurcație transcritică.
"""

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000014
md"""
Pragul Allee: $a = {}$ $(@bind a_alle PlutoUI.Slider(300:50:1500, default=800, show_value=true))
"""

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000015
let
    function alle!(du, u, p, t)
        a_val, k_val = p
        du[1] = 0.1 * u[1] * (1 - u[1] / k_val) * (u[1] / a_val - 1)
    end
    
    k0 = 800.0
    p_alle = [a_alle, k0]
    tspan = (0.0, 100.0)
    
    fig = Figure(size=(1000, 450))
    ax1 = Axis(fig[1, 1], title="Efectul Allee (a = $a_alle, k = $k0)", xlabel="Timp", ylabel="Populație x(t)")
    ax2 = Axis(fig[1, 2], title="Diagrama Echilibrelor", xlabel="x", ylabel="Rata de creștere f(x)")
    
    colors = Makie.wong_colors()
    inits = [[100.0], [400.0], [600.0], [1000.0], [1200.0]]
    
    for (i, u0) in enumerate(inits)
        prob = ODEProblem(alle!, u0, tspan, p_alle)
        sol = solve(prob, Tsit5(), saveat=0.1)
        lines!(ax1, sol.t, [u[1] for u in sol.u], color=colors[mod1(i, length(colors))], linewidth=2, label="x₀=$(Int(u0[1]))")
    end
    axislegend(ax1, position=:rb)
    
    # Plot f(x) for equilibrium analysis
    xs = LinRange(0, 1500, 500)
    ys = [0.1 * x * (1 - x / k0) * (x / a_alle - 1) for x in xs]
    lines!(ax2, xs, ys, linewidth=3, color=:blue)
    hlines!(ax2, [0], linestyle=:dash, color=:black, linewidth=1)
    
    # Mark equilibria
    eq_points = [0.0, a_alle, k0]
    scatter!(ax2, eq_points, zeros(length(eq_points)), markersize=12, color=:red, label="Echilibre")
    axislegend(ax2, position=:lb)
    
    fig
end

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000016
md"""
**Interpretare:** În panoul din dreapta, observăm cum echilibrele $a$ și $k$ se apropie. Când $a = k$, ele fuzionează (bifurcația), iar pentru $a > k$, $k$ devine instabil, iar $a$ devine stabil. Populația va tinde spre $a$ în loc de $k$.
"""

# ╔═╡ a1b2c3d4-0017-4a00-8000-000000000017
md"""
### 🐎 2.2. Bifurcația Șa-Nod (Lac Operon)

Modelul Lac Operon este:
$$x' = \frac{a_0 + x^2}{1 + x^2} - r x$$
Aici, $a_0$ este rata bazală, iar $r$ este rata de degradare. Echilibrele sunt intersecțiile dintre curba de producție (Hill) și dreapta de degradare.
"""

# ╔═╡ a1b2c3d4-0018-4a00-8000-000000000018
md"""
Rata de degradare: $r = {}$ $(@bind r_lac PlutoUI.Slider(0.1:0.05:1.5, default=0.5, show_value=true))
"""

# ╔═╡ a1b2c3d4-0019-4a00-8000-000000000019
let
    function lac!(du, u, p, t)
        a0, r_val = p
        du[1] = (a0 + u[1]^2) / (1 + u[1]^2) - r_val * u[1]
    end
    
    a0 = 0.1
    p_lac = [a0, r_lac]
    tspan = (0.0, 50.0)
    
    fig = Figure(size=(1000, 450))
    ax1 = Axis(fig[1, 1], title="Lac Operon - Serii temporale (r = $r_lac)", xlabel="Timp", ylabel="Concentrație x(t)")
    ax2 = Axis(fig[1, 2], title="Analiza Grafică a Echilibrelor", xlabel="x", ylabel="Rata de variație")
    
    colors = Makie.wong_colors()
    inits = [[0.1], [0.5], [1.0], [2.0]]
    
    for (i, u0) in enumerate(inits)
        prob = ODEProblem(lac!, u0, tspan, p_lac)
        sol = solve(prob, Tsit5(), saveat=0.1)
        lines!(ax1, sol.t, [u[1] for u in sol.u], color=colors[i], linewidth=2, label="x₀=$(u0[1])")
    end
    axislegend(ax1, position=:rt)
    
    xs = LinRange(0.0, 3.0, 500)
    prod = [(a0 + x^2) / (1 + x^2) for x in xs]
    deg = [r_lac * x for x in xs]
    
    lines!(ax2, xs, prod, linewidth=3, color=:blue, label="Producție: (a₀+x²)/(1+x²)")
    lines!(ax2, xs, deg, linewidth=3, color=:red, linestyle=:dash, label="Degradare: r·x")
    axislegend(ax2, position=:rt)
    
    fig
end

# ╔═╡ a1b2c3d4-0020-4a00-8000-000000000020
md"""
**Interpretare:** Pentru $r$ mare, există un singur echilibru stabil. Pe măsură ce $r$ scade, dreapta roșie devine mai puțin abruptă și apar **două noi echilibre** (unul stabil, unul instabil). Această "naștere" a unei perechi de echilibre este **bifurcația șa-nod**.
"""

# ╔═╡ a1b2c3d4-0021-4a00-8000-000000000021
md"""
### 💫 2.3. Bifurcația Hopf

Bifurcația Hopf este mecanismul principal prin care un punct de echilibru stabil se transformă într-un **ciclu limită**. Aceasta apare când valorile proprii ale matricei Jacobiene devin complexe conjugate cu partea reală care trece de la negativ la pozitiv.
"""

# ╔═╡ a1b2c3d4-0022-4a00-8000-000000000022
md"""
#### 🍷 2.3.1. Modelul Glicolizei (Sel'kov)

Glicoliza este procesul de conversie a glucozei în energie. Enzima PFK este activată de propriul său produs (ADP), creând un feedback pozitiv care, combinat cu degradarea, duce la oscilații.

$$\begin{aligned} S' &= 1 - c S P^2 \\ P' &= c S P^2 - k P \end{aligned}$$
unde $S$ este substratul (F6P) și $P$ este produsul (ADP).
"""

# ╔═╡ a1b2c3d4-0023-4a00-8000-000000000023
md"""
Rata de conversie: $c = {}$ $(@bind c_glico PlutoUI.Slider(0.1:0.05:2.0, default=0.5, show_value=true))
"""

# ╔═╡ a1b2c3d4-0024-4a00-8000-000000000024
let
    function glico!(du, u, p, t)
        c_val, k_val = p
        du[1] = 1.0 - c_val * u[1] * u[2]^2
        du[2] = c_val * u[1] * u[2]^2 - k_val * u[2]
    end
    
    p_glico = [c_glico, 1.0]
    u0_glico = [1.2, 1.2]
    tspan = (0.0, 100.0)
    prob = ODEProblem(glico!, u0_glico, tspan, p_glico)
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    fig = Figure(size=(1000, 450))
    ax1 = Axis(fig[1, 1], title="Portret de fază (c = $c_glico)", xlabel="S (Substrat)", ylabel="P (Produs)")
    lines!(ax1, [u[1] for u in sol.u], [u[2] for u in sol.u], linewidth=2, color=:teal)
    scatter!(ax1, [u0_glico[1]], [u0_glico[2]], markersize=12, color=:red, label="Start")
    axislegend(ax1, position=:rt)
    
    ax2 = Axis(fig[1, 2], title="Serii temporale", xlabel="Timp", ylabel="Concentrație")
    lines!(ax2, sol.t, [u[1] for u in sol.u], linewidth=2, color=:blue, label="S")
    lines!(ax2, sol.t, [u[2] for u in sol.u], linewidth=2, color=:orange, label="P")
    axislegend(ax2, position=:rt)
    
    fig
end

# ╔═╡ a1b2c3d4-0025-4a00-8000-000000000025
md"""
**Interpretare:** Pentru $c$ mare, sistemul converge la un echilibru stabil. Scăzând $c$ (de exemplu, sub 0.4), echilibrul devine instabil și apare un ciclu limită. Aceasta este **bifurcația Hopf**.
"""

# ╔═╡ a1b2c3d4-0026-4a00-8000-000000000026
md"""
#### 🐺 2.3.2. Modelul Prădător-Pradă Holling-Tanner

Un model ecologic mai realist decât Lotka-Volterra, care include capacitatea de suport a prăzii și un răspuns funcțional de tip II (saturație).

$$\begin{aligned} N' &= r_1 N \left(1-\frac{N}{k}\right) - \frac{w N}{N+d} P \\ P' &= r_2 P \left(1-\frac{j P}{N}\right) \end{aligned}$$
"""

# ╔═╡ a1b2c3d4-0027-4a00-8000-000000000027
md"""
Capacitatea de consum: $w = {}$ $(@bind w_ht PlutoUI.Slider(0.1:0.05:2.0, default=0.5, show_value=true))
"""

# ╔═╡ a1b2c3d4-0028-4a00-8000-000000000028
let
    function holling!(du, u, p, t)
        r1, r2, k, d, j, w_val = p
        du[1] = r1 * u[1] * (1 - u[1] / k) - w_val * u[1] * u[2] / (u[1] + d)
        du[2] = r2 * u[2] * (1 - j * u[2] / u[1])
    end
    
    p_ht = [1.0, 0.1, 7.0, 1.0, 1.0, w_ht]
    u0_ht = [2.0, 2.0]
    tspan = (0.0, 200.0)
    prob = ODEProblem(holling!, u0_ht, tspan, p_ht)
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    fig = Figure(size=(1000, 450))
    ax1 = Axis(fig[1, 1], title="Portret de fază (w = $w_ht)", xlabel="Pradă (N)", ylabel="Prădător (P)")
    lines!(ax1, [u[1] for u in sol.u], [u[2] for u in sol.u], linewidth=2, color=:teal)
    scatter!(ax1, [u0_ht[1]], [u0_ht[2]], markersize=12, color=:red, label="Start")
    axislegend(ax1, position=:rt)
    
    ax2 = Axis(fig[1, 2], title="Serii temporale", xlabel="Timp", ylabel="Densitate")
    lines!(ax2, sol.t, [u[1] for u in sol.u], linewidth=2, color=:blue, label="Pradă (N)")
    lines!(ax2, sol.t, [u[2] for u in sol.u], linewidth=2, color=:orange, label="Prădător (P)")
    axislegend(ax2, position=:rt)
    
    fig
end

# ╔═╡ a1b2c3d4-0029-4a00-8000-000000000029
md"""
**Interpretare:** Parametrul $w$ reprezintă capacitatea maximă de consum a prădătorului. Creșterea lui $w$ duce la o bifurcație Hopf: de la coexistență pașnică (echilibru stabil) la cicluri de "boom și bust" (ciclu limită), un fenomen des observat în ecologie (paradoxul îmbogățirii).
"""

# ╔═╡ a1b2c3d4-0030-4a00-8000-000000000030
md"""
---
## 📝 3. Concluzii și Sinteză (20 min)

### Rezumat al conceptelor cheie

| Concept | Definiție | Exemplu din lecție |
| :--- | :--- | :--- |
| **🌌 Atractor** | Mulțime de stări către care sistemul tinde (punct fix, ciclu limită). | Toate modelele studiate. |
| **🔄 Ciclu Limită** | Traiectorie închisă izolată care atrage alte traiectorii. | Clarinetul Rayleigh, Glicoliza |
| **🔄 Bifurcație Transcritică** | Două echilibre se ciocnesc și schimbă stabilitatea. | Efectul Allee ($a$ și $k$) |
| **🐎 Bifurcație Șa-Nod** | Apariția sau dispariția unei perechi de echilibre. | Lac Operon (variația lui $r$) |
| **💫 Bifurcație Hopf** | Un echilibru stabil devine instabil și naște un ciclu limită. | Glicoliza, Holling-Tanner |

### Întrebări de reflecție
1. De ce este necesară o "întârziere" (delay) sau cel puțin 3 variabile de stare pentru a obține oscilații stable prin feedback negativ?
2. În modelul Holling-Tanner, cum explici din punct de vedere ecologic "paradoxul îmbogățirii" (creșterea lui $w$ duce la oscilații mai ample, risc de extincție)?
3. Ce înseamnă din punct de vedere calitativ o **bifurcație** și cum ne ajută studiul ei să înțelegem tranzițiile bruște din sistemele biologice (de exemplu, trecerea de la coexistență pașnică la oscilații periculoase în modelul Holling-Tanner)?

---
**✨ Sfârșitul Cursului 6 ✨**
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╠═a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000002
# ╟─a1b2c3d4-0002-4a00-8000-000000000002
# ╟─a1b2c3d4-0002-4a00-8000-000000000003
# ╟─a1b2c3d4-0003-4a00-8000-000000000003
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
# ╟─a1b2c3d4-0016-4a00-8000-000000000016
# ╟─a1b2c3d4-0017-4a00-8000-000000000017
# ╟─a1b2c3d4-0018-4a00-8000-000000000018
# ╟─a1b2c3d4-0019-4a00-8000-000000000019
# ╟─a1b2c3d4-0020-4a00-8000-000000000020
# ╟─a1b2c3d4-0021-4a00-8000-000000000021
# ╟─a1b2c3d4-0022-4a00-8000-000000000022
# ╟─a1b2c3d4-0023-4a00-8000-000000000023
# ╟─a1b2c3d4-0024-4a00-8000-000000000024
# ╟─a1b2c3d4-0025-4a00-8000-000000000025
# ╟─a1b2c3d4-0026-4a00-8000-000000000026
# ╟─a1b2c3d4-0027-4a00-8000-000000000027
# ╟─a1b2c3d4-0028-4a00-8000-000000000028
# ╟─a1b2c3d4-0029-4a00-8000-000000000029
# ╠═a1b2c3d4-0030-4a00-8000-000000000030
