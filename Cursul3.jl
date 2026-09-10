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

# ╔═╡ 96931d62-a744-11f0-296e-dd13032fbfff
begin
    using Pkg
    Pkg.activate(".")
    using Makie, CairoMakie
    using PlutoUI
    using DifferentialEquations
    using OrdinaryDiffEq
    using OrdinaryDiffEqSDIRK
    using OrdinaryDiffEqFIRK
    using QuadGK
    using Printf
    TableOfContents()
end

# ╔═╡ 34ee0b2b-12f1-4d5e-9739-be3e8f9de99b
html"""
<style>
    @media screen {
        main {
            margin: 0 auto;
            max-width: 2000px;
        padding-left: max(283px, 10%);
        padding-right: max(383px, 10%);
        }
    }
</style>
"""

# ╔═╡ 4ed380d4-314a-426f-9ae3-7597edc47808
md"""
# Cursul III — Integrala și Ecuații Diferențiale

## 🎯 Obiectivele lecției:
1. Să înțelegem **integrarea** ca inversă a derivării și ca arie sub grafic.
2. Să aplicăm **metoda Euler (sume Riemann)** pentru a aproxima primitive și a rezolva numeric ecuații diferențiale.
3. Să rezolvăm **ecuații diferențiale** simple: creștere și descompunere exponențială.
4. Să studiem **punctele de echilibru** și **stabilitatea** lor (metoda punctului intermediar și analiza liniară).
5. Să aplicăm aceste concepte la modelarea populațiilor **(Efectul Allee)**.

---
"""

# ╔═╡ 1a2b3c4d-0001-4000-8000-000000000001
md"""
## ⨌ Partea I: Integrarea și Metoda Euler


### 1.1. De la derivată la integrală
După cum am văzut în cursul precedent, este destul de ușor să calculăm derivata unei funcții pornind de la câteva formule și aplicând reguli de calcul. Dar putem inversa procesul? Adică să determinăm ``f(x)`` dacă știm doar ``f'(x)``? Procedura care ne permite să facem aceasta se numește **integrare** (sau aflarea primitivei).

Pentru funcții simple nu este foarte greu: dacă ``f'(x) = 2x``, atunci orice funcție de forma ``f(x) = x^2 + c`` are derivata `2x` pentru orice constantă `c`. De asemenea, dacă ``f'(x) = e^x``, atunci ``f(x) = eˣ + c``. 

În general însă, este mult mai dificil să calculăm primitiva decât derivata. Adesea cunoaștem derivata funcției fără să știm funcția însăși: avem acces la înregistrările vitezometrului unei mașini (viteza `x'(t)`) și vrem să determinăm distanța parcursă `x(t)`.
"""

# ╔═╡ 95bf4b62-d4c3-4714-bcb2-c108cfafa53b
begin
    f(x) = x - cos(x)*sin(x)
    t1 = range(0, 10, length=15)
    y1 = [f(t) for t ∈ t1]
    fig1 = Figure(size = (800, 400))
    ax1 = Axis(fig1[1,1],
    title = "Exemplu de grafic obținut experimental (date discrete)",
    xlabel = "t", ylabel = "f(t)")
    lines!(t1, y1, color = :steelblue, linewidth = 2, label = "date")
    scatter!(t1, y1, color = :black, markersize = 8, label = "măsurători")
    axislegend(ax1, position = :lt)
    fig1
end

# ╔═╡ 2a0ac124-34d7-455c-bccc-1f875083981a
md"""
### 1.2. Metoda Euler și Sumele Riemann
Atunci când primitiva nu se poate exprima prin funcții elementare, avem nevoie de metode de aproximare. Ne amintim de definiția derivatei:

```math 
X'(t) \approx \frac{X(t+\Delta t) - X(t)}{\Delta t} 
```

De aici, putem exprima valoarea viitoare:
```math 
X(t+\Delta t) \approx X(t) + \Delta t \cdot X'(t) 
```

Aceasta este **Metoda Euler**. Repetăm construcția pas cu pas pentru $ t = 0, \Delta t, 2\Delta t, \dots $:

```math
\begin{align}
    X(\Delta t)      &= X(0) + \Delta t \cdot X'(0)  \\ 
    X(2\Delta t)     &= X(\Delta t) + \Delta t \cdot X'(\Delta t) \\
    X(3\Delta t)     &= X(0) + X'(0)\Delta t + X'(\Delta t)\Delta t + X'(2\Delta t)\Delta t
\end{align}
```math

Deci, cu $n$ pași, poziția la timpul $t = n \Delta t$ este:
```math 
X(t) \approx X(0) + \sum_{k=0}^{n-1} X'(k \Delta t) \cdot \Delta t 
```

Trecând la limită cu $\Delta t \to 0$, obținem **Teorema Fundamentală a Analizei (Leibniz–Newton)**:
```math
X(t) = X(0) + \int_0^t X'(s) ds 
```

**Exemplu practic.** Viteza unei mașini este $V(t) = 2\sqrt{1000 - t^3}$. Mașina pornește cu $\approx 63$ m/s ($226$ km/h) și se oprește după $10$ s. Primitiva lui $V$ nu se exprimă prin funcții elementare, deci vom calcula distanța folosind sume Riemann (Euler).
"""

# ╔═╡ 34cc7eec-c950-4491-bc7f-8ec29ffef73e
begin
    V(t) = 2*sqrt(1000 - t^3)
    t2 = range(0, 10, length=200)
    y2 = V.(t2)
    integral, err = quadgk(V, 0, 10)
    fig2 = Figure(size = (800, 420))
    ax2 = Axis(fig2[1,1],
    title = @sprintf("V(t) = 2√(1000 − t³)   —   ∫₀¹⁰ V dt ≈ %.2f m", integral),
    xlabel = "t  [s]", ylabel = "V(t)  [m/s]")
    band!(t2, zero.(t2), y2, color = (:steelblue, 0.25))
    lines!(t2, y2, color = :steelblue, linewidth = 2)
    hlines!(0, color = :black, linewidth = 0.8)
    fig2
end

# ╔═╡ 5300c0d9-36b7-421f-ba05-af81b6c786ed
md"""
Presupunem că $X(0)=0$ și $\Delta t=0.1$. Vom compara suma Riemann cu valoarea exactă obținută numeric (`quadgk`). Tabelul de mai jos afișează doar primii și ultimii pași.
"""

# ╔═╡ c6560f36-88ef-4160-b63a-12a96a0d726c
begin
    Δt = 0.1
    N  = 100
    Xs = zeros(N + 1)
    for k in 0:N-1
        Xs[k+2] = Xs[k+1] + V(k*Δt) * Δt
    end
    t_euler = (0:N) .* Δt
    println("------------------------------------------------------------")
    @printf("%6s | %12s | %14s\n", "t", "V(t)", "X(t)")
    println("------------------------------------------------------------")
    for k in vcat(0:3, N-2:N)
        @printf("%6.1f | %12.4f | %14.4f\n", k*Δt, V(k*Δt), Xs[k+1])
    end
    println("------------------------------------------------------------")
    @printf("Sumă Riemann (Δt = %.1f)      : %.4f\n", Δt, Xs[end])
    @printf("Valoare exactă (quadgk)      : %.4f\n", integral)
    @printf("Eroare absolută              : %.4f\n", abs(Xs[end] - integral))
end

# ╔═╡ 4124bec3-78dd-4f8c-b3a6-ab7046bdbd19
begin
    fig3 = Figure(size = (820, 420))
    ax3 = Axis(fig3[1,1],
    title = "Deplasarea X(t) — sumă Riemann vs. integrală exactă",
    xlabel = "t  [s]", ylabel = "X(t)  [m]")
    lines!(t_euler, Xs, linewidth = 3, color = :crimson, label = "Euler (Δt = 0.1)")
    lines!(t2, [quadgk(V, 0, t)[1] for t in t2],
    linestyle = :dash, linewidth = 2, color = :black, label = "exact (quadgk)")
    axislegend(ax3, position = :lt)
    fig3
end

# ╔═╡ f1728f88-5926-4ccb-b1b2-0424d0c3027e
md"""
**Ce am făcut de fapt (Geometria integralei):**
1. Am descompus intervalul $(0, t)$ în segmente egale $\Delta t$;
2. Am presupus viteza $V$ constantă pe fiecare subinterval;
3. Am adunat ariile dreptunghiurilor $V(k \Delta t) \cdot \Delta t$.

Integrala este exact **aria de sub grafic**:

```math
\int_a^b f(t) dt = F(b) - F(a), \quad \text{unde } F' = f 
```

**Formula de medie:** Valoarea medie a lui $f$ pe $[a,b]$ este:
```math
\bar{f} = \frac{1}{b-a} \int_a^b f(t) dt 
```
"""

# ╔═╡ 2b3c4d5e-0002-4000-8000-000000000002
md"""
---
## 🌀 Partea a II-a: Ecuații Diferențiale și Modele Populaționale


### 2.1. Modele Exponențiale
Ecuațiile diferențiale ordinare (EDO) de forma $ x' = f(x) $ descriu cum evoluează un sistem în timp. Cele mai simple cazuri sunt $ x' = kx $ și $ x' = -kx $ cu $ k > 0 $.

**Modelul exponențial:** Dacă $r = b + i - d - e$ este rata netă de creștere (nașteri + imigrare − decese − emigrare), populația evoluează după:

```math
X' = r X \implies X(t) = X(0) \cdot e^{rt} 
```
Curba pentru $r < 0$ se numește *scădere exponențială* (ex: dezintegrare radioactivă), iar pentru $r > 0$, *creștere exponențială*.
"""

# ╔═╡ ae8a653a-48fa-4247-859b-d2f4dd195966
begin
    f_decay(t) = 0.3exp(-0.5t)
    y_decay = f_decay.(t2)
    fig4 = Figure(size = (800, 400))
    ax4 = Axis(fig4[1,1],
    title = "Scădere exponențială: X(t) = 0.3·e^(−0.5t)",
    xlabel = "t", ylabel = "X(t)")
    lines!(t2, y_decay, linewidth = 3, color = :navy)
    hlines!(0, color = :gray, linestyle = :dot)
    fig4
end

# ╔═╡ 18b60c22-1365-4e86-90c5-466eef2c856d
begin
    f_growth(t) = 0.3exp(0.5t)
    y_growth = f_growth.(t2)
    fig5 = Figure(size = (800, 400))
    ax5 = Axis(fig5[1,1],
    title = "Creștere exponențială: X(t) = 0.3·e^(+0.5t)",
    xlabel = "t", ylabel = "X(t)")
    lines!(t2, y_growth, linewidth = 3, color = :darkgreen)
    fig5
end

# ╔═╡ 894ad5c0-b051-45c7-b465-98a14659dd9f
md"""
### 2.2. Modelul Logistic și Punctele de Echilibru
Modelarea exponențială este valabilă doar pe termen scurt (prezice o creștere nesfârșită, ceea ce e nerealist biologic). Trebuie să introducem **limitarea creșterii** din cauza resurselor finite, ceea ce duce la **ecuația logistică**.

Pentru $X' = f(X)$, **punctele de echilibru** sunt stările staționare, adică soluțiile ecuației $f(X_0) = 0$ (populația nu se mai schimbă).

Să studiem un model cu efect de prag (Allee):
```math
x' = r \left(1 - \frac{x}{k}\right) \left(\frac{x}{a} - 1\right) 
```
cu $r = 0.1$, $a = 5$ (pragul minim), $k = 100$ (capacitatea mediului) și diverse valori inițiale.
"""

# ╔═╡ 608f4fa5-37e6-465e-886b-5412fa1b3e5c
begin
    r, a, k = 0.1, 5.0, 100.0
    f(u, p, t) = r*(1 - u/k)*(u/a - 1)
    tspan = (0.0, 800.0)
    solve_ic(u0) = solve(ODEProblem(f, u0, tspan), RadauIIA5(), isoutofdomain=(u,p,t) -> u[1] < 0,
    reltol=1e-6, abstol=1e-6, saveat=0.01)
    sol1 = solve_ic(10.0)
    sol2 = solve_ic(20.0)
    sol3 = solve_ic(150.0)
    fig6 = Figure(size = (820, 440))
    ax6 = Axis(fig6[1,1],
        title = "x' = r(1 − x/k)(x/a − 1),   r = 0.1, a = 5, k = 100",
        xlabel = "t", ylabel = "x(t)")
    hlines!([a, k], color = :gray, linestyle = :dash, linewidth = 1)
    text!(0.0, k, text = " x = k", color = :gray, offset = (10, 5))
    text!(0.0, a, text = " x = a", color = :gray, offset = (10, 5))
    lines!(sol1.t, reduce(vcat, sol1.u), linewidth = 3, label = "x₀ = 10")
    lines!(sol2.t, reduce(vcat, sol2.u), linewidth = 3, color = :crimson, label = "x₀ = 20")
    lines!(sol3.t, reduce(vcat, sol3.u), linewidth = 3, color = :seagreen, label = "x₀ = 150")
    axislegend(ax6, position = :rc)
    fig6
end

# ╔═╡ df77753d-59db-4378-97ad-33551d0973de
begin
    sol4 = solve_ic(4.9)
    fig7 = Figure(size = (820, 440))
    ax7 = Axis(fig7[1,1],
    title = "Caz sub-prag: x₀ = 4.9 < a — soluția scade la 0 (extincție)",
    xlabel = "t", ylabel = "x(t)")
    hlines!([a, k], color = :gray, linestyle = :dash)
    lines!(sol4.t, reduce(vcat, sol4.u), linewidth = 4, color = :darkorange,
    label = "x₀ = 4.9")
    axislegend(ax7, position = :rt)
    fig7
end

# ╔═╡ e25589cd-579f-47f5-be8c-7573e76cc2de
@bind u1 html"<input type='range' style='width: 1000px' min='0' max='200' step='0.05'>"

# ╔═╡ 1630ec2f-6aa8-4e8f-80cd-4bcfe6026bbf
md"**``u_1`` = $(round(u1, digits=2))**"

# ╔═╡ 916800aa-7986-48ad-af69-934e115e5972
begin
    sol_var = solve_ic(u1)
    fig_var = Figure(size = (820, 440))
    ax_var = Axis(fig_var[1,1],
    title = @sprintf("Soluție pentru x₀ = %.2f", u1),
    xlabel = "t", ylabel = "x(t)")
    hlines!([a, k], color = :gray, linestyle = :dash, linewidth = 1)
    text!(0.0, k, text = " x = k", color = :gray, offset = (10, 5))
    text!(0.0, a, text = " x = a", color = :gray, offset = (10, 5))
    lines!(sol_var.t, reduce(vcat, sol_var.u), linewidth = 4, color = :purple)
    fig_var
end

# ╔═╡ f171e206-6e46-4525-9ecb-cf4b3a982b63
md"""
Dacă pornim exact dintr-un echilibru, soluția rămâne constantă. Să testăm $x_0 = k$ și $x_0 = a$.
"""

# ╔═╡ b0eacf7f-e6a9-4f98-9287-a58d87881397
begin
    sol5 = solve_ic(a)
    sol6 = solve_ic(k)
    fig8 = Figure(size = (820, 440))
    ax8 = Axis(fig8[1,1],
    title = "Soluții pornind exact din punctele de echilibru",
    xlabel = "t", ylabel = "x(t)")
    lines!(sol5.t, reduce(vcat, sol5.u), linewidth = 4,
    color = :steelblue, label = "x₀ = a = 5")
    lines!(sol6.t, reduce(vcat, sol6.u), linewidth = 4,
    color = :crimson, label = "x₀ = k = 100")
    axislegend(ax8, position = :rt)
    fig8
end

# ╔═╡ 7a5c562f-8c8e-41a1-88b2-dfaf2153bbdc
md"""
### 2.3. Stabilitatea Echilibrelor
Se observă comportamentul **bistabil**: 
- dacă $x_0 > a$, soluția tinde la $k$ (supraviețuire); 
- dacă $x_0 < a$, soluția tinde la $0$ (extincție). 

Punctele de echilibru sunt $x = 0, a, k$. Dar cum știm *fără să simulăm* dacă un echilibru este stabil (atracție) sau instabil (respingere)?

**Metoda 1: Punctul intermediar (Calitativ)**
Idee: Semnul lui $f(x)$ de o parte și de alta a echilibrului ne spune dacă soluțiile se apropie sau se îndepărtează. Fiind o funcție continuă, $f$ nu poate schimba semnul fără să treacă prin zero.
- Dacă $f(x) > 0$ la stânga și $f(x) < 0$ la dreapta $\implies x$ crește spre echilibru și scade spre el $\implies$ **Stabil**.
- Dacă $f(x) < 0$ la stânga și $f(x) > 0$ la dreapta $\implies x$ se îndepărtează $\implies$ **Instabil**.
"""

# ╔═╡ 6db1dbae-20f2-44c5-9a91-4ee14cdf744f
begin
    fig9 = Figure(size = (900, 220))
    ax9 = Axis(fig9[1,1], title = "Linia de fază: x = 0 instabil, x = k stabil (model logistic simplu)",
    xlabel = "x", yticks = [])
    hlines!(0, color = :black, linewidth = 1)
    for (x, lbl) in [(0.0, "x = 0"), (10.0, "x = k")]
        scatter!(ax9, x, 0, color = :black, markersize = 12)
        text!(ax9, x, 0, text = lbl, align = (:center, :bottom), offset = (0, 12))
    end
    arrows2d!(ax9, [2.0], [0.0], [3.0], [0.0],
    color = :crimson, tailwidth = 3, taillength= 14)
    arrows2d!(ax9, [8.0], [0.0], [-3.0], [0.0],
    color = :crimson, tailwidth = 3, taillength = 14)
    xlims!(ax9, -2, 14)
    fig9
end

# ╔═╡ 3ea69d1b-21ee-4da3-b10f-2cd5e4e858cd
md"""
**Exemplu cu pescuit.** Modelul logistic cu recoltare constantă:
```math
x' = 0.2 x (1 - x/1000) - 0.1 x
```
Echilibrele se găsesc rezolvând $x(0.2(1 - x/1000) - 0.1) = 0 \implies x = 0$ și $x = 500$. 
Semnul lui $x'$ pe intervalele $(0, 500)$ și $(500, \infty)$ se determină luând un punct-test (ex: $x=100 \implies x'>0$, deci populația crește spre 500).
"""

# ╔═╡ a5291a15-baac-4a5c-9334-e90dcd8468e7
md"""
**Metoda 2: Analiza Liniară (Criteriul derivatei)**
Pentru ecuația $x' = f(x)$, calculăm derivata $f'(x)$ în punctul de echilibru $x^*$.

*De ce funcționează intuitiv?* Dacă $f'(x^*) < 0$, funcția $f$ este descrescătoare în $x^*$. Asta înseamnă că dacă populația $x$ crește puțin peste $x^*$, $f(x)$ devine negativ, deci $x$ va scădea înapoi. Dacă $x$ scade sub $x^*$, $f(x)$ devine pozitiv, împingând $x$ din nou spre echilibru.

**Criteriu liniar (Hartman–Grobman, 1D):**
- `` f'(x^*) < 0 \implies x^* `` **stabil** (atractor)
- ``f'(x^*) > 0 \implies x^*`` **instabil** (respingător)

*Exemplu:* Pentru $x' = x(1 - x/k)$, avem $f(x) = x - x^2/k$, deci $f'(x) = 1 - 2x/k$.
- ``f'(0) = 1 > 0 \implies 0`` este instabil.
- ``f'(k) = -1 < 0 \implies k`` este stabil.
"""

# ╔═╡ 9654f5b8-c980-48ca-8dc9-99da68f331a9
begin
    k1 = 5
    f3(x) = x*(1 - x/k1)
    xs = range(0.0, 7.0, length = 200)
    fig10 = Figure(size = (900, 420))
    ax10 = Axis(fig10[1,1],
    title = "f(x) = x(1 − x/k),  k = 5",
    xlabel = "x", ylabel = "f(x)")
    hlines!(0, color = :black, linewidth = 1)
    lines!(xs, f3.(xs), linewidth = 3, color = :steelblue)
    for x in (0.0, Float64(k1))
        scatter!(ax10, x, 0, color = :black, markersize = 12)
    end
    text!(0.0, 0, text = "x = 0", align = (:center, :bottom), offset = (0, 10))
    text!(Float64(k1), 0, text = "x = k", align = (:center, :bottom), offset = (0, 10))
    xlims!(ax10, -0.8, 7.0)
    ylims!(ax10, -3.5, 3.5)
    fig10
end

# ╔═╡ c5bd9a17-9455-4a64-a7e6-d7ce4dae29b6
md"""
### 2.4. Efectul Allee
În anumite populații (ex: balene, plante care au nevoie de polenizare încrucișată) există un prag minim $a$ sub care populația nu se mai poate reproduce eficient și dispare. 

Modelul matematic:
```math
x' = r x \left(1 - \frac{x}{k}\right) \left(\frac{x}{a} - 1\right) 
```

Echilibre: $x = 0$, $x = a$, $x = k$. Să verificăm stabilitatea cu analiza liniară:
```math 
f'(0) = -r < 0 \quad \text{(stabil)} 
```

```math 
f'(a) =  r \left(1 - \frac{a}{k}\right) > 0 \quad \text{(instabil - prag critic)}
```

```math
 f'(k) =  r \left(1 - \frac{k}{a}\right) < 0 \quad \text{(stabil)} 
```
"""

# ╔═╡ d8f4e2a1-9b3c-4d5e-8f7a-1c2d3e4f5a6b
begin
    k2, a2, r2 = 7.0, 3.0, 5.0
    f4(x) = r2 * x * (1 - x/k2)*(x/a2 - 1)
    xs4 = range(0.0, 9.0, length = 300)
    fig11 = Figure(size = (900, 440))
    ax11 = Axis(fig11[1,1],
    title =  "Efectul Allee: f(x) = r x(1 − x/k)(x/a − 1) ",
    xlabel =  "x ", ylabel =  "f(x) ")
    hlines!(0, color = :black, linewidth = 1)
    for x in (0.0, a2, k2)
        vlines!(x, color = :gray, linestyle = :dash, linewidth = 1)
        scatter!(ax11, x, 0, color = :black, markersize = 12 )
    end
    text!(0.0, 0, text =  " x = 0 (stabil) ",   align = (:left,   :top))
    text!(a2,  0, text =  " x = a (instabil) ", align = (:center, :bottom), offset = (0, 8))
    text!(k2,  0, text =  " x = k (stabil) ",   align = (:center, :bottom), offset = (0, 8))
    lines!(xs4, f4.(xs4), linewidth = 3, color = :darkorange)
    xlims!(ax11, -0.5, 9.0)
    ylims!(ax11, -6.0, 8.0)
    fig11
end

# ╔═╡ e9f5d3b2-0c4d-5e6f-9a8b-2d3e4f5a6b7c
md"""
---
## Partea a III-a: Concluzii și Exerciții


### 📌 Concluzii ale lecției
1. **Integrarea** este operația inversă derivării și are o interpretare geometrică clară: **aria de sub curbă**. Când nu avem primitive analitice, folosim **Sume Riemann / Metoda Euler**.
2. **Ecuațiile diferențiale** ($x' = f(x)$) modelează rata de schimbare. Soluția lor ne arată evoluția în timp a sistemului.
3. **Punctele de echilibru** ($f(x^*) = 0$) reprezintă stările staționare (ex: capacitatea mediului, extincția).
4. **Stabilitatea** ne spune dacă sistemul revine la echilibru după o mică perturbație. O putem determina calitativ (linia de fază) sau cantitativ (semnul derivatei $f'(x^*)$).
5. **Modelele realiste** (precum cel Logistic sau Allee) introduc neliniarități care explică fenomene complexe precum extincția sub un anumit prag.

---

### 📝 Exerciții pentru lucrul individual
*Rezolvați următoarele probleme în celulele de cod de mai jos.*

**Exercițiul 1 (Integrare Numerică)**
Funcția $g(x) = e^{-x^2}$ nu are o primitivă elementară. 

a) Folosiți metoda Euler / Sume Riemann pentru a aproxima $\int_0^2 e^{-x^2} dx$ cu $N=10, 100, 1000$ pași.

b) Comparați rezultatul cu valoarea exactă obținută cu `quadgk`.
"""

# ╔═╡ f0a6e4c3-1d5e-6f7a-0b9c-3e4f5a6b7c8d
# EXERCIȚIUL 1: Completați codul aici
# g(x) = ...
# N_values = [10, 100, 1000]
# ...

# ╔═╡ a1b7f5d4-2e6f-7a8b-1c0d-4f5a6b7c8d9e
md"""
**Exercițiul 2 (Modelare Exponențială)**
Un izotop radioactiv are un timp de înjumătățire de $5$ ani. 

a) Scrieți ecuația diferențială $N' = -k N$ și găsiți constanta $k$.

b) Dacă avem inițial $100$ g de substanță, cât timp durează până când rămân doar $10$ g?
"""

# ╔═╡ b2c8a6e5-3f7a-8b9c-2d1e-5a6b7c8d9e0f
md"""
**Exercițiul 3 (Stabilitate și Analiză Liniară)**
Fie sistemul dinamic descris de ecuația:
```math
x' = x^3 - 4x
```
a) Găsiți toate punctele de echilibru.

b) Folosiți analiza liniară (calculând $f'(x)$) pentru a determina stabilitatea fiecărui punct de echilibru.

c) Desenați schema liniei de fază (ca în Figura 9) pe o foaie de hârtie.
"""

# ╔═╡ c3d9b7f6-4a8b-9c0d-3e2f-6b7c8d9e0f1a
# EXERCIȚIUL 3: Puteți folosi Julia pentru a verifica derivatele și echilibrele
# f_ex3(x) = ...
# ...

# ╔═╡ d4e0c8a7-5b9c-0d1e-4f3a-7b8c9d0e1f2a
md"""
**Exercițiul 4 (Logistic cu Recoltare - Provocare)**
Modelul pescuitului este $x' = 0.5 x (1 - x/100) - H$, unde $H$ este cantitatea recoltată constant.

a) Găsiți punctele de echilibru în funcție de $H$.

b) Ce se întâmplă cu populația dacă $H = 10$? Dar dacă $H = 15$?

c) Care este valoarea maximă a lui $H$ pentru care populația nu se stinge (fenomen cunoscut sub numele de *bifurcație*)?
"""

# ╔═╡ d95cbef1-b38e-4848-9bff-5e2192024ec4


# ╔═╡ 7aba95d1-dd96-49be-8a9e-16e990758291


# ╔═╡ Cell order:
# ╠═96931d62-a744-11f0-296e-dd13032fbfff
# ╟─34ee0b2b-12f1-4d5e-9739-be3e8f9de99b
# ╟─4ed380d4-314a-426f-9ae3-7597edc47808
# ╟─1a2b3c4d-0001-4000-8000-000000000001
# ╟─95bf4b62-d4c3-4714-bcb2-c108cfafa53b
# ╟─2a0ac124-34d7-455c-bccc-1f875083981a
# ╟─34cc7eec-c950-4491-bc7f-8ec29ffef73e
# ╟─5300c0d9-36b7-421f-ba05-af81b6c786ed
# ╟─c6560f36-88ef-4160-b63a-12a96a0d726c
# ╠═4124bec3-78dd-4f8c-b3a6-ab7046bdbd19
# ╟─f1728f88-5926-4ccb-b1b2-0424d0c3027e
# ╟─2b3c4d5e-0002-4000-8000-000000000002
# ╟─ae8a653a-48fa-4247-859b-d2f4dd195966
# ╟─18b60c22-1365-4e86-90c5-466eef2c856d
# ╟─894ad5c0-b051-45c7-b465-98a14659dd9f
# ╠═608f4fa5-37e6-465e-886b-5412fa1b3e5c
# ╠═df77753d-59db-4378-97ad-33551d0973de
# ╟─e25589cd-579f-47f5-be8c-7573e76cc2de
# ╟─1630ec2f-6aa8-4e8f-80cd-4bcfe6026bbf
# ╟─916800aa-7986-48ad-af69-934e115e5972
# ╟─f171e206-6e46-4525-9ecb-cf4b3a982b63
# ╟─b0eacf7f-e6a9-4f98-9287-a58d87881397
# ╟─7a5c562f-8c8e-41a1-88b2-dfaf2153bbdc
# ╟─6db1dbae-20f2-44c5-9a91-4ee14cdf744f
# ╟─3ea69d1b-21ee-4da3-b10f-2cd5e4e858cd
# ╟─a5291a15-baac-4a5c-9334-e90dcd8468e7
# ╟─9654f5b8-c980-48ca-8dc9-99da68f331a9
# ╟─c5bd9a17-9455-4a64-a7e6-d7ce4dae29b6
# ╟─d8f4e2a1-9b3c-4d5e-8f7a-1c2d3e4f5a6b
# ╟─e9f5d3b2-0c4d-5e6f-9a8b-2d3e4f5a6b7c
# ╟─f0a6e4c3-1d5e-6f7a-0b9c-3e4f5a6b7c8d
# ╟─a1b7f5d4-2e6f-7a8b-1c0d-4f5a6b7c8d9e
# ╟─b2c8a6e5-3f7a-8b9c-2d1e-5a6b7c8d9e0f
# ╟─c3d9b7f6-4a8b-9c0d-3e2f-6b7c8d9e0f1a
# ╟─d4e0c8a7-5b9c-0d1e-4f3a-7b8c9d0e1f2a
# ╠═d95cbef1-b38e-4848-9bff-5e2192024ec4
# ╟─7aba95d1-dd96-49be-8a9e-16e990758291
