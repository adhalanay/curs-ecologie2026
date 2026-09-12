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
    using Distributions
    using StatsBase
    TableOfContents()
end

# ╔═╡ a1b2c3d4-0000-4a00-8000-000000000001
WidthOverDocs()

# ╔═╡ a1b2c3d4-0001-4a00-8000-000000000001
md"""
# 🎲 Cursul 12 — Variabile Aleatoare și Teoreme Limită

**Durata estimată:** 100 minute

## 🎯 Obiectivele lecției
1. Formalizarea conceptului de **variabilă aleatoare** (discretă și continuă).
2. Înțelegerea funcțiilor de probabilitate (PMF/PDF), a repartiției cumulative (CDF) și a momentelor (Media, Varianța).
3. Explorarea distribuțiilor fundamentale: **Geometrică, Exponențială, Erlang și Weibull** (esențiale în analiza de supraviețuire).
4. Aplicarea **Legii Numerelor Mari** și a **Teoremei Limită Centrale** în contexte medicale și biologice.
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000002
md"""
---
## 📐 1. Variabile Aleatoare: Formalizare (15 min)

Adesea, în experimentele probabiliste, suntem mai interesați de o proprietate numerică a rezultatului decât de rezultatul în sine. De exemplu, dacă aruncăm două monede, nu ne interesează secvența exactă (ex: $HT$), ci *numărul* de steme obținute.

O **variabilă aleatoare** $X$ este o funcție care asociază un număr real fiecărui rezultat posibil din spațiul de selecție $\Omega$:
$$X: \Omega \to \mathbb{R}$$

În funcție de valorile pe… o obținem în medie pe termen lung:
$$\mu = E[X] = \sum_{i} x_i \cdot p(x_i)$$

**Varianța** (notată cu $\sigma^2$ sau $Var[X]$) măsoară dispersia valorilor în jurul mediei:
$$\sigma^2 = Var[X] = E[(X - \mu)^2] = \sum_{i} (x_i - \mu)^2 \cdot p(x_i)$$
*Formulă alternativă (foarte utilă în calcule):* $Var[X] = E[X^2] - (E[X])^2 = \left(\sum_{i} x_i^2 p(x_i)\right) - \mu^2$.

**Deviația standard** este $\sigma = \sqrt{Var[X]}$ și se exprimă în aceleași unități de măsură ca datele originale.
"""

# ╔═╡ 185738ea-79f9-4f10-a48f-6bfc3eaff4e0
md"""
### 📏 1.2. Variabile Aleatoare Continue
O variabilă aleatoare este **continuă** dacă ia valori într-un interval din $\mathbb{R}$. Deoarece există o infinitate de valori posibile, probabilitatea ca variabila să ia o valoare *exactă* este zero: $P(X = x) = 0$. Ne interesează doar probabilitatea ca $X$ să cadă într-un interval.

Pentru a o descrie, folosim:
1. **Funcția de Densitate de Probabilitate (PDF - Probability Density Function)**: $f(x)$.
   *Proprietăți:* $f(x) \geq 0$, $\int_{-\inf…X \leq x) = \int_{-\infty}^x f(t) dt$. 
   *Observație:* Prin urmare, $f(x) = F'(x)$.

#### 🎯 Media și Varianța (Cazul Continuu)
Trecem de la sume la integrale pentru a calcula momentele distribuției:

**Media (Valoarea Așteptată):**
$$\mu = E[X] = \int_{-\infty}^{\infty} x \cdot f(x) dx$$

**Varianța:**
$$\sigma^2 = Var[X] = \int_{-\infty}^{\infty} (x - \mu)^2 \cdot f(x) dx$$

*Formulă alternativă:* $Var[X] = E[X^2] - (E[X])^2 = \left(\int_{-\infty}^{\infty} x^2 f(x) dx\right) - \mu^2$.
"""

# ╔═╡ 80578626-232e-47b8-82b2-2a74369a0a8a
md"""
### 📊 1.3. Rezumat: Discret vs. Continuu

Pentru a fixa conceptele, iată o comparație directă a "uneltelor" matematice pentru cele două tipuri de variabile:

| Concept | Variabilă Discretă (Numărare) | Variabilă Continuă (Măsurare) |
| :--- | :--- | :--- |
| **Funcția de probabilitate** | Masa (PMF): $p(x_i) = P(X=x_i)$ | Densitatea (PDF): $f(x)$ |
| **Aria totală** | Suma: $\sum p(x_i) = 1$ | Integrala: $\int f(x) dx = 1$ |
| **Probabilitate pe interval** | Sumă: $\sum_{a \le x_i \le b} p(x_i)$ | Integrală: $\int_a^b f(x) dx$ |
| **Media ($\mu$)** | $\sum x_i \cdot p(x_i)$ | $\int x \cdot f(x) dx$ |
| **Varianța ($\sigma^2$)** | $\sum (x_i - \mu)^2 \cdot p(x_i)$ | $\int (x - \mu)^2 \cdot f(x) dx$ |

Acum că avem formalismul matematic, vom explora în secțiunile următoare cele mai importante distribuții de probabilitate utilizate în biologie și medicină, începând cu cele **discrete**.
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000003
md"""
---
## 🔢 2. Distribuții Discrete Fundamentale (25 min)

Vom analiza câteva distribuții discrete esențiale, care modelează evenimente de numărare.
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000004
md"""
### 🪙 2.1. Distribuția Geometrică
Modelează **numărul de eșecuri** (sau încercările) până la **primul succes** într-o succesiune de experimente Bernoulli independente.
$$P(X=k) = (1-p)^k p, \quad k \in \{0, 1, 2, \dots\}$$
*Exemplu biologic:* Numărul de diviziuni celulare până când apare prima mutație de rezistență la un medicament.
"""

# ╔═╡ 2f4a67fa-38d2-4a3e-8593-d2946184e31f
md"""
🧫 **Exemplu Biologic: Rezistența la Antibiotice și PCR**
Distribuția geometrică modelează „norocul” până la primul succes.

1. **Oncologie / Bacteriologie:** Care este numărul de diviziuni celulare necesare până când apare *prima* mutație critică ce conferă rezistență la un antibiotic? Dacă probabilitatea de mutație per diviziune este $p$, numărul de diviziuni "eșuate" urmează o distribuție geometrică.

2. **Biologie moleculară (qPCR):** În reacția de PCR cantitativ, câte cicluri de amplificare sunt necesare pentru ca semnalul de fluorescență să depășească pragul de detecție (Ct - threshold cycle)? Această valoare este invers proporțională cu 
încărcătura virală inițială.
"""

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000005
md"""
Probabilitatea de succes: $p = {}$ $(@bind p_geo PlutoUI.Slider(0.05:0.05:0.95, default=0.2, show_value=true))
"""

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000006
let
    d1 = Geometric(p_geo)
    max_trials = min(50, quantile(d1, 0.999) + 1)
    trials = 0:max_trials
    pmf_trials = [pdf(d1, k) for k in trials]
    
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1,1], title="Distribuția Geometrică (p = $p_geo)",
              xlabel="Număr de eșecuri până la primul succes (k)", 
              ylabel="P(X = k)")
    
    barplot!(ax, trials, pmf_trials, color=:teal, strokecolor=:black)
    vlines!(ax, [mean(d1)], color=:red, linestyle=:dash, linewidth=2, label="Media = $(round(mean(d1), digits=2))")
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000007
md"""
**Interpretare:** Pe măsură ce $p$ scade (evenimentul este rar), "coada" distribuției se lungește dramatic, iar media $\frac{1-p}{p}$ crește.
"""

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000008
md"""
### 🦠 2.2. Distribuția Poisson
Modelează numărul de evenimente rare într-un interval fix. Este limita distribuției Binomiale când $n \to \infty$ și $p \to 0$.
$$P(X=k) = e^{-\lambda}\frac{\lambda^k}{k!}, \quad \lambda > 0$$

"""

# ╔═╡ 9eb62253-326e-414d-ab9f-4dfecacdf693
md"""
⏱️ **Exemplu Biologic: Timpul de așteptare și Radioactivitatea**

Dacă Poisson numără *evenimentele*, Exponențiala măsoară *timpul* dintre ele.

1. **Neurofiziologie:** Timpul scurs între două potențiale de acțiune succesive (intervalul inter-spike) într-un neuron supus unui stimul constant.

2. **Medicină Nucleară:** Timpul de înjumătățire. Timpul până când un atom radioactiv (ex: Technetium-99m, folosit în scintigrafii) se dezintegrează urmează o distribuție exponențială.

3. **Epidemiologie:** Timpul de așteptare între sosirea a doi pacienți consecutivi la camera de gardă a unui spital.
"""

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000009
md"""
Rata medie: $\lambda = {}$ $(@bind λ_pois PlutoUI.Slider(0.5:0.5:15.0, default=4.0, show_value=true))
"""

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000010
let
    dist = Poisson(λ_pois)
    k_max = min(50, ceil(Int, λ_pois + 4*sqrt(λ_pois)))
    k_vals = 0:k_max
    pmf = pdf.(dist, k_vals)
    
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1,1], title="Distribuția Poisson (λ = $λ_pois)",
              xlabel="Evenimente (k)", ylabel="P(X = k)")
    
    barplot!(ax, k_vals, pmf, color=:purple, strokecolor=:black)
    vlines!(ax, [λ_pois], color=:red, linestyle=:dash, linewidth=2, label="Media = $λ_pois")
    
    # Modulul este floor(λ) pentru Poisson
    mode_val = floor(Int, λ_pois)
    scatter!(ax, [mode_val], [pdf(dist, mode_val)], color=:black, markersize=10, label="Modulul = $mode_val")
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000011
md"""
---
## 📏 3. Variabile Aleatoare Continue (30 min)

O variabilă aleatoare continuă ia valori într-un interval din $\mathbb{R}$. Probabilitatea ca $X$ să ia o valoare exactă este 0; ne interesează probabilitatea ca $X$ să cadă într-un interval:
$$P(a \leq X \leq b) = \int_a^b f(x) dx$$
Aria totală de sub curba $f(x)$ este egală cu 1.
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000012
md"""
### ⏳ 3.1. Distribuția Exponențială și Erlang (Gamma)
Dacă distribuția Poisson numără *evenimentele* într-un timp fix, distribuția Exponențială modelează **timpul de așteptare** până la *primul* eveniment.
$$f(x) = \lambda e^{-\lambda x}, \quad x \geq 0$$
Media este $1/\lambda$.

Dacă așteptăm până la **al $n$-lea eveniment**, obținem **Distribuția Erlang** (un caz particular al distribuției Gamma):
$$f_{n,\lambda}(x) = \frac{\lambda^n x^{n-1}}{(n-1)!} e^{-\lambda x}, \quad x \geq 0$$
"""

# ╔═╡ 56053b8c-ab04-482a-b817-58cbf22f5e4e
md"""
🔄 **Exemplu Biologic: Acumularea de evenimente și Recidive**
Erlang este timpul de așteptare până la *al $n$-lea* eveniment.

1. **Oncologie:** Timpul necesar pentru ca o celulă sănătoasă să acumuleze exact $n=3$ mutații „driver” (critice) pentru a deveni o celulă canceroasă malignă. Fiecare mutație apare aleatoriu (Poisson), deci timpul total până la a 3-a mutație urmează o distribuție Erlang/Gamma.

2. **Psihiatrie / Boli cronice:** Timpul de supraviețuire liber de boală până la *al doilea* episod de recidivă la pacienții cu scleroză multiplă sau depresie majoră.
"""

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000013
let
    xs = LinRange(0, 20, 200)
    λ0 = 0.3
    n0 = 3 # Erlang cu n=3 (timpul până la al 3-lea eveniment)
    
    # Densitatea Exponențială
    fe(x) = λ0 * exp(-λ0 * x)
    # Densitatea Erlang (Gamma cu shape=n0)
    ferlang(x) = (λ0^n0 * x^(n0-1) / factorial(n0-1)) * exp(-λ0 * x)
    
    ys_e = [fe(x) for x in xs]
    ys_er = [ferlang(x) for x in xs]
    
    fig = Figure(size=(900, 400))
    ax = Axis(fig[1,1], title="Timpul de așteptare: Exponențial vs Erlang",
              xlabel="Timp (x)", ylabel="Densitate f(x)")
    
    lines!(ax, xs, ys_e, color=:orange, linewidth=3, label="Exponențială (1-ul eveniment)")
    lines!(ax, xs, ys_er, color=:blue, linewidth=3, label="Erlang (al $n0-lea eveniment)")
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000014
md"""
### 📉 3.2. Distribuția Weibull (Analiza de Supraviețuire)
Este "regele" distribuțiilor în ingineria fiabilității și analiza de supraviețuire medicală. Poate modela rate de deces constante, crescătoare sau descrescătoare în funcție de parametrul de formă $k$.
$$f_{\lambda,k}(x) = \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}\exp\left(-\left(\frac{x}{\lambda}\right)^k\right), \quad x \geq 0$$
"""

# ╔═╡ f6652180-0c3c-4f4a-8ebc-ce17382e26f1
md"""
📉 **Exemplu Biologic: Analiza de Supraviețuire (Survival Analysis)**
Weibull este standardul de aur în medicina clinică și ingineria biomedicală pentru că poate modela riscuri care cresc, scad sau rămân constante.
1. **Implantologie:** Timpul de „supraviețuire” (funcționare) al unei proteze de șold sau al unui stimulator cardiac până la cedarea mecanică/bateriei.
2. **Oncologie clinică:** Timpul de supraviețuire al pacienților diagnosticați cu un anumit stadiu de cancer, de la momentul operației până la deces sau recidivă. 
   - Dacă $k < 1$: Mortalitate infantilă / eșecuri timpurii ale unui implant.
   - Dacă $k = 1$: Risc constant (echivalent cu Exponențiala).
   - Dacă $k > 1$: Risc crescut odată cu îmbătrânirea / uzura (ex: cancer pulmonar, unde riscul de deces crește pe măsură ce boala avansează).
"""

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000015
md"""
Parametrul de formă: $k = {}$ $(@bind k_wei PlutoUI.Slider(0.5:0.1:5.0, default=1.5, show_value=true))
"""

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000016
md"""
Parametrul de scară: $\lambda = {}$ $(@bind λ_wei PlutoUI.Slider(0.5:0.1:5.0, default=2.0, show_value=true))
"""

# ╔═╡ a1b2c3d4-0017-4a00-8000-000000000017
let
    xs = LinRange(0.01, 15, 300)
    
    function fw(x, k, λ)
        (k / λ) * (x / λ)^(k - 1) * exp(-(x / λ)^k)
    end
    
    ys = [fw(x, k_wei, λ_wei) for x in xs]
    
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1,1], title="Distribuția Weibull (k=$k_wei, λ=$λ_wei)",
              xlabel="Timp / Vârstă (x)", ylabel="Densitate f(x)")
    
    lines!(ax, xs, ys, color=:red, linewidth=3)
    
    # Interpretare parametru k
    if k_wei < 1
        text!(ax, "Rata de deces scade în timp (mortalitate infantilă)", position=(8, maximum(ys)*0.8), align=(:left, :top), color=:gray)
    elseif k_wei ≈ 1
        text!(ax, "Rata de deces constantă (proces Poisson / Exponențial)", position=(8, maximum(ys)*0.8), align=(:left, :top), color=:gray)
    else
        text!(ax, "Rata de deces crește în timp (îmbătrânire / uzură)", position=(8, maximum(ys)*0.8), align=(:left, :top), color=:gray)
    end
    
    fig
end

# ╔═╡ a1b2c3d4-0018-4a00-8000-000000000018
md"""
### 🔔 3.3. Aproximarea Normală a Binomialei
Când $n$ este mare, distribuția Binomială devine "clopotul" lui Gauss. Aceasta este o consecință directă a Teoremei Limită Centrale.
"""

# ╔═╡ a1b2c3d4-0019-4a00-8000-000000000019
let
    n1 = 40.0
    p2 = 0.6
    m = n1 * p2
    s = sqrt(n1 * p2 * (1 - p2))
    d2 = Binomial(n1, p2)
    
    k_vals = 0:Int(n1)
    pmf_vals = pdf.(d2, k_vals)
    
    xs = LinRange(0, 40, 400)
    ys_norm = [1/(s*sqrt(2*pi)) * exp(-0.5*((x-m)/s)^2) for x in xs]
    
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], title="Aproximarea Normală a Binomialei (n=$n1, p=$p2)",
              xlabel="Succese (k)", ylabel="Probabilitate")
    
    barplot!(ax, k_vals, pmf_vals, color=(:blue, 0.5), strokecolor=:blue, label="Binomială (Discretă)")
    lines!(ax, xs, ys_norm, color=:red, linewidth=3, label="Normală (Continuă)")
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0020-4a00-8000-000000000020
md"""
---
## ⚖️ 4. Legi Limită în Statistică (20 min)

De ce sunt atât de importante mediile și distribuția Normală în studiile clinice? Răspunsul stă în două teoreme fundamentale.
"""

# ╔═╡ a1b2c3d4-0021-4a00-8000-000000000021
md"""
### 🏥 4.1. Legea Numerelor Mari (LLN)
> *Cu cât eșantionul este mai mare, cu atât estimarea este mai aproape de realitate.*

**Exemplu Medical (Antigenul PSA):**
Nivelul antigenului specific prostatei (PSA) este folosit pentru a detecta cancerul. Un nivel peste 0.5 ng/ml indică celule canceroase rămase după operație.
- **Studiu Mare:** 429 de pacienți urmăriți 5 ani $\rightarrow$ 8% au avut nivel crescut.
- **Studiu Mic:** 30 de pacienți $\rightarrow$ 3 au avut nivel crescut (10%).

Care rezultat este mai credibil? Evident primul. LLN ne garantează că media eșantionului $\overline{X}_n$ converge în probabilitate către media teoretică $\mu$ pe măsură ce $n \to \infty$.
"""

# ╔═╡ a1b2c3d4-0022-4a00-8000-000000000022
md"""
### 🌌 4.2. Teorema Limită Centrală (CLT)
Acesta este cel mai important rezultat din teoria probabilităților. Ne spune că, dacă însumăm un număr mare de variabile aleatoare independente și identic distribuite (i.i.d.), suma lor va tinde către o **distribuție Normală**, *indiferent de forma distribuției originale*.

Dacă $S_n = \sum_{i=1}^n X_i$, atunci variabila standardizată:
$$Z = \frac{S_n - n\mu}{\sqrt{n\sigma^2}}$$
urmează o distribuție Normală Standard $\mathcal{N}(0, 1)$ când $n$ este mare.

**Implicație majoră:** Chiar dacă o variabilă biologică (ex: concentrația unei enzime) este puternic asimetrică (Exponențială sau Weibull), **media** unui eșantion de 100 de pacienți va fi distribuită Normal. Acest lucru permite folosirea testelor parametrice (ex: testul t, ANOVA) în medicină.
"""

# ╔═╡ a1b2c3d4-0023-4a00-8000-000000000023
md"""
---
## 📝 5. Concluzii și Sinteză (10 min)

### 📊 Tabel Rezumtiv: Distribuții și Utilizări

| Distribuție | Tip | Parametri | Utilizare Biologică / Medicală |
| :--- | :--- | :--- | :--- |
| **Geometrică** | Discretă | $p$ (prob. succes) | Timpul (nr. de generații) până la prima mutație de rezistență. |
| **Poisson** | Discretă | $\lambda$ (rata) | Numărul de colonii de bacterii pe o placă Petri. |
| **Exponențială** | Continuă | $\lambda$ (rata) | Timpul de așteptare între două sosiri de pacienți la UPU. |
| **Erlang / Gamma** | Continuă | $n, \lambda$ | Timpul de așteptare până la al $n$-lea eveniment (ex: recidivă). |
| **Weibull** | Continuă | $k$ (formă), $\lambda$ (scară) | Analiza de supraviețuire; modelarea timpului până la deces. |

### 💡 Întrebări de reflecție
1. Dacă distribuția Weibull are parametrul de formă $k < 1$, ce ne spune acest lucru despre riscul de deces al unui pacient pe măsură ce trece timpul? (Gândiți-vă la nou-născuți vs. vârstnici).
2. De ce Distribuția Erlang este considerată o "punte" între procesul Poisson (numărare) și distribuțiile de timp de așteptare (continue)?
3. În exemplul cu antigenul PSA, de ce un studiu pe 30 de pacienți care raportează 10% rată de recidivă nu este suficient pentru a contrazice studiul pe 429 de pacienți care raportează 8%? Cum intervine Legea Numerelor Mari aici?

---
**✨ Sfârșitul Cursului 12 ✨**
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╠═a1b2c3d4-0002-4a00-8000-000000000002
# ╟─185738ea-79f9-4f10-a48f-6bfc3eaff4e0
# ╟─80578626-232e-47b8-82b2-2a74369a0a8a
# ╟─a1b2c3d4-0003-4a00-8000-000000000003
# ╟─a1b2c3d4-0004-4a00-8000-000000000004
# ╟─2f4a67fa-38d2-4a3e-8593-d2946184e31f
# ╟─a1b2c3d4-0005-4a00-8000-000000000005
# ╠═a1b2c3d4-0006-4a00-8000-000000000006
# ╟─a1b2c3d4-0007-4a00-8000-000000000007
# ╠═a1b2c3d4-0008-4a00-8000-000000000008
# ╠═9eb62253-326e-414d-ab9f-4dfecacdf693
# ╟─a1b2c3d4-0009-4a00-8000-000000000009
# ╠═a1b2c3d4-0010-4a00-8000-000000000010
# ╟─a1b2c3d4-0011-4a00-8000-000000000011
# ╟─a1b2c3d4-0012-4a00-8000-000000000012
# ╟─56053b8c-ab04-482a-b817-58cbf22f5e4e
# ╠═a1b2c3d4-0013-4a00-8000-000000000013
# ╟─a1b2c3d4-0014-4a00-8000-000000000014
# ╟─f6652180-0c3c-4f4a-8ebc-ce17382e26f1
# ╟─a1b2c3d4-0015-4a00-8000-000000000015
# ╟─a1b2c3d4-0016-4a00-8000-000000000016
# ╠═a1b2c3d4-0017-4a00-8000-000000000017
# ╟─a1b2c3d4-0018-4a00-8000-000000000018
# ╠═a1b2c3d4-0019-4a00-8000-000000000019
# ╟─a1b2c3d4-0020-4a00-8000-000000000020
# ╟─a1b2c3d4-0021-4a00-8000-000000000021
# ╟─a1b2c3d4-0022-4a00-8000-000000000022
# ╟─a1b2c3d4-0023-4a00-8000-000000000023
