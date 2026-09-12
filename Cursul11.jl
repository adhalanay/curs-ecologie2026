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
# 📊 Cursul 11 — Statistică, Probabilități și Distribuții

**Durata estimată:** 100 minute

## 🎯 Obiectivele lecției
1. Extragerea de informații din date brute folosind **statistica descriptivă**.
2. Analiza **relațiilor** între variabile categoriale și numerice.
3. Înțelegerea conceptelor de **combinatorică** și **probabilitate condiționată** în contexte biologice.
4. Modelarea fenomenelor aleatoare folosind **distribuții clasice** (Binomială, Poisson, Normală) și **legile limită**.
"""

# ╔═╡ a1b2c3d4-0002-4a00-8000-000000000002
md"""
---
## 📏 1. Statistică Descriptivă (20 min)

Statistica descriptivă ne permite să sintetizăm un set mare de date măsurate pe un eșantion de indivizi. Variabilele măsurate pot fi **categoriale** (ex: genotip, prezența unei boli) sau **numerice** (ex: concentrația unui hormon, numărul de paraziți).

Pentru variabilele numerice, căutăm două tipuri de informații:
1. **Tendința centrală**: Unde se concentrează datele? (Media, Mediana, Modulul)
2. **Dispersia**: Cât de mult se împrăștie datele față de centru? (Varianța, Deviația standard, IQR)
"""

# ╔═╡ 28496985-d587-466f-ba69-8db8a67635e9
md"""
### 🧬 1.1. Tipuri de Variabile Statistice

Înainte de a calcula orice indicator, trebuie să clasificăm datele pe care le colectăm. Obiectele studiate se numesc *indivizi* (ex: un pacient, o libelulă, o celulă), iar proprietățile măsurate se numesc **variabile**.

Variabilele se împart în două mari categorii:

**1. Variabile Categoriale (Calitative)**
Descriu apartenența la o anumită clasă sau categorie. Nu au sens operațiile aritmetice (nu putem face "media" a două grupe de sânge).
*   **… sens operațiile aritmetice.
*   **Discrete:** Rezultă din *numărare*. Pot lua doar valori întregi (adesea din $\mathbb{N}$).
    *   *Exemple biologice:* Numărul de paraziți pe o gazdă, numărul de mutații într-o secvență de ADN, numărul de celule albe dintr-un volum de sânge.
*   **Continue:** Rezultă din *măsurători*. Pot lua orice valoare reală într-un interval.
    *   *Exemple biologice:* Concentrația de glucoză din sânge, tensiunea arterială, înălțimea, rata de creștere a unei colonii.
"""

# ╔═╡ 0684c3db-fdd9-4c91-9801-f3bac7bfa73a
md"""
### 📐 1.2. Descriptori Numerici: Tendința Centrală și Dispersia

Pentru **variabilele numerice**, scopul statisticii descriptive este să rezume un tabel mare de date (un eșantion de $n$ valori: $x_1, x_2, \dots, x_n$) folosind câțiva indicatori cheie.

#### 🎯 A. Măsuri ale Tendinței Centrale (Unde se află "centrul" datelor?)
1.  **Media ($\mu$ sau $\overline{x}$):** Suma tuturor valorilor împărțită la numărul lor.
    $$\overline{x} = \frac{x_1 + x_2 + \dots + x_n}{n}$$
    *Atenție:* …le.
3.  **Cuartile și IQR (Interquartile Range):**
    *   $Q_1$ (Sfertul 1): 25% din date sunt mai mici decât această valoare.
    *   $Q_2$ (Sfertul 2): Este chiar **Mediana** (50%).
    *   $Q_3$ (Sfertul 3): 75% din date sunt mai mici decât această valoare.
    *   **$IQR = Q_3 - Q_1$**: Măsoară întinderea celor 50% din mijlocul datelor. Este folosit pentru a identifica matematic **outliers** (orice valoare mai mică decât $Q_1 - 1.5 \times IQR$ sau mai mare decât $Q_3 + 1.5 \times IQR$).
"""

# ╔═╡ a1b2c3d4-0003-4a00-8000-000000000003
md"""
### 🧫 1.3. Exemplu Biologic: Paraziți pe libelule

Să aplicăm acești descriptori pe un set de date reale. Variabila noastră este **numerică discretă** (număr de paraziți). 
Să considerăm numărul de paraziți pe un eșantion de 20 de libelule:
`{5, 0, 0, 54, 5, 12, 27, 24, 36, 5, 56, 43, 15, 42, 12, 62, 36, 34, 58, 23}`

Putem calcula acești indicatori direct în Julia folosind pachetul `StatsBase`:
"""

# ╔═╡ a1b2c3d4-0004-4a00-8000-000000000004
let
    lst = [5,0,0,54,5,12,27,24,36,5,56,43,15,42,12,62,36,34,58,23]
    
    println("📊 Sumar statistic pentru datele despre paraziți:")
    println("--------------------------------------------------")
    println("Media (μ):            ", round(mean(lst), digits=2))
    println("Mediana:              ", median(lst))
    println("Modulul:              ", mode(lst))
    println("Deviația standard:    ", round(std(lst), digits=2))
    println("Sfertul 1 (Q1):       ", percentile(lst, 25))
    println("Sfertul 3 (Q3):       ", percentile(lst, 75))
    println("IQR (Q3 - Q1):        ", iqr(lst))
end

# ╔═╡ a1b2c3d4-0005-4a00-8000-000000000005
md"""
**Interpretare:** Observăm că media (27.45) este mai mare decât mediana (25.5), ceea ce sugerează o **asimetrie pozitivă** (datele sunt trase la dreapta de câteva valori foarte mari). Aceste valori extreme se numesc **outliers** și pot fi identificate folosind regula $1.5 \times IQR$.
"""

# ╔═╡ a1b2c3d4-0006-4a00-8000-000000000006
let
    lst = [5,0,0,54,5,12,27,24,36,5,56,43,15,42,12,62,36,34,58,23]
    
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], title="Histogramă: Numărul de paraziți pe libelule",
              xlabel="Număr de paraziți", ylabel="Frecvență (Probabilitate)")
    
    hist!(ax, lst, bins=10, color=:lightblue, strokecolor=:black, normalization=:probability)
    fig
end

# ╔═╡ a1b2c3d4-0007-4a00-8000-000000000007
md"""
### 📊 1.2. Vizualizarea Datelor Categoriale
Pentru variabilele categoriale (ex: localizarea unei boli), folosim grafice de bare sau diagrame circulare.
"""

# ╔═╡ a1b2c3d4-0008-4a00-8000-000000000008
let
    locatii = ["Col uterin", "Gură", "Organe genitale"]
    counts = [500000, 10000, 50000]
    colors = Makie.wong_colors()
    
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], title="Localizarea cancerelor induse de HPV în SUA",
              xlabel="Localizare", ylabel="Număr de cazuri")
    
    x_pos = 1:length(locatii)
    barplot!(ax, x_pos, counts, color=colors[1:3], width=0.6)
    ax.xticks = (x_pos, locatii)
    fig
end

# ╔═╡ a1b2c3d4-0009-4a00-8000-000000000009
let
    sizes = [43, 349, 299, 37, 33]
    labels = ["Canada", "China", "Hong Kong", "Taiwan", "Singapore"]
    colors = Makie.wong_colors()
    
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1,1], title="Distribuția globală a deceselor SARS (2002-2003)", aspect=DataAspect())
    
    pie!(ax, sizes, color=colors[1:5], strokecolor=:black, strokewidth=1)
    fig
end

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000010
md"""
---
## 🔗 2. Relații între Variabile (15 min)

Cum determinăm dacă două variabile măsurate pe același eșantion sunt corelate sau dacă aparțin unor categorii dependente?
"""

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000011
md"""
### 🦠 2.1. Variabile Categoriale: Tabele de Contingență
Pentru a studia efectele antiviralului Rimantadină au fost luați 76 de copii împărțiți în două grupuri. Tabelul de contingență rezultată este:

|  | Infectați | Neinfectați | Total |
| :--- | :--- | :--- | :--- |
| **Placebo** | 20 | 21 | 41 |
| **Rimantadină** | 1 | 34 | 35 |

Tabelul sugerează puternic că infecția se corelează cu lipsa medicamentului.
"""

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000012
let
    categorii = ["Placebo", "Rimantadină"]
    infectati = [20, 1]
    neinfectati = [21, 34]
    
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], title="Efectul Rimantadinei asupra Infecției",
              xlabel="Tratament", ylabel="Număr de pacienți")
    
    x_pos = 1:2
    width = 0.35
    
    barplot!(ax, x_pos .- width/2, infectati, width=width, color=:red, label="Infectați")
    barplot!(ax, x_pos .+ width/2, neinfectati, width=width, color=:green, label="Neinfectați")
    
    ax.xticks = (x_pos, categorii)
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0010-4a00-8000-000000000013
md"""
### 📈 2.2. Variabile Numerice: Scatter Plot și Regresia Liniară
Un studiu în Marea Britanie a măsurat consumul mediu de țigări și probabilitatea de deces prin cancer pulmonar pentru 25 de grupuri ocupaționale.
"""

# ╔═╡ a1b2c3d4-0011-4a00-8000-000000000014
let
    x = [77,102,115,137,91,105,117,104,87,94,107,91,116,112,100,102,113,76,111,110,66,93,125,88,133]
    y = [84,88,128,116,104,115,123,129,79,128,86,85,155,96,120,101,144,60,118,139,51,113,113,104,146]
    
    mx, my = mean(x), mean(y)
    a = (mean(x .* y) - mx * my) / (mean(x .^ 2) - mx^2)
    b = my - a * mx
    y_fit(x_val) = a * x_val + b
    
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1,1], title="Corelația Fumat - Cancer Pulmonar",
              xlabel="Consum mediu de țigări", ylabel="Rata mortalității (indice)")
    
    scatter!(ax, x, y, color=:red, markersize=8, label="Date observate")
    lines!(ax, range(minimum(x)-5, maximum(x)+5, length=100), 
           y_fit.(range(minimum(x)-5, maximum(x)+5, length=100)), 
           color=:blue, linewidth=3, linestyle=:dash, label="Regresie liniară")
    axislegend(ax, position=:lt)
    fig
end

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000015
md"""
---
## 🔢 3. Recapitulare Combinatorică (10 min)

Pentru a calcula probabilități, trebuie să știm să numărăm rezultatele posibile.
- **Permutări**: Aranjări ordonate. $A_n^k = \frac{n!}{(n-k)!}$
- **Combinări**: Submulțimi neordonate. $C_n^k = \frac{n!}{k!(n-k)!}$

### 🧬 Exemplu Biologic: Diversitatea HIV
După infectarea cu HIV, apar numeroase genotipuri. La transmitere, doar câteva se transmit. Să presupunem că o persoană infectată posedă **150** de genotipuri, dar transmite doar **10**. Câte posibilități pentru transmiterea genotipurilor există?
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000016
let
    N = 150
    k = 10
    nr_combinatii = binomial(big(N), big(k))
    println("Numărul de combinații posibile (C_150^10) este:")
    println(nr_combinatii)
end

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000017
md"""
---
## 🎲 4. Elemente de Teoria Probabilităților (20 min)

### 📜 4.1. Axiome și Probabilități Condiționate
O probabilitate $P$ asociază oricărui eveniment un număr între 0 și 1.
Dacă avem informații suplimentare, calculăm **probabilitatea condiționată**:
$$P(E | F) = \frac{P(E \cap F)}{P(F)}$$

### 🩸 Exemplu: Boala Huntington
Boala Huntington este cauzată de o mutație dominantă. Indivizii pot fi $AA$ (homozigot mutant), $Aa$ (heterozigot) sau $aa$ (normal). Atât $AA$ cât și $Aa$ fac boala.
Dacă încrucișăm doi părinți heterozigoți ($Aa \times Aa$), obținem o plantă/individ cu flori roz (bolnav). Care este probabilitatea ca acesta să fie homozigot mutant ($AA$)?
- Spațiul reduș (doar bolnavii): $\{AA, Aa, aA\}$ (3 cazuri).
- Cazul favorabil ($AA$): 1 caz.
$$P(AA | \text{Bolnav}) = \frac{1}{3}$$
"""

# ╔═╡ a1b2c3d4-0012-4a00-8000-000000000018
md"""
### 🔄 4.2. Independența și Legea Probabilității Totale
Două evenimente sunt **independente** dacă $P(E \cap F) = P(E) \cdot P(F)$.

Dacă putem partiționa spațiul de selecție în evenimente disjuncte $F_1, F_2, \dots$, atunci:
$$P(E) = P(E|F_1)P(F_1) + P(E|F_2)P(F_2) + \dots$$

**Exercițiu:** Un bărbat cu Huntington (genotip $AA$ cu prob. 5% sau $Aa$ cu prob. 95%) și o femeie normală ($aa$) au un copil. Care e probabilitatea ca copilul să fie bolnav?
- Dacă tatăl e $AA$ ($F_1$), copilul primește $A$ $\implies$ sigur bolnav ($P(E|F_1) = 1$).
- Dacă tatăl e $Aa$ ($F_2$), copilul primește $A$ sau $a$ $\implies$ 50% bolnav ($P(E|F_2) = 0.5$).
$$P(E) = 1 \cdot 0.05 + 0.5 \cdot 0.95 = 0.525$$
"""

# ╔═╡ a1b2c3d4-0013-4a00-8000-000000000019
md"""
---
## 📦 5. Variabile Aleatoare și Distribuții (25 min)

O **variabilă aleatoare** $X$ asociază un număr real fiecărui rezultat al unui experiment.
"""

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000020
md"""
### 🪙 5.1. Distribuția Binomială
Modelează numărul de succese în $n$ experimente independente, fiecare cu probabilitatea $p$.
"""

# ╔═╡ 014d2112-0077-4ad5-812a-07ed451574af
md"""
🧬 **Exemplu Biologic: Genetica Mendeliană și Epidemiologia**
Imaginați-vă doi părinți care sunt purtători sănătoși (heterozigoți, $Aa$) pentru o boală recesivă (ex: fibroza chistică). Probabilitatea ca un copil să fie bolnav (homozigot $aa$) este $p = 0.25$. 
Dacă familia are $n = 4$ copii, distribuția binomială ne răspunde la întrebarea: *Care este probabilitatea ca exact $k=1$ copil să fie bolnav?*
$$P(X=1) = \binom{4}{1} 0.25^1 \cdot 0.75^3 \approx 0.42$$

**Alt exemplu:** Dintr-un eșantion de $n=100$ de pacienți vaccinați, dacă eficacitatea vaccinului este de 95% ($p=0.95$), binomiala prezice câți pacienți vor fi protejați și câți se vor infecta totuși.
"""

# ╔═╡ a1b2c3d4-0014-4a00-8000-000000000021
md"""
Numărul de experimente: $n = {}$ $(@bind n_bin PlutoUI.Slider(1:50, default=20, show_value=true))
"""

# ╔═╡ a1b2c3d4-0015-4a00-8000-000000000022
md"""
Probabilitatea de succes: $p = {}$ $(@bind p_bin PlutoUI.Slider(0.01:0.01:0.99, default=0.5, show_value=true))
"""

# ╔═╡ a1b2c3d4-0016-4a00-8000-000000000023
let
    d = Binomial(n_bin, p_bin)
    k_vals = 0:n_bin
    pmf_vals = pdf.(d, k_vals)
    
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1,1], title="Distribuția Binomială (n=$n_bin, p=$p_bin)",
              xlabel="Număr de succese (k)", ylabel="P(X = k)")
    
    barplot!(ax, k_vals, pmf_vals, color=:teal, strokecolor=:black)
    fig
end

# ╔═╡ a1b2c3d4-0017-4a00-8000-000000000024
md"""
### 🦠 5.2. Distribuția Poisson
Modelează numărul de evenimente rare într-un interval fix (ex: substituții amino-acizi). Parametrul $\lambda$ este rata medie.
"""

# ╔═╡ 79dd87de-3d21-476c-a95b-1eec06f8a8ab
md"""
🦠 **Exemplu Biologic: Mutații și Microbiologie**
Distribuția Poisson este „regele” evenimentelor rare în biologie. 

1. **Genetica moleculară:** Dacă rata medie de mutație pentru o genă este de $\lambda = 2$ mutații per $10^6$ diviziuni celulare, Poisson ne spune probabilitatea de a găsi exact $k=0$ (niciun mutant) sau $k=3$ (trei mutanți) într-o cultură de bacterii.

2. **Microbiologie:** Numărul de colonii de bacterii care cresc pe o placă Petri după ce ai împrăștiat un volum fix de $10 \mu L$ dintr-o cultură diluată.

3. **Neurobiologie:** Numărul de potențiale de acțiune (spike-uri) pe care un neuron le generează într-o fereastră de timp foarte scurtă (ex: 10 milisecunde) sub un stimul constant.
"""

# ╔═╡ a1b2c3d4-0018-4a00-8000-000000000025
md"""
Rata medie: $\lambda = {}$ $(@bind λ_pois PlutoUI.Slider(0.1:0.1:20.0, default=5.0, show_value=true))
"""

# ╔═╡ a1b2c3d4-0019-4a00-8000-000000000026
let
    d = Poisson(λ_pois)
    k_max = max(20, ceil(Int, λ_pois + 4*sqrt(λ_pois)))
    k_vals = 0:k_max
    pmf_vals = pdf.(d, k_vals)
    
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1,1], title="Distribuția Poisson (λ=$λ_pois)",
              xlabel="Evenimente (k)", ylabel="P(X = k)")
    
    barplot!(ax, k_vals, pmf_vals, color=:purple, strokecolor=:black)
    vlines!(ax, [λ_pois], color=:red, linestyle=:dash, linewidth=2, label="Media = $λ_pois")
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0020-4a00-8000-000000000027
md"""
### 📏 5.3. Distribuția Normală (Gaussiană)
Definită de media $\mu$ și deviația standard $\sigma$. Este fundamentul Teoremei Limită Centrale.
"""

# ╔═╡ a481b8e4-c36c-4384-bd8c-dbd2bc195035
md"""
📏 **Exemplu Biologic: Trăsături Poligenice și Erori de Măsurare**
Datorită Teoremei Limită Centrale, Normala apare oriunde un fenotip este rezultatul sumei a zeci de factori mici și independenți.
1. **Fiziologie umană:** Tensiunea arterială sistolică a unei populații adulte sănătoase, greutatea la naștere a nou-născuților sau înălțimea adulților. Acestea sunt *trăsături poligenice* (influențate de sute de gene și factori de mediu).
2. **Biochimie:** Concentrația unui metabolit (ex: glucoza sau colesterolul) în sângele unei populații. Valorile extreme (pe „cozile” clopotului) reprezintă indivizii bolnavi (ex: diabet, hipercolesterolemie).
3. **Laborator:** Erorile de măsurare ale unui spectrofotometru când măsoară absorbția unei probe. Media erorilor este zero, iar ele se distribuie normal în jurul valorii reale.
"""

# ╔═╡ a1b2c3d4-0021-4a00-8000-000000000028
md"""
Media: $\mu = {}$ $(@bind μ_norm PlutoUI.Slider(-10.0:0.5:10.0, default=0.0, show_value=true))
"""

# ╔═╡ a1b2c3d4-0022-4a00-8000-000000000029
md"""
Deviația standard: $\sigma = {}$ $(@bind σ_norm PlutoUI.Slider(0.5:0.1:5.0, default=1.0, show_value=true))
"""

# ╔═╡ a1b2c3d4-0023-4a00-8000-000000000030
let
    d = Normal(μ_norm, σ_norm)
    xs = LinRange(μ_norm - 4*σ_norm, μ_norm + 4*σ_norm, 500)
    ys = pdf.(d, xs)
    
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1,1], title="Distribuția Normală (μ=$μ_norm, σ=$σ_norm)",
              xlabel="x", ylabel="Densitate f(x)")
    
    lines!(ax, xs, ys, color=:orange, linewidth=3)
    
    vspan!(ax, μ_norm - 2*σ_norm, μ_norm + 2*σ_norm; color=(:blue, 0.15), label="95.4% (±2σ)")
    vspan!(ax, μ_norm - σ_norm, μ_norm + σ_norm; color=(:orange, 0.3), label="68.2% (±1σ)")
    
    axislegend(ax, position=:rt)
    fig
end

# ╔═╡ a1b2c3d4-0024-4a00-8000-000000000031
md"""
---
## ⚖️ 6. Legi Limită (10 min)

### 📈 6.1. Legea Numerelor Mari (LLN)
Media eșantionului $\overline{X}_n$ tinde către media teoretică $\mu$ pe măsură ce $n \to \infty$.
*Exemplu:* Un studiu pe 429 de pacienți dă o estimare mult mai precisă a prevalenței unei boli decât un studiu pe 30 de pacienți.

### 🌌 6.2. Teorema Limită Centrală (CLT)
Chiar dacă populația originală **nu** este normală, distribuția mediilor de eșantionare va fi aproximativ **Normală**. Aceasta este "puntea" care leagă statistica descriptivă de inferența statistică.
"""

# ╔═╡ a1b2c3d4-0025-4a00-8000-000000000032
md"""
---
## 📝 7. Concluzii și Sinteză

### 📊 Tabel Rezumtiv

| Concept | Definiție / Utilizare |
| :--- | :--- |
| **Tabel de contingență** | Compararea a două variabile categoriale (ex: Tratament vs Infecție). |
| **Combinări ($C_n^k$)** | Numărul de moduri de a alege $k$ elemente din $n$ (fără ordine). |
| **Probabilitate condiționată** | Probabilitatea lui $E$ știind că $F$ s-a produs. |
| **Binomială / Poisson** | Distribuții discrete pentru succese sau evenimente rare. |
| **Normală (CLT)** | Distribuție continuă, universală pentru mediile eșantioanelor. |

### 💡 Întrebări de reflecție
1. De ce este importantăIndependența evenimentelor atunci când aplicăm Legea Probabilității Totale în studiile clinice?
2. Cum ne ajută Distribuția Poisson să modelăm mutațiile genetice rare, comparativ cu Binomiala?
3. În ce condiții un grafic de regresie liniară nu este suficient pentru a demonstra o relație cauzală între două variabile numerice?

---
**✨ Sfârșitul Cursului 11 ✨**
"""

# ╔═╡ Cell order:
# ╟─a1b2c3d4-0000-4a00-8000-000000000000
# ╟─a1b2c3d4-0000-4a00-8000-000000000001
# ╟─a1b2c3d4-0001-4a00-8000-000000000001
# ╠═a1b2c3d4-0002-4a00-8000-000000000002
# ╟─28496985-d587-466f-ba69-8db8a67635e9
# ╟─0684c3db-fdd9-4c91-9801-f3bac7bfa73a
# ╟─a1b2c3d4-0003-4a00-8000-000000000003
# ╠═a1b2c3d4-0004-4a00-8000-000000000004
# ╟─a1b2c3d4-0005-4a00-8000-000000000005
# ╠═a1b2c3d4-0006-4a00-8000-000000000006
# ╟─a1b2c3d4-0007-4a00-8000-000000000007
# ╠═a1b2c3d4-0008-4a00-8000-000000000008
# ╠═a1b2c3d4-0009-4a00-8000-000000000009
# ╟─a1b2c3d4-0010-4a00-8000-000000000010
# ╟─a1b2c3d4-0010-4a00-8000-000000000011
# ╠═a1b2c3d4-0010-4a00-8000-000000000012
# ╟─a1b2c3d4-0010-4a00-8000-000000000013
# ╠═a1b2c3d4-0011-4a00-8000-000000000014
# ╟─a1b2c3d4-0012-4a00-8000-000000000015
# ╠═a1b2c3d4-0012-4a00-8000-000000000016
# ╟─a1b2c3d4-0012-4a00-8000-000000000017
# ╟─a1b2c3d4-0012-4a00-8000-000000000018
# ╟─a1b2c3d4-0013-4a00-8000-000000000019
# ╟─a1b2c3d4-0014-4a00-8000-000000000020
# ╟─014d2112-0077-4ad5-812a-07ed451574af
# ╟─a1b2c3d4-0014-4a00-8000-000000000021
# ╟─a1b2c3d4-0015-4a00-8000-000000000022
# ╠═a1b2c3d4-0016-4a00-8000-000000000023
# ╟─a1b2c3d4-0017-4a00-8000-000000000024
# ╟─79dd87de-3d21-476c-a95b-1eec06f8a8ab
# ╟─a1b2c3d4-0018-4a00-8000-000000000025
# ╠═a1b2c3d4-0019-4a00-8000-000000000026
# ╟─a1b2c3d4-0020-4a00-8000-000000000027
# ╟─a481b8e4-c36c-4384-bd8c-dbd2bc195035
# ╟─a1b2c3d4-0021-4a00-8000-000000000028
# ╟─a1b2c3d4-0022-4a00-8000-000000000029
# ╠═a1b2c3d4-0023-4a00-8000-000000000030
# ╟─a1b2c3d4-0024-4a00-8000-000000000031
# ╟─a1b2c3d4-0025-4a00-8000-000000000032
