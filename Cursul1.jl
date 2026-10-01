### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> title = "Cursul I"
#> description = "Sisteme cu feedback și funcții"
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

# ╔═╡ a8e79b63-4624-4f3c-a985-5bf6a1f51d43
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoTeachingTools
	using PlutoUI
end

# ╔═╡ 3c523909-557c-4473-9736-4e271c3e958a
TableOfContents()

# ╔═╡ ccc7e71a-1ce3-4bf0-92fb-9906287c5176


# ╔═╡ b1c2d3e4-5f6a-7b8c-9d0e-1f2a3b4c5d6e
HTML("""
<style>
    :root {
        --bg-color: #282a36;
        --text-color: #f8f8f2;
        --accent-color: #8be9fd;
        --code-bg: #44475a;
        --code-text: #ff79c6;
        --quote-bg: #343746;
        --quote-border: #50fa7b;
        --quote-text: #f1fa8c;
    }
    body {
        background-color: var(--bg-color) !important;
        color: var(--text-color) !important;
    }
    main {
        max-width: 1900px !important;
    }
    .pluto-output {
        background-color: transparent !important;
        color: var(--text-color) !important;
    }
    .pluto-output h1, .pluto-output h2, .pluto-output h3 {
        color: var(--accent-color) !important;
        border-bottom: 2px solid #44475a !important;
    }
    .pluto-output h1 { border-bottom: 3px solid var(--accent-color) !important; }
    .pluto-output blockquote {
        border-left: 5px solid var(--quote-border) !important;
        background: var(--quote-bg) !important;
        color: var(--quote-text) !important;
    }
    .pluto-output code {
        background: var(--code-bg) !important;
        color: var(--code-text) !important;
    }
    .pluto-output a { color: #bd93f9 !important; }
    .pluto-output img {
        border-radius: 12px;
        box-shadow: 0 8px 16px rgba(0,0,0,0.5);
        margin: 20px 0;
        transition: transform 0.2s;
        max-width: 100%;
    }
    .pluto-output img:hover { transform: scale(1.02); }
    .pluto-output ul li::marker { color: var(--accent-color); }
    .pluto-cell { background-color: transparent !important; }
</style>
""")

# ╔═╡ 1855e998-fc81-4de3-b833-e903669b6f03
md"""---
### 📝 Preliminarii-Cum se va face notarea
- **20 puncte** pentru activitate (răspunsuri, teme etc.);
- **30 puncte** un proiect (poate fi făcut în echipe de maxim 2 persoane);
- **70 puncte** examenul final.

> 💡 *Pentru cei care doresc va fi un examen parțial în săptămîna a 8-a.*
> **Condiții de trecere:** cel puțin 50 de puncte în total și cel puțin 35 de puncte la examen.
"""

# ╔═╡ 215a9db0-9667-11f0-3e8b-db501bdbd2d5
md"""
# 🌿 Cursul I: Sisteme cu feedback și funcții
**andrei.halanay@unibuc.ro**
[📓 Notebook-uri](https://github.com/adhalanay/curs_ecologie2026)

"""

# ╔═╡ 8c1d2e3f-4a5b-4c6d-8e7f-9a0b1c2d3e4f
md"""
## 🎯 Obiective
- 🔄 Diferența dintre feedback pozitiv și feedback negativ;
- 🧬 Exemple de sisteme biologice cu feedback;
- 📐 Noțiunile de funcție, domeniu, codomeniu și imagine;
- 🔗 Compunerea a două funcții;
- 📊 Grafice obișnuite, semi-logaritmice și logaritmice;
- 🌌 Spațiul stărilor și dimensiunea acestuia;
- ➕ Operații cu vectori
"""

# ╔═╡ 11efb0cf-73b4-4d19-b0f7-612a5abec2d3
md"""
## 🔄 1. Sisteme cu feedback
Modelarea matematică a proceselor biologice a început să fie folosită sistematic de la mijlocul secolului XIX. Astăzi matematica este esențială în:
- 🛡️ imunologie și boli autoimune;
- 💊 farmacologie și proiectarea medicamentelor;
- 🦠 dinamica populațiilor;
- 🌍 epidemiologie;
- 🧠 neuroștiință;
- 🩻 imagistică medicală.

Un exemplu clasic de studiu ecologic este oscilația dintre:
- 🐈 **râsul canadian** (*Lynx canadensis*);
- 🐇 **iepurele de zăpadă** (*Lepus americanus*).

Râsul se hrănește aproape exclusiv cu iepurele de zăpadă.
"""

# ╔═╡ 9f282b4f-da9c-49b1-94ef-03cbe39ca1f9
PlutoUI.LocalResource("Canadalynx.jpg", :width => 300)

# ╔═╡ 4d08da5d-e63f-49d9-a25a-05c9ff0c2cd4
md"""
**Iepurele de zăpadă:**
"""

# ╔═╡ 9968fa47-4085-4417-8100-7c8385175aca
PlutoUI.LocalResource("SnowshoeHare.jpg", :width => 300)

# ╔═╡ ccd56781-7fbe-48ed-9617-6c88829384c0
md"""
Datele istorice provin din pieile colectate de Hudson's Bay Company timp de aproape 100 de ani. Numărul de râși este cunoscut, numărul de iepuri a fost estimat.
"""

# ╔═╡ d42b9ecf-086a-47fe-b2e8-742a3b0ed214
PlutoUI.LocalResource("fourrures.jpg")

# ╔═╡ 474aa1f8-42b2-44d7-a7ef-441d0ddffa36
md"""
**Pe termen scurt:**
"""

# ╔═╡ a6e282d9-dd41-440d-a8c6-bc068f23720a
PlutoUI.LocalResource("rîși.png")

# ╔═╡ 391dd0e7-66e6-42cb-a967-0345457f85f7
md"""
### 🔍 Se observă:
- oscilații cu perioadă de aproximativ 10 ani;
- populația de prădători crește/scade cu întârziere față de cea a prăzii.

Pentru a explica aceste oscilații avem nevoie de un model matematic.
**Ideea centrală:** fiecare populație determină valoarea celeilalte populații.

---
### 📈 Feedback pozitiv
Avem feedback pozitiv dacă o valoare pozitivă a unei variabile determină creșterea acelei variabile.
**Exemple:**
- 💰 banii investiți pot produce mai mulți bani;
- 🐣 animalele fac pui, ceea ce crește populația, deci numărul de pui;
- 🌡️ creșterea nivelului de $CO_2$ crește temperatura, ceea ce accelerează descompunerea materiei organice și produce mai mult $CO_2$.

### 📉 Feedback negativ
Avem feedback negativ dacă o valoare pozitivă a unei variabile determină scăderea ei, iar o valoare negativă determină creșterea ei.
**Exemplu:** aerul condiționat.
Fie `T_0` temperatura setată și `C` temperatura curentă. Definim: $T = C - T_0$.
- dacă `T > 0`, sistemul răcește;
- dacă `T < 0`, sistemul încălzește.

Un alt exemplu clasic: 🩸 **glucoza și insulina**.
Când glucoza crește, insulina crește, iar glicemia scade.
"""

# ╔═╡ 9dd94e6d-04e5-4669-adff-414b25cb1dd3
PlutoUI.LocalResource("insulină.png")

# ╔═╡ 291381ae-569e-4335-9110-b3af9bee2c9f
md"""
Comportamentul este asemănător cu sistemul râși/iepuri.
În epidemiologie, contactele dintre infectați și susceptibili cresc numărul de infectați și scad numărul de susceptibili. Modelele epidemiologice ajută la elaborarea strategiilor de intervenție.

### ⚠️ Comportamente neintuitive
Sistemele cu feedback au adesea efecte neintuitive.

> **Exemplu:** vrem să reducem numărul de prădători. Dacă scoatem prădători din mediu, prada crește, ceea ce poate duce la creșterea numărului de prădători peste valoarea inițială. Acesta este fenomenul de *rebound*.

Rezultatul unei intervenții depinde de **faza ciclului** în care are loc intervenția.

**Concluzie:** sisteme aparent simple pot avea comportamente foarte neintuitive. De aceea modelarea matematică este indispensabilă.
"""

# ╔═╡ 839ddca1-d3ba-4d99-b15f-a795b040c161
md"""
## 📊 2. Funcții și graficele lor
Graficul unei funcții $f : D \to Y$ este mulțimea tuturor perechilor (input, output):

$$G_f = \{(x, f(x)) \;:\; x \in D\} \subset \mathbb{R}^2.$$

Graficul obișnuit folosește două axe perpendiculare:
- axa orizontală $Ox$: variabila independentă $x$ (input-ul);
- axa verticală $Oy$: valorile funcției $y=f(x)$ (output-ul).

Fiecare punct al graficului este o pereche $(x, f(x))$: mergem pe orizontală pînă la $x$, apoi pe verticală pînă la înălțimea $f(x)$.

---
### 🔢 De la tabel la grafic
Reprezentarea numerică (tabelul) și cea vizuală (graficul) sunt legate direct: fiecare linie din tabel devine un punct în plan.

**Exemplu:** $f(x) = x^5$

| $x$ | $-1.5$ | $-1$ | $-0.5$ | $0$ | $0.5$ | $1$ | $1.5$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $f(x)$ | $-7.59$ | $-1$ | $-0.03$ | $0$ | $0.03$ | $1$ | $7.59$ |

Cu puține puncte ghicim doar forma curbei. Cu cît calculăm mai multe valori, cu atît punctele se unesc într-o curbă continuă. Așa desenează și calculatorul graficele: calculează sute de puncte și le unește prin segmente.
"""

# ╔═╡ 778d276d-9c3c-49c7-9d44-14baf3ae16e7
begin
	culoare_fundal = "#282a36"
	stil_ax = (
		backgroundcolor = culoare_fundal,
		xgridcolor = :grey30, ygridcolor = :grey30,
		xtickcolor = :white, ytickcolor = :white,
		xticklabelcolor = :white, yticklabelcolor = :white,
		xlabelcolor = :white, ylabelcolor = :white,
		titlecolor = :white, subtitlecolor = :grey80,
		leftspinecolor = :grey60, rightspinecolor = :grey60,
		topspinecolor = :grey60, bottomspinecolor = :grey60,
	)
	stil_leg = (backgroundcolor = (:black, 0.3), framecolor = :grey60, labelcolor = :white)
	figura(; kw...) = Figure(; backgroundcolor = culoare_fundal, kw...)
end

# ╔═╡ e4d39336-7ddd-4619-9288-6169aadd6a0e
begin
	f(x) = x^5
	xs = range(-1.5, 1.5, length=300)
	ys = f.(xs)
	xs_tabel = -1.5:0.5:1.5
	fig = figura(size=(1000, 400))
	ax_a = Axis(fig[1,1]; stil_ax..., xlabel="x", ylabel="f(x)",
		title="7 puncte din tabel")
	scatter!(ax_a, xs_tabel, f.(xs_tabel), color="#ff79c6", markersize=14,
		strokewidth=2, strokecolor=:white)
	ax = Axis(fig[1,2]; stil_ax..., xlabel="x", ylabel="f(x)",
		title="300 de puncte unite: graficul lui f(x) = x⁵")
	lines!(ax, xs, ys, color="#8be9fd", linewidth=3)
	scatter!(ax, xs_tabel, f.(xs_tabel), color="#ff79c6", markersize=10)
	linkaxes!(ax_a, ax)
	fig
end

# ╔═╡ 8c24d3d5-bc4f-4f41-a7da-c8d4a4be8077
md"""
### 👀 Ce putem citi pe un grafic
Graficul arată dintr-o privire proprietăți care din formulă se văd greu:
- **monotonia:** funcția este *crescătoare* dacă graficul urcă de la stînga la dreapta și *descrescătoare* dacă coboară;
- **zerourile:** punctele în care graficul taie axa $Ox$, adică soluțiile ecuației $f(x)=0$;
- **semnul:** $f(x)>0$ unde graficul este deasupra axei $Ox$, $f(x)<0$ unde este dedesubt;
- **extremele:** punctele cele mai înalte (maxime) sau cele mai joase (minime);
- **simetria:** față de axa $Oy$ (funcție *pară*, $f(-x)=f(x)$) sau față de origine (funcție *impară*, $f(-x)=-f(x)$).

Pentru $f(x)=x^5$: funcția este crescătoare, are un singur zero ($x=0$), este negativă pentru $x<0$ și pozitivă pentru $x>0$, iar graficul este simetric față de origine, deci $f$ este impară.

> 💡 Observați că pe intervalul $[-0.5, 0.5]$ graficul pare aproape lipit de axă: valorile $x^5$ sunt foarte mici acolo ($0.5^5 \approx 0.03$). Un grafic depinde mult de **scara** aleasă; revenim la asta la scara logaritmică.

---
### 📐 Domeniul și imaginea pe grafic
Pe grafic, domeniul și imaginea se văd ca **proiecții**:
- **domeniul** este „umbra” graficului pe axa $Ox$;
- **imaginea** este „umbra” graficului pe axa $Oy$.

**Exemplu:** $f : [-1, 2] \to \mathbb{R}$, $f(x) = x^2$. Domeniul este $[-1,2]$, codomeniul este $\mathbb{R}$, dar imaginea este doar $[0, 4]$: valorile negative nu sînt atinse, iar cea mai mare valoare este $f(2)=4$.
"""

# ╔═╡ e4d67807-c4cb-4573-a8a9-2b7ba7e59fc9
let
	g(x) = x^2
	xg = range(-1, 2, length=300)
	fig5 = figura(size=(600, 450))
	ax5 = Axis(fig5[1,1]; stil_ax..., xlabel="x", ylabel="f(x)",
		title="f(x) = x² pe [-1, 2]",
		subtitle="domeniu [-1, 2] (verde) · imagine [0, 4] (portocaliu)")
	# proiecțiile punctelor extreme pe axe
	for (x0, y0) in [(-1, 1), (2, 4)]
		lines!(ax5, [x0, x0], [0, y0], color=:grey70, linestyle=:dash)
	end
	lines!(ax5, [0, 2], [4, 4], color=:grey70, linestyle=:dash)
	lines!(ax5, xg, g.(xg), color="#8be9fd", linewidth=3)
	# domeniul pe Ox și imaginea pe Oy
	lines!(ax5, [-1, 2], [0, 0], color="#50fa7b", linewidth=8)
	lines!(ax5, [0, 0], [0, 4], color="#ffb86c", linewidth=8)
	scatter!(ax5, [-1, 2], [1, 4], color="#ff79c6", markersize=12)
	limits!(ax5, -1.6, 2.4, -0.6, 4.6)
	fig5
end

# ╔═╡ 928d4fd2-74ad-48a5-8c83-b24156dca337
md"""
---
### ✋ Testul dreptei verticale
O curbă din plan este graficul unei funcții $y = f(x)$ dacă și numai dacă **orice dreaptă verticală o intersectează în cel mult un punct**. Motivul: dreapta verticală $x = x_0$ taie graficul exact în punctele $(x_0, y)$, iar o funcție asociază fiecărui $x_0$ un singur $y$.

- Parabola $y = x^2$ trece testul: fiecare verticală o taie o singură dată.
- Relația $y^2 = x^2$ nu trece testul: dreapta $x = 1$ o taie în $(1, 1)$ și în $(1, -1)$. Deci $y^2=x^2$ **nu** definește o funcție $y = f(x)$.
"""

# ╔═╡ 05ecb7e4-89fc-4265-8b43-0fc6208d7e8e
let
	xg = range(-2, 2, length=300)
	fig6 = figura(size=(1000, 400))
	ax6a = Axis(fig6[1,1]; stil_ax..., xlabel="x", ylabel="y",
		title="y = x²: este funcție ✔")
	lines!(ax6a, xg, xg .^ 2, color="#8be9fd", linewidth=3)
	vlines!(ax6a, 1, color="#f1fa8c", linestyle=:dash, linewidth=2)
	scatter!(ax6a, [1], [1], color="#50fa7b", markersize=16,
		strokewidth=2, strokecolor=:white)
	ax6b = Axis(fig6[1,2]; stil_ax..., xlabel="x", ylabel="y",
		title="y² = x²: nu este funcție ✘")
	lines!(ax6b, xg, xg, color="#8be9fd", linewidth=3)
	lines!(ax6b, xg, -xg, color="#8be9fd", linewidth=3)
	vlines!(ax6b, 1, color="#f1fa8c", linestyle=:dash, linewidth=2)
	scatter!(ax6b, [1, 1], [1, -1], color="#ff5555", markersize=16,
		strokewidth=2, strokecolor=:white)
	fig6
end

# ╔═╡ 284cd1ac-987b-46dc-a066-e5ac244684e1
md"""
---
### 🧰 Funcții elementare: explorare interactivă
Majoritatea modelelor din acest curs sînt construite din cîteva funcții de bază. Alegeți o familie și modificați parametrul $a$ pentru a vedea cum se schimbă graficul, domeniul și imaginea.

**Familia:** $(@bind familie PlutoUI.Select([
	"liniara" => "Liniară: f(x) = a·x",
	"putere" => "Putere: f(x) = xᵃ (x > 0)",
	"exp" => "Exponențială: f(x) = aˣ",
	"log" => "Logaritmică: f(x) = logₐ(x)",
]))

**Parametrul a:** $(@bind param_a PlutoUI.Slider(0.1:0.1:3, default=2, show_value=true))

De urmărit:
- la exponențială, ce se întîmplă pentru $a < 1$ față de $a > 1$? Dar pentru $a = 1$?
- la putere, cum diferă graficele pentru $a < 1$ și $a > 1$? Toate trec prin punctul $(1, 1)$. De ce?
- exponențiala și logaritmul cu aceeași bază $a$ sînt funcții *inverse*: graficele lor sînt simetrice față de dreapta $y = x$.
"""

# ╔═╡ e4ff1401-fe97-43b5-bb9a-85df5b0310a4
let
	a = param_a
	if familie == "liniara"
		g = x -> a * x
		xg = range(-3, 6, length=400)
		titlu = "f(x) = $(a)·x"
		info = "domeniu ℝ · imagine ℝ"
	elseif familie == "putere"
		g = x -> x^a
		xg = range(0.001, 6, length=400)
		titlu = "f(x) = x^$(a)"
		info = "domeniu (0, ∞) · imagine (0, ∞)"
	elseif familie == "exp"
		g = x -> a^x
		xg = range(-3, 6, length=400)
		titlu = "f(x) = $(a)ˣ"
		info = a == 1 ? "domeniu ℝ · imagine {1} (funcție constantă)" :
			"domeniu ℝ · imagine (0, ∞)"
	else
		g = x -> a == 1 ? NaN : log(x) / log(a)
		xg = range(0.001, 6, length=400)
		titlu = "f(x) = logₐ(x),  a = $(a)"
		info = a == 1 ? "baza a = 1 nu este permisă" :
			"domeniu (0, ∞) · imagine ℝ"
	end
	fig7 = figura(size=(700, 450))
	ax7 = Axis(fig7[1,1]; stil_ax..., xlabel="x", ylabel="f(x)",
		title=titlu, subtitle=info)
	hlines!(ax7, 0, color=:grey60)
	vlines!(ax7, 0, color=:grey60)
	lines!(ax7, xg, g.(xg), color="#8be9fd", linewidth=3)
	limits!(ax7, -3, 6, -4, 8)
	fig7
end

# ╔═╡ 6dcfeafd-2cd9-41a0-83c6-d297f1fa2298
md"""
---
### 🔀 Transformări ale graficului
Pornind de la graficul lui $f$, putem obține graficul lui

$$g(x) = A \cdot f(x - h) + k$$

prin operații geometrice simple:

| Parametru | Efect asupra graficului |
|:---:|:---|
| $h > 0$ / $h < 0$ | translație la **dreapta** / la **stînga** cu $\lvert h \rvert$ |
| $k > 0$ / $k < 0$ | translație în **sus** / în **jos** cu $\lvert k \rvert$ |
| $\lvert A \rvert > 1$ / $\lvert A \rvert < 1$ | **dilatare** / **comprimare** pe verticală |
| $A < 0$ | **reflexie** față de axa $Ox$ |

> ⚠️ Atenție la semn: $f(x - 2)$ mută graficul spre **dreapta**, nu spre stînga. Valoarea pe care $f$ o lua în $0$ este acum luată în $x = 2$.

Luăm $f(x) = x^2$ (linia punctată) și modificăm parametrii:

``A`` = $(@bind transf_A PlutoUI.Slider(-2:0.25:2, default=1, show_value=true))
``h`` = $(@bind transf_h PlutoUI.Slider(-3:0.5:3, default=0, show_value=true))
``k`` = $(@bind transf_k PlutoUI.Slider(-3:0.5:3, default=0, show_value=true))
"""

# ╔═╡ 787799ae-e91c-409b-aba7-58aaa683c0cc
let
	xg = range(-5, 5, length=400)
	g(x) = transf_A * (x - transf_h)^2 + transf_k
	fig8 = figura(size=(700, 450))
	ax8 = Axis(fig8[1,1]; stil_ax..., xlabel="x", ylabel="y",
		title="g(x) = $(transf_A)·(x − $(transf_h))² + $(transf_k)",
		subtitle="vîrful parabolei: ($(transf_h), $(transf_k))")
	hlines!(ax8, 0, color=:grey60)
	vlines!(ax8, 0, color=:grey60)
	lines!(ax8, xg, xg .^ 2, color=:grey70, linestyle=:dash, linewidth=2,
		label="f(x) = x²")
	lines!(ax8, xg, g.(xg), color="#ff79c6", linewidth=3, label="g(x)")
	scatter!(ax8, [transf_h], [transf_k], color="#f1fa8c", markersize=14,
		strokewidth=2, strokecolor=:white)
	limits!(ax8, -5, 5, -6, 8)
	axislegend(ax8; stil_leg..., position=:rb)
	fig8
end

# ╔═╡ 21314352-320b-4b2a-a925-2a6ece098140
md"""
---
### 🧬 Funcții întîlnite des în biologie
Cîteva forme de grafic apar în aproape toate modelele din acest curs. E util să le recunoașteți „din vedere”.

**1. Creșterea exponențială** $N(t) = N_0\, e^{rt}$: fiecare individ produce noi indivizi, deci populația crește cu atît mai repede cu cît e mai mare. Este semnătura **feedback-ului pozitiv**.

**2. Creșterea logistică** $N(t) = \dfrac{K}{1 + \frac{K - N_0}{N_0} e^{-rt}}$: la început arată ca o exponențială, apoi resursele limitate frînează creșterea și populația se stabilizează la **capacitatea de suport** $K$. Aici intervine un **feedback negativ**.

**3. Cinetica Michaelis–Menten** $v(S) = \dfrac{V_{max}\, S}{K_m + S}$: viteza unei reacții enzimatice în funcție de concentrația substratului $S$. Crește aproape liniar pentru $S$ mic și se saturează la $V_{max}$; pentru $S = K_m$ viteza este exact $V_{max}/2$.

**4. Funcția Hill** $h(S) = \dfrac{S^n}{K^n + S^n}$: pentru $n > 1$ graficul are formă de **S** (sigmoidă) și se comportă ca un „comutator”: răspunsul trece rapid de la aproape $0$ la aproape $1$ în jurul lui $S = K$. Exemplu: legarea oxigenului de hemoglobină ($n \approx 2.8$).
"""

# ╔═╡ 4a964e73-1e3d-4fd8-adf9-ec0830c00b25
let
	fig9 = figura(size=(1300, 420))

	N0, r, K = 10, 0.5, 500
	t = range(0, 16, length=300)
	ax9a = Axis(fig9[1,1]; stil_ax..., xlabel="timp t", ylabel="N(t)",
		title="Exponențial vs. logistic")
	lines!(ax9a, t, N0 .* exp.(r .* t), color="#ff79c6", linewidth=3,
		label="exponențial")
	lines!(ax9a, t, K ./ (1 .+ (K - N0) / N0 .* exp.(-r .* t)),
		color="#50fa7b", linewidth=3, label="logistic")
	hlines!(ax9a, K, color=:grey70, linestyle=:dash)
	text!(ax9a, 0.3, K, text="K = $K", color=:grey80, align=(:left, :bottom))
	ylims!(ax9a, 0, 1.6K)
	axislegend(ax9a; stil_leg..., position=:lt)

	Vmax, Km = 1.0, 2.0
	S = range(0, 20, length=300)
	ax9b = Axis(fig9[1,2]; stil_ax..., xlabel="substrat S", ylabel="v(S)",
		title="Michaelis–Menten")
	lines!(ax9b, S, Vmax .* S ./ (Km .+ S), color="#8be9fd", linewidth=3)
	hlines!(ax9b, Vmax, color=:grey70, linestyle=:dash)
	lines!(ax9b, [Km, Km, 0], [0, Vmax / 2, Vmax / 2], color="#f1fa8c",
		linestyle=:dot, linewidth=2)
	text!(ax9b, Km, 0.02, text=" Kₘ", color="#f1fa8c", align=(:left, :bottom))
	text!(ax9b, 19.5, Vmax, text="Vₘₐₓ", color=:grey80, align=(:right, :bottom))
	ylims!(ax9b, 0, 1.15)

	Kh = 2.0
	S2 = range(0, 6, length=300)
	ax9c = Axis(fig9[1,3]; stil_ax..., xlabel="S", ylabel="h(S)",
		title="Funcții Hill")
	for (n, c) in [(1, "#8be9fd"), (2, "#bd93f9"), (4, "#ff79c6"), (8, "#ffb86c")]
		lines!(ax9c, S2, S2 .^ n ./ (Kh^n .+ S2 .^ n), color=c, linewidth=3,
			label="n = $n")
	end
	vlines!(ax9c, Kh, color=:grey70, linestyle=:dash)
	axislegend(ax9c; stil_leg..., position=:rb)

	fig9
end

# ╔═╡ 5477506e-4926-41d0-a31c-7bf486040402
md"""
---
### 📏 Scara logaritmică
Uneori valorile lui $y$ acoperă multe ordine de mărime (de la $1$ la $1\,000\,000$). Pe o scară obișnuită valorile mici devin invizibile, așa cum am văzut la $x^5$ lîngă $0$. Folosim atunci scara logaritmică pentru una dintre axe sau pentru amîndouă.

Pe o axă logaritmică **distanțe egale corespund unor rapoarte egale**: de la $1$ la $10$ este aceeași distanță ca de la $10$ la $100$ sau de la $1000$ la $10\,000$. Practic, în loc de $y$ desenăm $\log_{10} y$.

**Exemple de scări logaritmice:**
- 🧪 pH: $\text{pH} = -\log_{10}[H^+]$, deci o unitate de pH înseamnă de 10 ori mai mulți ioni $H^+$;
- 🌍 scara Richter;
- 🔊 decibeli.

**De ce ne interesează?** Scările logaritmice transformă două tipuri importante de legi în **drepte**, iar dreptele se recunosc și se măsoară ușor.

| Lege | Formula | Logaritmăm | Grafic în care apare dreaptă |
|:---|:---:|:---:|:---:|
| exponențială | $y = a \cdot b^x$ | $\log y = \log a + x \log b$ | **semi-log** ($y$ logaritmic), panta $\log b$ |
| putere | $y = a \cdot x^k$ | $\log y = \log a + k \log x$ | **log-log** (ambele logaritmice), panta $k$ |

Comparăm $y = 2^x$ (exponențială) și $y = x^5$ (putere) pe cele trei tipuri de grafic:
"""

# ╔═╡ 019adb7d-8418-4149-a81b-ec99341efe50
let
	xs_log = range(1, 20, length=300)
	fig1 = figura(size=(1300, 420))
	scari = [
		("Scară liniară", identity, identity),
		("Semi-log: 2ˣ devine dreaptă", identity, log10),
		("Log-log: x⁵ devine dreaptă", log10, log10),
	]
	for (i, (titlu, sx, sy)) in enumerate(scari)
		ax1 = Axis(fig1[1,i]; stil_ax..., xscale=sx, yscale=sy,
			xlabel=sx == log10 ? "x (scară log)" : "x",
			ylabel=sy == log10 ? "y (scară log)" : "y", title=titlu,
			xticks=sx == log10 ? [1, 2, 5, 10, 20] : Makie.automatic)
		lines!(ax1, xs_log, 2 .^ xs_log, color="#ff79c6", linewidth=3,
			label="y = 2ˣ")
		lines!(ax1, xs_log, xs_log .^ 5, color="#50fa7b", linewidth=3,
			label="y = x⁵")
		i == 1 && axislegend(ax1; stil_leg..., position=:lt)
	end
	fig1
end

# ╔═╡ 803b1d6b-e26b-44eb-917d-d0265c228b14
md"""
**Observații:**
- pe scara liniară ambele curbe „explodează” și nu putem spune care lege este care;
- dacă graficul **semi-logaritmic** este o dreaptă, datele satisfac o **lege exponențială** (creștere de tip populație, dezintegrare radioactivă, eliminarea unui medicament);
- dacă graficul **log-log** este o dreaptă, datele satisfac o **lege putere**, iar panta dreptei este exponentul $k$.

---
### 🐘 Exemplu: legea lui Kleiber
Rata metabolică bazală $B$ a mamiferelor (energia consumată în repaus) depinde de masa corporală $M$ după o lege putere:

$$B \approx 70 \cdot M^{3/4} \quad (\text{kcal/zi},\ M \text{ în kg}).$$

Masele merg de la cîteva grame (șoarece) la cîteva tone (elefant), deci pe o scară obișnuită animalele mici se îngrămădesc toate în colțul din stînga jos. Pe graficul log-log legea devine o dreaptă cu panta $3/4$.
"""

# ╔═╡ 4f5a7238-8855-44b1-8afa-b774958caaf0
let
	kleiber(M) = 70 * M^(3/4)
	animale = ["șoarece" => 0.025, "pisică" => 4.0, "om" => 70.0,
		"cal" => 500.0, "elefant" => 4000.0]
	nume = first.(animale)
	M = last.(animale)
	Mg = 10 .^ range(-2, 4, length=300)
	fig2 = figura(size=(1100, 450))
	for (i, (sx, titlu)) in enumerate([(identity, "Scară liniară"),
			(log10, "Scară log-log: dreaptă cu panta 3/4")])
		ax2 = Axis(fig2[1,i]; stil_ax..., xscale=sx, yscale=sx,
			xlabel="masa M (kg)", ylabel="B (kcal/zi)", title=titlu)
		lines!(ax2, Mg, kleiber.(Mg), color="#50fa7b", linewidth=3)
		scatter!(ax2, M, kleiber.(M), color="#ff79c6", markersize=14,
			strokewidth=2, strokecolor=:white)
		# pe scara liniară animalele mici se suprapun: le etichetăm doar pe cele mari
		vizibile = i == 1 ? (M .>= 500) : trues(length(M))
		text!(ax2, M[vizibile], kleiber.(M[vizibile]), text=nume[vizibile],
			color=:white, align=(:left, :top), offset=(8, -6))
		i == 2 && (xlims!(ax2, 5e-3, 1e5); ax2.xticks = (10.0 .^ (-2:4), ["10⁻²", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴"]))
	end
	fig2
end

# ╔═╡ 39adb7af-b1f0-4cf0-93ea-be4063654c55
md"""
### ✅ Verificare rapidă
1. Pe un grafic semi-logaritmic, numărul de bacterii dintr-o cultură apare ca o dreaptă care urcă cu o diviziune (adică de 10 ori) la fiecare 3 ore. Ce formulă are $N(t)$? În cît timp se dublează populația?
"""



# ╔═╡ 17c910ec-7256-4ba0-8ee6-d29c1642703e
hint(md"Dreapta în semi-log înseamnă ``\log_{10} N = \log_{10} N_0 + t/3``, deci ``N(t) = N_0 \cdot 10^{t/3}``")

# ╔═╡ d8ac2f4c-cc9f-4e12-9275-c429dc65010b
md"""
Dublarea: $10^{t/3} = 2 \Rightarrow t = 3\log_{10} 2 \approx 0.9$ ore, adică aproximativ 54 de minute."))

2. Graficul funcției $g(x) = -(x+1)^2 + 3$ se obține din graficul lui $x^2$ prin ce transformări? Care este imaginea lui $g$?
"""


# ╔═╡ 1549ca87-3cb5-497b-b812-b66fc8e4d73a
hint(md"Translație cu 1 la stînga (``h=-1``), reflexie față de ``Ox`` (``A=-1``), translație cu 3 în sus (``k=3``). Vîrful este în ``(-1, 3)`` și parabola este cu deschiderea în jos, deci imaginea este ``(-\infty, 3]``.")

# ╔═╡ 7d5f4320-933e-4358-8322-dcf613a8e702
md"""
## 🌌 3. Spațiul stărilor
Variabilele de stare descriu cantitativ sistemul la un moment dat.
**Exemplu:** starea unei populații poate fi mărimea ei. Putem fi interesați și de:
- ⚥ raportul între sexe;
- 🎂 distribuția pe clase de vârstă;
- 🗺️ distribuția teritorială.

Alegerea variabilelor de stare este cea mai importantă și dificilă parte a modelării.
Variabilele de stare sunt funcții de timp.
Mulțimea tuturor valorilor posibile ale variabilelor de stare se numește **spațiul stărilor**.
"""

# ╔═╡ 920b2efd-f83f-4300-b34f-b9695cf1079b
md"""
## 📏 4. Sisteme de dimensiune mai mare
Dacă avem două variabile de stare, de exemplu râsi `R` și iepuri `I`, spațiul stărilor este format din perechi: $(R, I)$.

Putem defini:
- **adunarea pe componente:** $(R_1,I_1)+(R_2,I_2)=(R_1+R_2,I_1+I_2)$;
- **înmulțirea cu scalari:** $a(R,I)=(aR,aI)$.

Putem desena câte o axă pentru fiecare componentă.
"""

# ╔═╡ 5b1a6088-d669-4892-9ab4-1068dc68ec3f
begin
	fig3 = Figure(size=(600, 400))
	ax3 = Axis(fig3[1,1], 
		xlabel="Râsi (R)", 
		ylabel="Iepuri (I)", 
		title="Spațiul Stărilor: Râsi vs. Iepuri",
		backgroundcolor=:transparent,
		xgridcolor=:grey40, ygridcolor=:grey40,
		xtickcolor=:white, ytickcolor=:white,
		xlabelcolor=:white, ylabelcolor=:white, titlecolor=:white
	)
	scatter!([2,1,3,5], [3,2,4,1], color="#bd93f9", markersize=15, strokewidth=2, strokecolor=:white)
	fig3
end

# ╔═╡ 5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
md"""
O pereche de numere se numește **2-vector**. Fiecare număr este o componentă.
Un vector poate fi reprezentat și ca o săgeată:
"""

# ╔═╡ 4a3e5f6a-7b8c-4d9e-8f0a-1b2c3d4e5f6a
begin
	fig4 = Figure(size=(600, 400))
	ax4 = Axis(fig4[1,1], 
		xlabel="x", 
		ylabel="y", 
		title="Reprezentarea geometrică a unui vector",
		backgroundcolor=:transparent,
		xgridcolor=:grey40, ygridcolor=:grey40,
		xtickcolor=:white, ytickcolor=:white,
		xlabelcolor=:white, ylabelcolor=:white, titlecolor=:white
	)
	arrows2d!([0], [0], [1.985], [2.985], color="green", tailwidth=3, taillength=15)
	scatter!(1.985, 2.985, markersize=10, color="red", strokewidth=2, strokecolor=:white)
	fig4
end

# ╔═╡ 9d2e3f4a-5b6c-4d7e-8f9a-0b1c2d3e4f5a
md"""
## 📝 5. Recapitulare și exerciții

### 🧠 Recapitulare
- 🔄 **Feedback pozitiv:** variabila își amplifică creșterea.
- 🛑 **Feedback negativ:** variabila își limitează creșterea.
- ⚠️ Sistemele cu feedback pot avea *rebound* și dependență de fază.
- 🎯 O funcție asociază fiecărui input exact un output.
- 🗺️ Domeniul, codomeniul și imaginea sunt noțiuni diferite.
- 🎨 Funcțiile pot fi reprezentate verbal, numeric, vizual și algebric.
- 📏 Scara logaritmică este utilă când valorile diferă mult ca ordin de mărime.
- 📈 În grafic semi-log legile exponențiale devin drepte; în grafic log-log legile putere devin drepte.
- ⏳ Variabilele de stare sunt funcții de timp.
- 🌌 Spațiul stărilor este mulțimea valorilor posibile ale variabilelor de stare.
- ➡️ Pentru mai multe variabile folosim vectori în $\mathbb{R}^n$.

---
### ✍️ Exerciții
1. Dați un exemplu de feedback pozitiv și unul de feedback negativ din biologie.
2. Explicați fenomenul de *rebound* într-un sistem prădător-pradă.
3. Pentru $f(x)=\sqrt{x}$ și $g(x)=x^2+1$, calculați $f\circ g$ și $g\circ f$. Precizați domeniile.
4. Desenați graficul lui $f(x)=2^x$ și apoi graficul semi-logaritmic. Ce observați?
5. Care este spațiul stărilor pentru:
   - masa unui organism;
   - concentrația unei substanțe;
   - un procent?
6. Calculați: $(2,3)+(-1,4), \quad 3(2,-1)$.
7. Dați exemplu de doi vectori care nu pot fi adunați.

---
### 🏠 Temă
Alegeți un sistem biologic cu feedback, descrieți:
- variabilele de stare;
- tipul de feedback;
- o posibilă intervenție;
- un efect neintuitiv posibil.
"""

# ╔═╡ Cell order:
# ╠═3c523909-557c-4473-9736-4e271c3e958a
# ╠═a8e79b63-4624-4f3c-a985-5bf6a1f51d43
# ╠═ccc7e71a-1ce3-4bf0-92fb-9906287c5176
# ╟─b1c2d3e4-5f6a-7b8c-9d0e-1f2a3b4c5d6e
# ╟─1855e998-fc81-4de3-b833-e903669b6f03
# ╟─215a9db0-9667-11f0-3e8b-db501bdbd2d5
# ╟─8c1d2e3f-4a5b-4c6d-8e7f-9a0b1c2d3e4f
# ╟─11efb0cf-73b4-4d19-b0f7-612a5abec2d3
# ╟─9f282b4f-da9c-49b1-94ef-03cbe39ca1f9
# ╟─4d08da5d-e63f-49d9-a25a-05c9ff0c2cd4
# ╟─9968fa47-4085-4417-8100-7c8385175aca
# ╟─ccd56781-7fbe-48ed-9617-6c88829384c0
# ╟─d42b9ecf-086a-47fe-b2e8-742a3b0ed214
# ╟─474aa1f8-42b2-44d7-a7ef-441d0ddffa36
# ╟─a6e282d9-dd41-440d-a8c6-bc068f23720a
# ╟─391dd0e7-66e6-42cb-a967-0345457f85f7
# ╠═9dd94e6d-04e5-4669-adff-414b25cb1dd3
# ╟─291381ae-569e-4335-9110-b3af9bee2c9f
# ╟─839ddca1-d3ba-4d99-b15f-a795b040c161
# ╟─778d276d-9c3c-49c7-9d44-14baf3ae16e7
# ╟─e4d39336-7ddd-4619-9288-6169aadd6a0e
# ╟─8c24d3d5-bc4f-4f41-a7da-c8d4a4be8077
# ╟─e4d67807-c4cb-4573-a8a9-2b7ba7e59fc9
# ╟─928d4fd2-74ad-48a5-8c83-b24156dca337
# ╟─05ecb7e4-89fc-4265-8b43-0fc6208d7e8e
# ╟─284cd1ac-987b-46dc-a066-e5ac244684e1
# ╟─e4ff1401-fe97-43b5-bb9a-85df5b0310a4
# ╟─6dcfeafd-2cd9-41a0-83c6-d297f1fa2298
# ╟─787799ae-e91c-409b-aba7-58aaa683c0cc
# ╟─21314352-320b-4b2a-a925-2a6ece098140
# ╟─4a964e73-1e3d-4fd8-adf9-ec0830c00b25
# ╟─5477506e-4926-41d0-a31c-7bf486040402
# ╟─019adb7d-8418-4149-a81b-ec99341efe50
# ╟─803b1d6b-e26b-44eb-917d-d0265c228b14
# ╟─4f5a7238-8855-44b1-8afa-b774958caaf0
# ╟─39adb7af-b1f0-4cf0-93ea-be4063654c55
# ╠═17c910ec-7256-4ba0-8ee6-d29c1642703e
# ╠═d8ac2f4c-cc9f-4e12-9275-c429dc65010b
# ╠═1549ca87-3cb5-497b-b812-b66fc8e4d73a
# ╟─7d5f4320-933e-4358-8322-dcf613a8e702
# ╟─920b2efd-f83f-4300-b34f-b9695cf1079b
# ╟─5b1a6088-d669-4892-9ab4-1068dc68ec3f
# ╟─5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
# ╟─4a3e5f6a-7b8c-4d9e-8f0a-1b2c3d4e5f6a
# ╟─9d2e3f4a-5b6c-4d7e-8f9a-0b1c2d3e4f5a
