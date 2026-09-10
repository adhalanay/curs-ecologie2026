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

# ╔═╡ a8e79b63-4624-4f3c-a985-5bf6a1f51d43
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
end

# ╔═╡ 9d43ec2a-5d67-4219-829e-fc747d8d4057
using PlutoUI

# ╔═╡ 3c523909-557c-4473-9736-4e271c3e958a
TableOfContents()

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

# ╔═╡ 215a9db0-9667-11f0-3e8b-db501bdbd2d5
md"""
# 🌿 Cursul I: Sisteme cu feedback și funcții
**andrei.halanay@unibuc.ro**
[📓 Notebook-uri](https://github.com/adhalanay/curs_ecologie2026)

---
### 📝 Cum se va face notarea
- **20 puncte** pentru activitate (răspunsuri, teme etc.);
- **30 puncte** un proiect (poate fi făcut în echipe de maxim 2 persoane);
- **70 puncte** examenul final.

> 💡 *Pentru cei care doresc va fi un examen parțial în săptămîna a 8-a.*
> **Condiții de trecere:** cel puțin 50 de puncte în total și cel puțin 35 de puncte la examen.
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
## 📊 3. Reprezentarea grafică
Graficul obișnuit folosește două axe perpendiculare:
- axa orizontală $Ox$: variabila independentă $x$;
- axa verticală $Oy$: valorile funcției $y=f(x)$.

**Exemplu:**
"""

# ╔═╡ e4d39336-7ddd-4619-9288-6169aadd6a0e
begin
	f(x) = x^5
	xs = range(-1.5, 1.5, length=300)
	ys = f.(xs)
	fig = Figure(size=(600, 400))
	ax = Axis(fig[1,1], 
		xlabel="x", 
		ylabel="f(x)", 
		title="Graficul funcției f(x) = x⁵",
		backgroundcolor=:transparent,
		xgridcolor=:grey40, ygridcolor=:grey40,
		xtickcolor=:white, ytickcolor=:white,
		xlabelcolor=:white, ylabelcolor=:white, titlecolor=:white
	)
	lines!(xs, ys, color="#8be9fd", linewidth=3)
	fig
end

# ╔═╡ 5477506e-4926-41d0-a31c-7bf486040402
md"""
### 📏 Scara logaritmică
Uneori valorile lui $y$ sunt mult mai mari decât cele ale lui $x$. Folosim atunci scara logaritmică pentru una dintre axe.
**Exemple de scări logaritmice:**
- 🧪 pH;
- 🌍 scara Richter;
- 🔊 decibeli.

Într-un grafic semi-logaritmic reprezentăm perechile $(x, \log_{10} y)$.
"""

# ╔═╡ 019adb7d-8418-4149-a81b-ec99341efe50
begin
	xs_slog = range(0.1, 5, length=300)
	ys_slog = f.(xs_slog)
	fig1 = Figure(size=(600, 400))
	ax1 = Axis(fig1[1,1], 
		yscale=log10, 
		xlabel="x", 
		ylabel="log₁₀ f(x)", 
		title="Reprezentare Semi-logaritmică",
		backgroundcolor=:transparent,
		xgridcolor=:grey40, ygridcolor=:grey40,
		xtickcolor=:white, ytickcolor=:white,
		xlabelcolor=:white, ylabelcolor=:white, titlecolor=:white
	)
	lines!(xs_slog, ys_slog, color="#ff79c6", linewidth=3)
	fig1
end

# ╔═╡ 803b1d6b-e26b-44eb-917d-d0265c228b14
md"""
Dacă graficul semi-logaritmic este o dreaptă, datele satisfac o lege exponențială.
Într-un grafic log-log folosim scară logaritmică pentru ambele axe.
"""

# ╔═╡ 4f5a7238-8855-44b1-8afa-b774958caaf0
begin
	fig2 = Figure(size=(600, 400))
	ax2 = Axis(fig2[1,1], 
		xscale=log10, 
		yscale=log10, 
		xlabel="log₁₀ x", 
		ylabel="log₁₀ f(x)", 
		title="Reprezentare Log-Log",
		backgroundcolor=:transparent,
		xgridcolor=:grey40, ygridcolor=:grey40,
		xtickcolor=:white, ytickcolor=:white,
		xlabelcolor=:white, ylabelcolor=:white, titlecolor=:white
	)
	lines!(xs_slog, ys_slog, color="#50fa7b", linewidth=3)
	fig2
end

# ╔═╡ 7d5f4320-933e-4358-8322-dcf613a8e702
md"""
## 🌌 4. Spațiul stărilor
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
## 📏 5. Sisteme de dimensiune mai mare
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
## 📝 6. Recapitulare și exerciții

### 🧠 Recapitulare
- 🔄 **Feedback pozitiv:** variabila își amplifică creșterea.
- 🛑 **Feedback negativ:** variabila își limitează creșterea.
- ⚠️ Sistemele cu feedback pot avea *rebound* și dependență de fază.
- 🎯 O funcție asociază fiecărui input exact un output.
- 🗺️ Domeniul, codomeniul și imaginea sunt noțiuni diferite.
- 🎨 Funcțiile pot fi reprezentate verbal, numeric, vizual și algebric.
- 📏 Scara logaritmică este utilă când valorile diferă mult ca ordin de mărime.
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
# ╟─a8e79b63-4624-4f3c-a985-5bf6a1f51d43
# ╟─9d43ec2a-5d67-4219-829e-fc747d8d4057
# ╟─b1c2d3e4-5f6a-7b8c-9d0e-1f2a3b4c5d6e
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
# ╠═e4d39336-7ddd-4619-9288-6169aadd6a0e
# ╟─5477506e-4926-41d0-a31c-7bf486040402
# ╟─019adb7d-8418-4149-a81b-ec99341efe50
# ╟─803b1d6b-e26b-44eb-917d-d0265c228b14
# ╟─4f5a7238-8855-44b1-8afa-b774958caaf0
# ╟─7d5f4320-933e-4358-8322-dcf613a8e702
# ╟─920b2efd-f83f-4300-b34f-b9695cf1079b
# ╟─5b1a6088-d669-4892-9ab4-1068dc68ec3f
# ╟─5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
# ╟─4a3e5f6a-7b8c-4d9e-8f0a-1b2c3d4e5f6a
# ╟─9d2e3f4a-5b6c-4d7e-8f9a-0b1c2d3e4f5a
