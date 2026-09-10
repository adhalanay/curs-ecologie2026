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

# ╔═╡ 215a9db0-9667-11f0-3e8b-db501bdbd2d5
md"""
# Cursul I

andrei.halanay@unibuc.ro

[Notebook-uri](https://github.com/adhalanay/curs_ecologie2026)

## Cum se va face notarea

- 20 de puncte pentru activitate (răspunsuri, teme etc.);
- 30 de puncte un proiect (poate fi făcut în echipe de maxim 2 persoane);
- 70 de puncte examenul final.

Pentru cei care doresc va fi un examen parțial în săptămîna a 8-a.

Condiții de trecere: cel puțin 50 de puncte în total și cel puțin 35 de puncte la examen.
"""

# ╔═╡ 8c1d2e3f-4a5b-4c6d-8e7f-9a0b1c2d3e4f
md"""
# Obiective


1. Diferența dintre feedback pozitiv și feedback negativ;
2. Exemple de sisteme biologice cu feedback;
3. Noțiunile de funcție, domeniu, codomeniu și imagine;
4. Compunerea a două funcții;
5. Grafice obișnuite, semi-logaritmice și logaritmice;
7. Spațiul stărilor și dimensiunea acestuia;
8. Operații cu vectori
"""

# ╔═╡ 11efb0cf-73b4-4d19-b0f7-612a5abec2d3
md"""
# 1. Sisteme cu feedback

Modelarea matematică a proceselor biologice a început să fie folosită sistematic de la mijlocul secolului XIX. Astăzi matematica este esențială în:

- imunologie și boli autoimune;
- farmacologie și proiectarea medicamentelor;
- dinamica populațiilor;
- epidemiologie;
- neuroștiință;
- imagistică medicală.

Un exemplu clasic de studiu ecologic este oscilația dintre:

- râsul canadian (*Lynx canadensis*);
- iepurele de zăpadă (*Lepus americanus*).

Râsul se hrănește aproape exclusiv cu iepurele de zăpadă.
"""

# ╔═╡ 9f282b4f-da9c-49b1-94ef-03cbe39ca1f9
PlutoUI.LocalResource("Canadalynx.jpg", :width => 300)

# ╔═╡ 4d08da5d-e63f-49d9-a25a-05c9ff0c2cd4
md"""
Iepurele de zăpadă:
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
Pe termen scurt:
"""

# ╔═╡ a6e282d9-dd41-440d-a8c6-bc068f23720a
PlutoUI.LocalResource("rîși.png")

# ╔═╡ 391dd0e7-66e6-42cb-a967-0345457f85f7
md"""
Se observă:

- oscilații cu perioadă de aproximativ 10 ani;
- populația de prădători crește/scade cu întârziere față de cea a prăzii.

Pentru a explica aceste oscilații avem nevoie de un model matematic.

**Ideea centrală:** fiecare populație determină valoarea celeilalte populații.

## Feedback pozitiv

Avem feedback pozitiv dacă o valoare pozitivă a unei variabile determină creșterea acelei variabile.

Exemple:

- banii investiți pot produce mai mulți bani;
- animalele fac pui, ceea ce crește populația, deci numărul de pui;
- creșterea nivelului de $CO_2$ crește temperatura, ceea ce accelerează descompunerea materiei organice și produce mai mult $CO_2$.

## Feedback negativ

Avem feedback negativ dacă o valoare pozitivă a unei variabile determină scăderea ei, iar o valoare negativă determină creșterea ei.

Exemplu: aerul condiționat.

Fie ``T_0`` temperatura setată și ``C`` temperatura curentă. Definim:

```math
T = C - T_0.
```

- dacă ``T > 0``, sistemul răcește;
- dacă ``T < 0``, sistemul încălzește.

Un alt exemplu clasic: glucoza și insulina.

Când glucoza crește, insulina crește, iar glicemia scade.
"""

# ╔═╡ 9dd94e6d-04e5-4669-adff-414b25cb1dd3
PlutoUI.LocalResource("insulină.png")

# ╔═╡ 291381ae-569e-4335-9110-b3af9bee2c9f
md"""
Comportamentul este asemănător cu sistemul râși/iepuri.

În epidemiologie, contactele dintre infectați și susceptibili cresc numărul de infectați și scad numărul de susceptibili. Modelele epidemiologice ajută la elaborarea strategiilor de intervenție.

## Comportamente neintuitive

Sistemele cu feedback au adesea efecte neintuitive.

Exemplu: vrem să reducem numărul de prădători. Dacă scoatem prădători din mediu, prada crește, ceea ce poate duce la creșterea numărului de prădători peste valoarea inițială. Acesta este fenomenul de **rebound**.

Rezultatul unei intervenții depinde de faza ciclului în care are loc intervenția.

**Concluzie:** sisteme aparent simple pot avea comportamente foarte neintuitive. De aceea modelarea matematică este indispensabilă.

---

# 2. Funcții

Un concept fundamental este cel de **funcție**.

O funcție este o relație între o mulțime de date de intrare și una de date de ieșire astfel încât fiecărui input îi corespunde **exact** un output.

Funcțiile pot fi definite:

- printr-un tabel cu două coloane;
- printr-o formulă;
- printr-un grafic;
- verbal.

Exemplu:

```math
f(x) = x^2.
```

Două aspecte importante:

1. Nu orice formulă definește o funcție. De exemplu, ``y^2 = x^2`` nu definește o funcție deoarece pentru $x=2$ avem $y=2$ și $y=-2$.
2. Nu orice serie temporală poate fi descrisă printr-o formulă explicită.

## Domeniu, codomeniu, imagine

- **Domeniul** este mulțimea valorilor de intrare.
- **Codomeniul** este mulțimea în care funcția ia valori.
- **Imaginea** este mulțimea valorilor efectiv atinse.

Pentru ``f(x)=x^2``:

- domeniul poate fi ``\mathbb{R}``;
- codomeniul poate fi ``\mathbb{R}``;
- imaginea este ``\mathbb{R}_+ = [0,\infty)``.

Notăm o funcție cu domeniu $X$ și codomeniu $Y$ astfel:

```math
f : X \to Y.
```

## Exemplu din biologia moleculară

ADN-ul este compus din bazele ``A, C, G, T``. La transcriere, ``T`` este înlocuit cu ``U`` în ARN mesager.

Avem funcția:

```math
\text{transcriere} : \{A,C,G,T\} \to \{A,C,G,U\}.
```

Transcrierea duce un codon ADN într-un codon ARN. Apoi translația duce un codon ARN într-un aminoacid:

```math
\text{translație} : \text{codoni ARN} \to \text{aminoacizi}.
```

Expresia genetică este compunerea celor două:

```math
\text{expresie genetică} : \text{codoni ADN} \to \text{aminoacizi}.
```

## Compunerea funcțiilor

Dacă ``f(x)=2x^2+1`` și ``g(x)=\sqrt{x}``, atunci:

```math
(g \circ f)(x) = g(f(x)) = \sqrt{2x^2+1},
```

```math
(f \circ g)(x) = f(g(x)) = 2x+1.
```

**Exercițiu.** Presupunem că pe Marte codonul ``AUC`` produce în 60% din cazuri izoleucină și în 40% treonină. Mai este expresia genetică o funcție?

*Răspuns:* Nu, deoarece același input poate da două output-uri diferite.

## Patru reprezentări ale unei funcții

1. verbal;
2. numeric, printr-un tabel;
3. vizual, printr-un grafic;
4. algebric, printr-o formulă.
"""

# ╔═╡ 839ddca1-d3ba-4d99-b15f-a795b040c161
md"""
# 3. Reprezentarea grafică

Graficul obișnuit folosește două axe perpendiculare:

- axa orizontală $Ox$: variabila independentă $x$;
- axa verticală $Oy$: valorile funcției $y=f(x)$.

Exemplu:
"""

# ╔═╡ e4d39336-7ddd-4619-9288-6169aadd6a0e
begin
	f(x) = x^5
	xs = range(-1.5, 1.5, length=300)
	ys = f.(xs)

	fig = Figure()
	ax = Axis(fig[1,1], xlabel="x", ylabel="f(x)")
	lines!(xs, ys)
	fig
end

# ╔═╡ 5477506e-4926-41d0-a31c-7bf486040402
md"""
## Scara logaritmică

Uneori valorile lui $y$ sunt mult mai mari decât cele ale lui $x$. Folosim atunci scara logaritmică.

Exemple de scări logaritmice:

- pH;
- scara Richter;
- decibeli.

Într-un grafic **semi-logaritmic** reprezentăm perechile $(x, \log_{10} y)$.
"""

# ╔═╡ 019adb7d-8418-4149-a81b-ec99341efe50
begin
	xs_slog = range(0.1, 5, length=300)
	ys_slog = f.(xs_slog)

	fig1 = Figure()
	ax1 = Axis(fig1[1,1], yscale=log10, xlabel="x", ylabel="log10 f(x)")
	lines!(xs_slog, ys_slog)
	fig1
end

# ╔═╡ 803b1d6b-e26b-44eb-917d-d0265c228b14
md"""
Dacă graficul semi-logaritmic este o dreaptă, datele satisfac o lege exponențială.

Într-un grafic **log-log** folosim scară logaritmică pentru ambele axe.
"""

# ╔═╡ 4f5a7238-8855-44b1-8afa-b774958caaf0
begin
	
	fig2 = Figure()
	ax2 = Axis(fig2[1,1], xscale=log10, yscale=log10, xlabel="log10 x", ylabel="log10 f(x)")
	lines!(xs_slog, ys_slog)
	fig2
end

# ╔═╡ 7d5f4320-933e-4358-8322-dcf613a8e702
md"""
# 4. Spațiul stărilor

Variabilele de stare descriu cantitativ sistemul la un moment dat.

Exemplu: starea unei populații poate fi mărimea ei. Putem fi interesați și de:

- raportul între sexe;
- distribuția pe clase de vârstă;
- distribuția teritorială.

Alegerea variabilelor de stare este cea mai importantă și dificilă parte a modelării.

**Variabilele de stare sunt funcții de timp.**

Mulțimea tuturor valorilor posibile ale variabilelor de stare se numește **spațiul stărilor**.

Presupunem de obicei:

- spațiul stărilor este o mulțime de numere;
- variabilele pot lua valori intermediare între anumite limite: ipoteza de continuitate.

## Sisteme de dimensiune 1

Dacă avem o singură variabilă de stare, spațiul stărilor este mulțimea numerelor reale sau un interval.

Exemple:

- voltaj: $(-\infty,\infty)$;
- mărimea populației: $[0,\infty)$;
- temperatură în Kelvin: $[0,\infty)$;
- procentul unui grup: $[0,1]$.

**Exercițiu.** Care este spațiul stărilor pentru temperatura măsurată în $^\circ C$?

*Răspuns:* fizic, $[-273.15,\infty)$; matematic, putem lua $(-273.15,\infty)$ sau chiar $\mathbb{R}$ dacă ignorăm limita fizică.

O variabilă de stare desenată cu timpul pe axa orizontală și valoarea pe axa verticală se numește **serie temporală**.
"""

# ╔═╡ 920b2efd-f83f-4300-b34f-b9695cf1079b
md"""
# 5. Sisteme de dimensiune mai mare

Dacă avem două variabile de stare, de exemplu râși ``R`` și iepuri ``I``, spațiul stărilor este format din perechi:

```math
(R, I).
```

Putem defini:

- adunarea pe componente:

```math
(R_1,I_1)+(R_2,I_2)=(R_1+R_2,I_1+I_2);
```

- înmulțirea cu scalari:

```math
a(R,I)=(aR,aI).
```

Putem desena câte o axă pentru fiecare componentă.
"""

# ╔═╡ 5b1a6088-d669-4892-9ab4-1068dc68ec3f
begin
	fig3 = Figure()
	ax3 = Axis(fig3[1,1], xlabel="R", ylabel="I")
	scatter!([2,1,3,5], [3,2,4,1])
	fig3
end

# ╔═╡ 5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
md"""
O pereche de numere se numește **2-vector**. Fiecare număr este o **componentă**.

Un vector poate fi reprezentat și ca o săgeată:
"""

# ╔═╡ 4a3e5f6a-7b8c-4d9e-8f0a-1b2c3d4e5f6a
begin
	fig4 = Figure()
	ax4 = Axis(fig4[1,1], xlabel="x", ylabel="y")
	arrows2d!([0], [0], [1.985], [2.985])
	scatter!(1.985, 2.985, markersize=7)
	fig4
end

# ╔═╡ 9d2e3f4a-5b6c-4d7e-8f9a-0b1c2d3e4f5a
md"""
Înmulțirea cu scalari are interpretare geometrică:

- dacă ``|a|<1``, vectorul se scurtează;
- dacă ``|a|>1``, vectorul se lungește;
- dacă ``a<0``, vectorul își schimbă sensul.

Pentru mai multe variabile folosim vectori cu ``n`` componente.

Dacă:

```math
\mathbf{a}=(a_1,\dots,a_n), \quad \mathbf{b}=(b_1,\dots,b_n),
```

atunci:

```math
\mathbf{a}+\mathbf{b}=(a_1+b_1,\dots,a_n+b_n).
```

**Ca să adunăm doi vectori, ei trebuie să aibă aceeași dimensiune.**

Dacă ``\mathbf{a}=(a_1,\dots,a_n)``, atunci:

```math
\alpha \mathbf{a}=(\alpha a_1,\dots,\alpha a_n).
```

Mulțimea ``n``-vectorilor se notează cu ``\mathbb{R}^n``.

---

# 6. Recapitulare și exerciții

## Recapitulare

1. Feedback pozitiv: variabila își amplifică creșterea.
2. Feedback negativ: variabila își limitează creșterea.
3. Sistemele cu feedback pot avea rebound și dependență de fază.
4. O funcție asociază fiecărui input exact un output.
5. Domeniul, codomeniul și imaginea sunt noțiuni diferite.
6. Funcțiile pot fi reprezentate verbal, numeric, vizual și algebric.
7. Scara logaritmică este utilă când valorile diferă mult ca ordin de mărime.
8. Variabilele de stare sunt funcții de timp.
9. Spațiul stărilor este mulțimea valorilor posibile ale variabilelor de stare.
10. Pentru mai multe variabile folosim vectori în $\mathbb{R}^n$.

## Exerciții

1. Dați un exemplu de feedback pozitiv și unul de feedback negativ din biologie.
2. Explicați fenomenul de rebound într-un sistem prădător-pradă.
3. Pentru $f(x)=\sqrt{x}$ și $g(x)=x^2+1$, calculați ``f\circ g`` și ``g\circ f``. Precizați domeniile.
4. Desenați graficul lui ``f(x)=2^x`` și apoi graficul semi-logaritmic. Ce observați?
5. Care este spațiul stărilor pentru:
   - masa unui organism;
   - concentrația unei substanțe;
   - un procent?
6. Calculați:
   ```math
   (2,3)+(-1,4), \quad 3(2,-1).
   ```
7. Dați exemplu de doi vectori care nu pot fi adunați.

## Temă

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
# ╠═4f5a7238-8855-44b1-8afa-b774958caaf0
# ╟─7d5f4320-933e-4358-8322-dcf613a8e702
# ╟─920b2efd-f83f-4300-b34f-b9695cf1079b
# ╠═5b1a6088-d669-4892-9ab4-1068dc68ec3f
# ╟─5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
# ╠═4a3e5f6a-7b8c-4d9e-8f0a-1b2c3d4e5f6a
# ╟─9d2e3f4a-5b6c-4d7e-8f9a-0b1c2d3e4f5a
