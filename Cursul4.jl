### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> title = "Cursul IV"
#> description = "Sisteme dinamice discrete și haos"
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

# ╔═╡ e7e9e376-adb0-11f1-91f0-f7cc7081a2ec
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using PlutoTeachingTools
	using Printf

	# Setăm o temă vizuală curată și modernă pentru grafice
	CairoMakie.activate!(type = "svg")
	set_theme!(theme_minimal())

	TableOfContents()
end

# ╔═╡ ce44eb4e-4b93-43fa-9105-1a10e9015ccd
WidthOverDocs()

# ╔═╡ f612f2c6-5756-406d-ae1d-f4d55b03e87f
begin
	# cele trei modele din curs, cu capacitatea de suport normalizată la K = 1
	f_logistic(x, r) = r * x * (1 - x)
	f_ricker(x, r) = x * exp(r * (1 - x))
	f_beverton(x, r) = r * x / (1 + (r - 1) * x)

	# primii n pași ai sistemului x_{k+1} = f(x_k), pornind din x₀
	function itereaza(f, x0, n)
		xs = [float(x0)]
		for k in 1:n
			push!(xs, f(xs[end]))
		end
		xs
	end

	# diagrama pînză de păianjen: graficul lui f, prima bisectoare și n pași ai iterației
	function panza!(ax, f, x0, n, xmin, xmax; culoare = :darkorange)
		xs = range(xmin, xmax, length = 400)
		lines!(ax, xs, xs, color = :gray, linestyle = :dash, linewidth = 1.5)
		lines!(ax, xs, f.(xs), color = :teal, linewidth = 3)
		traj = itereaza(f, x0, n)
		px, py = [traj[1]], [0.0]                       # pornim de pe axa Ox
		for k in 1:n
			push!(px, traj[k], traj[k+1])               # vertical pînă la grafic,
			push!(py, traj[k+1], traj[k+1])             # apoi orizontal pînă la bisectoare
		end
		lines!(ax, px, py, color = (culoare, 0.8), linewidth = 1.5)
		scatter!(ax, [x0], [0.0], color = culoare, markersize = 12)
		xlims!(ax, xmin, xmax)
		traj
	end
	nothing
end

# ╔═╡ 5a49c8fe-6d17-4ce3-8865-3c0ba68af9c8
md"""
# 🌿 Cursul 4 — Sisteme dinamice discrete și haos

**Durata estimată:** 100 minute
"""

# ╔═╡ e07e94e7-938f-4546-b9cc-71f297ba1234
md"""
## 🎯 Obiectivele lecției
1. Să modelăm populațiile cu **generații care nu se suprapun** prin ecuații cu diferențe ``x_{n+1} = f(x_n)``.
2. Să rezolvăm modelul de **creștere geometrică** și să îl legăm de creșterea exponențială și de metoda Euler.
3. Să găsim **punctele de echilibru** ``x^* = f(x^*)`` și să le determinăm stabilitatea, grafic (**diagrama pînză de păianjen**) și analitic (**criteriul ``\lvert f'(x^*)\rvert < 1``**).
4. Să urmărim, în **modelul logistic discret**, drumul de la echilibru la **cicluri periodice** și la **haos** (diagrama de bifurcație, exponentul Liapunov).
5. Să comparăm modelele **Ricker** și **Beverton–Holt**, folosite în gestionarea pescuitului.
"""

# ╔═╡ 9f47681d-543a-40f2-a663-cc685f395057
md"""
## 🔁 0. Ce știm deja
În Cursul 3 am studiat ecuații în **timp continuu**, ``x' = f(x)``:
- echilibrele sînt soluțiile ecuației ``f(x^*) = 0``;
- un echilibru este stabil dacă ``f'(x^*) < 0`` și instabil dacă ``f'(x^*) > 0``;
- soluțiile sînt **monotone**: nu pot traversa un echilibru, deci nu pot oscila.

Multe specii însă **nu se reproduc continuu**, ci în pulsuri, o dată pe sezon:

> 💡 **Exemple din natură:**
> * 🐟 **Somonul de Pacific** (specie semelpară): se reproduce o singură dată în viață, apoi moare.
> * 🦗 **Magicicada**: trăiește 13 sau 17 ani sub pămînt ca nimfă, apoi toți indivizii ies în același an pentru a se reproduce.
> * 🌸 **Planta secolului** (*Agave americana*): înflorește după zeci de ani, produce semințe și apoi moare.
> * 🦋 Multe **insecte din zona temperată**: o singură generație pe an, adulții mor înainte să apară urmașii.

Pentru astfel de specii este natural să numărăm populația o dată pe generație: ``x_0, x_1, x_2, \dots``. Timpul nu mai „curge”, ci **sare** din generație în generație, iar modelul devine o **ecuație cu diferențe** (sau **sistem dinamic discret**):
```math
x_{n+1} = f(x_n), \qquad n = 0, 1, 2, \dots
```
Funcția ``f`` spune cîți indivizi vor fi în generația următoare, știind cîți sînt acum. Cunoscînd ``x_0``, întregul viitor este determinat:
```math
x_1 = f(x_0), \quad x_2 = f(f(x_0)), \quad \dots, \quad x_n = \underbrace{f(f(\cdots f}_{n\ \text{ori}}(x_0))) = f^{n}(x_0).
```
Atenție: aici ``f^{n}`` înseamnă **compunerea** lui ``f`` cu ea însăși de ``n`` ori, nu puterea ``n``.

Am mai întîlnit un astfel de sistem: **metoda Euler** ``x_{n+1} = x_n + \Delta t\, f(x_n)`` este o ecuație cu diferențe.

> 🎬 **Surpriza cursului:** ecuația logistică, care în timp continuu are un comportament liniștit (orice soluție tinde monoton la ``K``), poate produce în timp discret **oscilații**, **cicluri** și **haos**.
"""

# ╔═╡ 1e558223-6779-4a42-b75c-39e7c1cfb046
md"""
## 📈 1. Creșterea geometrică

**Cel mai simplu model:** o populație de căprioare crește cu ``5\%`` pe an. Dacă ``x_n`` este numărul de căprioare din anul ``n``, atunci
```math
x_{n+1} = x_n + 0.05\,x_n = 1.05\,x_n .
```
În general, dacă fiecare individ lasă în medie ``\lambda`` urmași (incluzîndu-l pe el însuși, dacă supraviețuiește), obținem modelul de **creștere geometrică**
```math
x_{n+1} = \lambda\, x_n \quad\Longrightarrow\quad x_n = \lambda^n x_0 .
```
- ``\lambda < 1``: populația scade geometric către ``0`` (extincție);
- ``\lambda = 1``: populația rămîne constantă;
- ``\lambda > 1``: populația crește nelimitat.

**Legătura cu timpul continuu.** Scriind ``\lambda = e^{r}``, adică ``r = \ln \lambda``, obținem ``x_n = x_0 e^{rn}``: exact soluția ecuației ``x' = rx`` din Cursul 3, citită la momentele ``t = 0, 1, 2, \dots`` Timpul de dublare este deci
```math
T_2 = \frac{\ln 2}{\ln \lambda}.
```
Pe de altă parte, metoda Euler cu pasul ``\Delta t`` aplicată ecuației ``x' = rx`` dă ``x_{n+1} = (1 + r\Delta t)\,x_n``, tot o creștere geometrică, cu ``\lambda = 1 + r\Delta t``.

### 🎛️ Experimentați
Modificați rata de creștere ``\lambda`` și populația inițială ``x_0`` și observați evoluția populației în 50 de ani.
"""

# ╔═╡ aa081a08-4f61-48a9-865b-f3e487b5d231
md"""
Rata de creștere ``\lambda``: $(@bind λ_geo PlutoUI.Slider(0.80:0.01:1.20, show_value=true, default=1.05))

Populația inițială ``x_0``: $(@bind x0_geo PlutoUI.Slider(0.5:0.5:5.0, show_value=true, default=1.0))
"""

# ╔═╡ 09fd950e-89f8-4627-adfd-ae595cacc545
let
	xs = itereaza(x -> λ_geo * x, x0_geo, 50)
	subtitlu = λ_geo > 1 ? @sprintf("timp de dublare T₂ = ln 2 / ln λ ≈ %.1f ani", log(2) / log(λ_geo)) :
		λ_geo < 1 ? @sprintf("timp de înjumătățire ln 2 / |ln λ| ≈ %.1f ani", log(2) / abs(log(λ_geo))) :
		"populație constantă"

	fig = Figure(size = (900, 400))
	ax = Axis(fig[1, 1],
		title = "Creștere geometrică: xₙ = λⁿ x₀,  λ = $(λ_geo)",
		subtitle = subtitlu,
		xlabel = "Timp (ani)",
		ylabel = "Populație (xₙ)",
		limits = (0, 50, 0, nothing)
	)
	lines!(ax, 0:50, xs, color = :teal, linewidth = 2)
	scatter!(ax, 0:50, xs, color = :teal, markersize = 8, strokewidth = 1, strokecolor = :white)
	fig
end

# ╔═╡ 182c62ba-84ab-4233-b7d1-12a534582095
md"""
### 📌 Concluzii parțiale
* Creșterea geometrică este varianta discretă a creșterii exponențiale.
* Creșterea absolută ``x_{n+1} - x_n = (\lambda - 1)\,x_n`` este **proporțională cu populația**: cu cît populația este mai mare, cu atît crește mai mult.
* Modelul este realist doar pe termen scurt: nicio populație nu crește la nesfîrșit.

### ✅ Verificare rapidă
1. Dacă dublăm populația inițială ``x_0`` (de la 1 la 2), ce se întîmplă cu graficul pentru ``\lambda = 1.1``? Verificați cu slider-ele.
2. O populație de căprioare crește cu ``5\%`` pe an. În cîți ani se dublează?
"""

# ╔═╡ 2c8f0921-8e71-49a9-b910-9940abcd25bd
hint(md"1. Toate valorile se dublează, deoarece ``x_n = \lambda^n x_0`` este proporțional cu ``x_0`` (modelul este **liniar**); forma graficului nu se schimbă. 2. ``T_2 = \ln 2 / \ln 1.05 \approx 0.693 / 0.0488 \approx 14.2`` ani.")

# ╔═╡ dabc61af-b537-48ac-9e20-caca44da4b65
md"""
## 🌾 2. Modelul logistic discret

Ca în timp continuu, presupunem că mediul poate susține cel mult ``K`` indivizi (**capacitatea de suport**) și că rata de reproducere scade cînd resursele se împuținează:
* o populație ``x_n`` consumă fracțiunea ``\frac{x_n}{K}`` din resurse;
* rămîne liberă fracțiunea ``1 - \frac{x_n}{K}``;
* fiecare individ lasă un număr de urmași proporțional cu resursele libere, cu factorul ``r``.

Obținem funcția
```math
 f(x) = r\, x \left(1 - \frac{x}{K}\right) .
```
Măsurînd populația ca fracțiune din capacitatea de suport (adică normalizînd ``K = 1``), obținem celebrul **model logistic discret**, popularizat de ecologul **Robert May** în 1976:
```math
x_{n+1} = r\, x_n (1 - x_n), \qquad 0 \le x_n \le 1 .
```
Graficul lui ``f`` este o parabolă cu maximul ``f(1/2) = r/4``. Ca populația să rămînă în intervalul ``[0, 1]``, trebuie ca ``0 \le r \le 4``.

### 🎛️ Experimentați
Mutați slider-ul pentru ``r`` și urmăriți ce se întîmplă. Încercați, pe rînd, ``r = 0.8,\ 1.8,\ 2.8,\ 3.2,\ 3.5,\ 3.9``.
"""

# ╔═╡ 2a3b4bbb-7891-455c-99e0-d820d1b4ffe3
md"""
Valoarea inițială ``x_0``: $(@bind x0_log PlutoUI.Slider(0.01:0.01:0.99, show_value=true, default=0.1))

Parametrul ``r``: $(@bind r_log PlutoUI.Slider(0.1:0.01:4.0, show_value=true, default=2.8))
"""

# ╔═╡ b46b9ded-0cd2-4fec-a51c-17baaaf52db9
let
	n = 100
	xs = itereaza(x -> f_logistic(x, r_log), x0_log, n)

	fig = Figure(size = (900, 400))
	ax = Axis(fig[1, 1],
		title = "Modelul logistic discret: xₙ₊₁ = r xₙ (1 − xₙ),  r = $(r_log)",
		xlabel = "Generația (n)",
		ylabel = "Populație normalizată (xₙ)",
		limits = (0, n, 0, 1.05)
	)
	if r_log > 1
		hlines!(ax, [1 - 1 / r_log], color = :purple, linestyle = :dash, linewidth = 2,
			label = "x* = 1 − 1/r")
		axislegend(ax, position = :rt)
	end
	lines!(ax, 0:n, xs, color = :darkorange, linewidth = 1.5)
	scatter!(ax, 0:n, xs, color = :darkorange, markersize = 6)
	fig
end

# ╔═╡ 4ce78424-8bc8-4cc4-986b-3ad350e627e1
md"""
### 🤔 Ce am observat?
* Pentru ``r < 1`` populația dispare, iar pentru ``1 < r < 3`` ea se stabilizează la nivelul ``x^* = 1 - \frac{1}{r}`` (linia mov).
* Pentru ``2 < r < 3`` populația **oscilează** în jurul echilibrului înainte de a se stabiliza. Pentru ecuația logistică continuă acest lucru era imposibil!
* Pentru ``r`` puțin mai mare decît ``3`` populația alternează între **două valori**, iar pentru ``r`` aproape de ``4`` pare să se comporte **la întîmplare**.

**De ce oscilează?** În timp continuu populația „simte” imediat cînd se apropie de ``K`` și încetinește. În timp discret, întreaga generație următoare este stabilită de densitatea din generația curentă: reacția vine **cu o întîrziere de o generație**. Dacă ``r`` este mare, o populație mică produce prea mulți urmași, depășește capacitatea mediului, iar generația următoare se prăbușește (**supracompensare**).

> 🧮 Același fenomen apare cînd aplicăm metoda Euler ecuației logistice continue cu un pas ``\Delta t`` prea mare: soluția numerică începe să oscileze, deși soluția exactă este monotonă (Exercițiul 3).

În secțiunile următoare explicăm matematic toate aceste comportamente.
"""

# ╔═╡ 01fa65ac-bed0-41de-b98c-1c1ca5ea5622
md"""
## 🕸️ 3. Echilibre, diagrama pînză de păianjen și stabilitate

### 3.1. Puncte de echilibru
Un **punct de echilibru** (punct fix) este o stare care se reproduce identic de la o generație la alta:
```math
x^* = f(x^*) .
```
Grafic, echilibrele sînt punctele în care graficul lui ``f`` intersectează **prima bisectoare** ``y = x``.

Pentru modelul logistic: ``x = rx(1-x) \iff x\big(r(1-x) - 1\big) = 0``, deci
```math
x^* = 0 \qquad \text{și} \qquad x^* = 1 - \frac{1}{r} = \frac{r-1}{r} \quad (\text{în } (0,1) \text{ doar dacă } r > 1).
```

### 3.2. Diagrama pînză de păianjen
Iterația se poate urmări direct pe graficul lui ``f``:
1. pornim din ``x_0`` pe axa ``Ox`` și urcăm **vertical** pînă la graficul lui ``f``: ajungem în punctul ``(x_0, x_1)``;
2. mergem **orizontal** pînă la bisectoarea ``y = x``: ajungem în ``(x_1, x_1)``, adică am „mutat” ``x_1`` pe axa ``Ox``;
3. repetăm: vertical pînă la grafic, orizontal pînă la bisectoare, ...

Linia frîntă obținută seamănă cu o pînză de păianjen (în engleză, *cobweb diagram*).
"""

# ╔═╡ 7128b3da-1907-4946-a117-41cbfa3381e5
md"""
Parametrul ``r``: $(@bind r_pz PlutoUI.Slider(0.5:0.01:4.0, show_value=true, default=2.8))

Valoarea inițială ``x_0``: $(@bind x0_pz PlutoUI.Slider(0.01:0.01:0.99, show_value=true, default=0.1))

Numărul de pași ``n``: $(@bind n_pz PlutoUI.Slider(1:1:100, show_value=true, default=15))
"""

# ╔═╡ 37806866-5fb1-4f6e-8d76-0133ee350ebb
let
	f(x) = f_logistic(x, r_pz)
	fig = Figure(size = (1100, 450))
	subtitlu = r_pz > 1 ?
		@sprintf("x* = 1 − 1/r = %.3f,   f'(x*) = 2 − r = %.2f", 1 - 1 / r_pz, 2 - r_pz) :
		@sprintf("x* = 0,   f'(0) = r = %.2f", r_pz)
	ax1 = Axis(fig[1, 1], title = "Diagrama pînză de păianjen", subtitle = subtitlu,
		xlabel = "xₙ", ylabel = "xₙ₊₁", aspect = 1)
	traj = panza!(ax1, f, x0_pz, n_pz, 0, 1)
	ylims!(ax1, 0, 1)
	if r_pz > 1
		scatter!(ax1, [1 - 1 / r_pz], [1 - 1 / r_pz], color = :purple, markersize = 14)
	end

	ax2 = Axis(fig[1, 2], title = "Evoluția în timp", xlabel = "Generația (n)", ylabel = "xₙ",
		limits = (0, n_pz, 0, 1))
	r_pz > 1 && hlines!(ax2, [1 - 1 / r_pz], color = :purple, linestyle = :dash)
	scatterlines!(ax2, 0:n_pz, traj, color = :darkorange, markersize = 7)
	fig
end

# ╔═╡ ef07c739-7314-48ac-bad0-0206cf48f2fe
md"""
### 3.3. Criteriul de stabilitate
Ca în Cursul 3, un echilibru este **stabil** (atractor) dacă traiectoriile care pornesc aproape de el se apropie de el și **instabil** (repulsor) dacă se îndepărtează.

Fie ``u_n = x_n - x^*`` abaterea de la echilibru. Aproximînd ``f`` prin tangenta sa în ``x^*`` (Cursul 2):
```math
u_{n+1} = f(x^* + u_n) - x^* \approx \underbrace{f(x^*) - x^*}_{=\,0} + f'(x^*)\,u_n .
```
Abaterea satisface deci o ecuație de **creștere geometrică** (secțiunea 1) cu factorul ``\lambda = f'(x^*)``, deci ``u_n \approx \big(f'(x^*)\big)^n u_0``.

> **Teorema de stabilitate.** Fie ``x^*`` un echilibru al sistemului ``x_{n+1} = f(x_n)``.
> * Dacă ``\lvert f'(x^*)\rvert < 1``, atunci ``x^*`` este **stabil**.
> * Dacă ``\lvert f'(x^*)\rvert > 1``, atunci ``x^*`` este **instabil**.
> * Dacă ``\lvert f'(x^*)\rvert = 1``, criteriul **nu decide**.

**Semnul** derivatei arată cum se apropie (sau se depărtează) traiectoria:
* ``f'(x^*) > 0``: abaterea își păstrează semnul ⟹ evoluție **monotonă** (pînza de păianjen arată ca o scară);
* ``f'(x^*) < 0``: abaterea își schimbă semnul la fiecare pas ⟹ **oscilații** în jurul echilibrului (pînza de păianjen este o spirală).
"""

# ╔═╡ c1b2337e-73a6-4842-8e62-7bfc97e8b519
let
	fig = Figure(size = (1100, 420))
	cazuri = [(0.7, "stabil, monoton"), (-0.7, "stabil, cu oscilații"),
		(1.15, "instabil, monoton"), (-1.15, "instabil, cu oscilații")]
	for (i, (λ, descriere)) in enumerate(cazuri)
		ax = Axis(fig[1, i], title = "f'(x*) = $λ", subtitle = descriere,
			xlabel = "n", ylabel = i == 1 ? "abaterea uₙ" : "", limits = (-0.5, 15.5, -3.5, 3.5))
		hlines!(ax, 0, color = :purple, linestyle = :dash)
		culoare = abs(λ) < 1 ? :seagreen : :crimson
		scatterlines!(ax, 0:15, λ .^ (0:15), color = culoare, markersize = 7)
	end
	fig
end

# ╔═╡ bc222929-5b58-48b3-8b1a-00aa502ed3e5
md"""
> ⚠️ **Diferența față de timpul continuu.** Pentru ``x' = f(x)`` condiția de stabilitate era ``f'(x^*) < 0``: oricît de negativă ar fi derivata, echilibrul rămîne stabil. În timp discret o derivată **prea negativă** (``f'(x^*) < -1``) destabilizează echilibrul: corecția este atît de puternică încît sistemul „sare” de cealaltă parte, tot mai departe. Aceasta este supracompensarea din secțiunea 2.

### 3.4. Aplicație: modelul logistic
Avem ``f(x) = rx(1-x)``, deci ``f'(x) = r(1 - 2x)``.
* În ``x^* = 0``: ``f'(0) = r``. Extincția este stabilă pentru ``r < 1`` și instabilă pentru ``r > 1``.
* În ``x^* = 1 - \frac{1}{r}``: ``f'(x^*) = r\left(1 - 2 + \frac{2}{r}\right) = 2 - r``. Deci:

| Interval pentru ``r`` | ``f'(x^*) = 2 - r`` | Comportament lîngă ``x^* = 1 - 1/r`` |
| :--- | :---: | :--- |
| ``1 < r < 2`` | ``(0, 1)`` | 🎯 stabil, apropiere monotonă |
| ``2 < r < 3`` | ``(-1, 0)`` | 🎯 stabil, oscilații amortizate |
| ``r > 3`` | ``< -1`` | ❌ instabil: ce se întîmplă atunci? |

Verificați tabelul cu diagrama pînză de păianjen de mai sus!

### ✅ Verificare rapidă
Pentru ``x_{n+1} = 2.5\,x_n(1 - x_n)``, găsiți echilibrul nenul și stabilitatea lui. Apropierea de echilibru este monotonă sau cu oscilații?
"""

# ╔═╡ 3b9199fc-11d0-4908-beaa-9da23c27ffe8
hint(md"``x^* = 1 - 1/2.5 = 0.6`` și ``f'(0.6) = 2.5\,(1 - 1.2) = -0.5``. Cum ``\lvert -0.5\rvert < 1``, echilibrul este stabil; derivata fiind negativă, apropierea se face cu oscilații, iar abaterea se înjumătățește (în valoare absolută) la fiecare generație.")

# ╔═╡ a815fcaa-7fc1-4d0b-997e-6e6fe63158dd
md"""
## 🔄 4. Ciclul de perioadă 2

Pentru ``r > 3`` ambele echilibre sînt instabile, dar populația rămîne în ``[0, 1]``. Simulările arată că ea ajunge să **alterneze între două valori** ``p`` și ``q``: un an bun, un an slab, un an bun, ... Un astfel de comportament se observă la unele populații de insecte și de rozătoare.

Dacă ``f(p) = q`` și ``f(q) = p``, atunci ``f(f(p)) = p``: punctele ``p`` și ``q`` sînt **echilibre ale funcției compuse** ``f^2 = f \circ f``, adică ale sistemului citit din două în două generații. Ecuația ``f(f(x)) = x`` este de gradul 4. Două rădăcini le cunoaștem deja (``0`` și ``1 - 1/r``, care sînt echilibre și pentru ``f^2``); împărțind la ele rămîne o ecuație de gradul 2, cu soluțiile
```math
p,\,q = \frac{(r+1) \pm \sqrt{(r+1)(r-3)}}{2r},
```
reale și distincte **exact cînd ``r > 3``**: ciclul de perioadă 2 se naște în momentul în care echilibrul își pierde stabilitatea.

**Stabilitatea ciclului.** Aplicăm criteriul din secțiunea 3 funcției ``f^2``. Cu regula de derivare a funcțiilor compuse (Cursul 2):
```math
(f^2)'(p) = f'(f(p))\cdot f'(p) = f'(q)\, f'(p) = -r^2 + 2r + 4 .
```
Condiția ``\lvert -r^2 + 2r + 4\rvert < 1`` dă ``3 < r < 1 + \sqrt{6} \approx 3.449``. Pentru ``r > 1 + \sqrt 6`` și ciclul de perioadă 2 devine instabil, iar în locul lui apare un ciclu de **perioadă 4** (echilibre ale lui ``f^4``), apoi unul de perioadă 8, ... Fenomenul se numește **dublarea perioadei**.

Parametrul ``r``: $(@bind r_p2 PlutoUI.Slider(2.5:0.01:4.0, show_value=true, default=3.2))
"""

# ╔═╡ 4f5b0ffe-0db5-485e-94fc-6949b6416e62
let
	f(x) = f_logistic(x, r_p2)
	xs = range(0, 1, length = 400)
	xe = 1 - 1 / r_p2
	fig = Figure(size = (1100, 470))

	subtitlu = r_p2 > 3 ?
		@sprintf("f'(p)·f'(q) = −r² + 2r + 4 = %.3f  ⟹  ciclu %s", -r_p2^2 + 2r_p2 + 4,
			abs(-r_p2^2 + 2r_p2 + 4) < 1 ? "stabil" : "instabil") :
		"r ≤ 3: nu există ciclu de perioadă 2"
	ax1 = Axis(fig[1, 1], title = "Graficele lui f și f∘f", subtitle = subtitlu,
		xlabel = "x", ylabel = "y", aspect = 1, limits = (0, 1, 0, 1))
	lines!(ax1, xs, xs, color = :gray, linestyle = :dash, linewidth = 1.5, label = "y = x")
	lines!(ax1, xs, f.(xs), color = :teal, linewidth = 2.5, label = "f(x)")
	lines!(ax1, xs, f.(f.(xs)), color = :purple, linewidth = 2.5, label = "f(f(x))")
	scatter!(ax1, [xe], [xe], color = :black, markersize = 12, label = "x* = 1 − 1/r")
	if r_p2 > 3
		Δ = sqrt((r_p2 + 1) * (r_p2 - 3))
		pq = [(r_p2 + 1 + Δ) / (2r_p2), (r_p2 + 1 - Δ) / (2r_p2)]
		scatter!(ax1, pq, pq, color = :crimson, markersize = 14, label = "p, q")
	end
	axislegend(ax1, position = :cb)

	ax2 = Axis(fig[1, 2], title = "Pînza de păianjen pentru f (x₀ = 0.2, 100 de pași)",
		xlabel = "xₙ", ylabel = "xₙ₊₁", aspect = 1)
	panza!(ax2, f, 0.2, 100, 0, 1)
	ylims!(ax2, 0, 1)
	fig
end

# ╔═╡ 6eb206ba-143d-423a-a1e4-0f1151ea024b
md"""
Observați că pînza de păianjen se „înfășoară” în jurul unui **pătrat**, cu colțurile pe grafic și pe bisectoare: acesta este ciclul ``p \to q \to p``. Pentru ``r > 3.449`` pătratul se dedublează.

### ✅ Verificare rapidă
Pentru ``r = 3.2``, calculați ``p`` și ``q`` și verificați că ``f(p) = q``.
"""

# ╔═╡ d1bdcb58-19f6-4759-8c16-49eb24fcff8b
hint(md"``\sqrt{4.2 \cdot 0.2} = \sqrt{0.84} \approx 0.9165``, deci ``p \approx (4.2 + 0.9165)/6.4 \approx 0.7995`` și ``q \approx (4.2 - 0.9165)/6.4 \approx 0.5130``. Verificare: ``3.2 \cdot 0.7995 \cdot 0.2005 \approx 0.5130``. Populația alternează între aproximativ 80% și 51% din capacitatea de suport.")

# ╔═╡ 62ff91b1-fadf-435f-8a4e-0f2a146cf252
md"""
## 🌪️ 5. Diagrama de bifurcație și haosul

### 5.1. Diagrama de bifurcație
Pentru fiecare valoare a lui ``r`` iterăm modelul 500 de generații (ca să dispară regimul tranzitoriu) și desenăm următoarele 150 de valori ale lui ``x_n``. Pe axa orizontală avem parametrul ``r``, iar pe cea verticală **comportamentul pe termen lung** al populației:
* un singur punct deasupra lui ``r`` ⟹ echilibru stabil;
* două puncte ⟹ ciclu de perioadă 2; patru puncte ⟹ perioadă 4, ...
* o „ceață” de puncte ⟹ haos.

Diagrama este „amprenta” întregului model. Graficul de dedesubt arată **exponentul Liapunov** ``\lambda(r)``, despre care vorbim în secțiunea 5.3. Folosiți slider-ele pentru a mări o porțiune a diagramei.

``r`` minim: $(@bind r_min PlutoUI.Slider(1.0:0.01:3.95, show_value=true, default=2.5))

``r`` maxim: $(@bind r_max PlutoUI.Slider(1.05:0.01:4.0, show_value=true, default=4.0))
"""

# ╔═╡ 05a9b88d-8511-4a57-9665-a56b88b64a8b
let
	r_sus = max(r_max, r_min + 0.01)
	rs = range(r_min, r_sus, length = 1200)
	n_tranz, n_desen, n_liap = 500, 150, 1000
	r_pct, x_pct = Float64[], Float64[]
	liapunov = zeros(length(rs))
	for (i, r) in enumerate(rs)
		x = 0.3
		for _ in 1:n_tranz
			x = f_logistic(x, r)
		end
		s = 0.0
		for k in 1:n_liap
			s += log(abs(r * (1 - 2x)) + 1e-15)
			x = f_logistic(x, r)
			if k <= n_desen
				push!(r_pct, r); push!(x_pct, x)
			end
		end
		liapunov[i] = s / n_liap
	end

	fig = Figure(size = (1400, 800))
	ax1 = Axis(fig[1, 1], title = "Diagrama de bifurcație a modelului logistic",
		ylabel = "Comportamentul pe termen lung (xₙ)")
	# bifurcațiile de dublare a perioadei și începutul haosului
	r_bif = filter(r -> r_min < r < r_sus, [3.0, 1 + sqrt(6), 3.5441, 3.5699])
	# Makie nu poate desena vlines! cu un vector gol (fereastră fără bifurcații marcate)
	isempty(r_bif) || vlines!(ax1, r_bif, color = (:crimson, 0.5), linestyle = :dash)
	scatter!(ax1, r_pct, x_pct, markersize = 1.5, color = (:navy, 0.25), rasterize = 2)

	ax2 = Axis(fig[2, 1], title = "Exponentul Liapunov",
		xlabel = "Parametrul r", ylabel = "λ(r)")
	hlines!(ax2, 0, color = :black, linewidth = 1)
	band!(ax2, rs, zeros(length(rs)), max.(liapunov, 0), color = (:crimson, 0.3))
	lines!(ax2, rs, max.(liapunov, -2.5), color = :navy, linewidth = 1.5)
	ylims!(ax2, -2.5, 1)

	linkxaxes!(ax1, ax2)
	xlims!(ax2, r_min, r_sus)
	hidexdecorations!(ax1, grid = false)
	rowsize!(fig.layout, 2, Relative(0.3))
	fig
end

# ╔═╡ 188dfdbb-4e1c-4dec-9741-108142504561
md"""
### 🗺️ Harta comportamentelor modelului logistic

| Interval pentru ``r`` | Comportamentul sistemului |
| :--- | :--- |
| ``0 < r \le 1`` | 📉 **Extincție**: orice populație tinde la ``0``. |
| ``1 < r \le 2`` | 🎯 **Echilibru stabil** ``1 - \frac{1}{r}``, atins monoton. |
| ``2 < r \le 3`` | 🎯 **Echilibru stabil**, atins prin oscilații amortizate. |
| ``3 < r < 3.449`` | 🔄 **Ciclu de perioadă 2**. |
| ``3.449 < r < 3.544`` | 🔄 **Ciclu de perioadă 4**. |
| ``3.544 < r < 3.5699`` | 🔄 **Perioade 8, 16, 32, …** (cascada dublărilor de perioadă). |
| ``3.5699 < r \le 4`` | 🌪️ **Haos**, întrerupt de „ferestre” periodice. |

Liniile roșii punctate din diagramă marchează valorile ``r_1 = 3``, ``r_2 \approx 3.449``, ``r_3 \approx 3.544`` și limita ``r_\infty \approx 3.5699``. Intervalele dintre dublări se micșorează de fiecare dată de aproximativ ``4.669`` ori:
```math
\lim_{k \to \infty} \frac{r_k - r_{k-1}}{r_{k+1} - r_k} = \delta \approx 4.669 .
```
Numărul ``\delta`` (**constanta lui Feigenbaum**) este **universal**: apare la fel pentru modelul Ricker și pentru orice funcție cu o singură „cocoașă”. De aceea dublarea perioadei este considerată „drumul tipic spre haos”.

### 📝 Exercițiul de la tablă: vînătoare de ferestre periodice
> ❓ Deși pentru ``r > 3.57`` sistemul este în general haotic, există mici „ferestre” de ordine. Măriți diagrama de bifurcație (de exemplu, ``r`` între ``3.8`` și ``3.9``) și găsiți valoarea lui ``r`` de la care apare un **ciclu de perioadă 3** (trei benzi orizontale clare).
>
> *Valoarea găsită: ``r \approx`` _______*
>
> Cum arată exponentul Liapunov în această fereastră?
"""

# ╔═╡ 36c2d7fd-4356-4af2-a394-ae46993a7a3d
hint(md"Fereastra de perioadă 3 începe la ``r = 1 + \sqrt{8} \approx 3.828``. În interiorul ei exponentul Liapunov este negativ: ciclul este stabil. Celebra teoremă a lui Li și Yorke (1975), intitulată *Period three implies chaos*, arată că un sistem care are un ciclu de perioadă 3 are cicluri de **orice** perioadă.")

# ╔═╡ 28c40066-0753-4f87-aed7-9e5a164fba7a
md"""
### 5.2. Sensibilitatea la condițiile inițiale
Să urmărim două populații care pornesc din valori aproape identice, ``x_0 = 0.2`` și ``x_0 = 0.2 + \varepsilon``. Diferența ``\varepsilon`` poate fi o eroare de numărare: nimeni nu poate număra exact o populație de insecte.

Parametrul ``r``: $(@bind r_fl PlutoUI.Slider(2.5:0.01:4.0, show_value=true, default=3.9))

Diferența inițială ``\varepsilon = 10^{k}``, cu ``k`` = $(@bind k_eps PlutoUI.Slider(-12:1:-1, show_value=true, default=-6))
"""

# ╔═╡ a2198010-22d4-4e79-82f9-18f5a8ed9dbe
let
	f(x) = f_logistic(x, r_fl)
	n = 60
	ε = 10.0^k_eps
	a = itereaza(f, 0.2, n)
	b = itereaza(f, 0.2 + ε, n)
	dist = max.(abs.(a .- b), 1e-17)

	# exponentul Liapunov, estimat pe o traiectorie lungă
	lung = itereaza(f, 0.2, 10_500)[501:end]
	λ = sum(x -> log(abs(r_fl * (1 - 2x)) + 1e-15), lung) / length(lung)

	fig = Figure(size = (1100, 650))
	ax1 = Axis(fig[1, 1], title = @sprintf("Două traiectorii, r = %.2f,  ε = 10^%d", r_fl, k_eps),
		ylabel = "xₙ", limits = (0, n, 0, 1))
	scatterlines!(ax1, 0:n, a, color = :teal, markersize = 6, label = "x₀ = 0.2")
	scatterlines!(ax1, 0:n, b, color = :darkorange, markersize = 6, linestyle = :dash,
		label = "x₀ = 0.2 + ε")
	axislegend(ax1, position = :lb, orientation = :horizontal)

	subtitlu = λ > 0 ?
		@sprintf("λ ≈ %.3f > 0: eroarea crește ca ε·e^(λn); predicția devine inutilă după n ≈ ln(1/ε)/λ ≈ %.0f generații", λ, log(1 / ε) / λ) :
		@sprintf("λ ≈ %.3f < 0: diferența se stinge, traiectoriile se apropie", λ)
	ax2 = Axis(fig[2, 1], title = "Distanța dintre traiectorii (scară logaritmică)", subtitle = subtitlu,
		xlabel = "Generația (n)", ylabel = "|xₙ − yₙ|", yscale = log10, limits = (0, n, 1e-17, 10))
	lines!(ax2, 0:n, dist, color = :purple, linewidth = 2)
	if λ > 0
		lines!(ax2, 0:n, min.(ε .* exp.(λ .* (0:n)), 1), color = :crimson, linestyle = :dash,
			label = "ε·e^(λn)")
		axislegend(ax2, position = :rb)
	end
	fig
end

# ╔═╡ 83cd0f6c-d829-47f4-b679-d5ac5076ac81
md"""
Pentru ``r = 3.9`` cele două traiectorii coincid la început, dar după cîteva zeci de generații devin complet diferite. Micșorînd ``\varepsilon`` de un milion de ori cîștigăm doar cîteva generații de predicție. Pentru ``r < 3.57`` (alegeți, de exemplu, ``r = 3.2``) diferența se stinge.

### 5.3. Exponentul Liapunov
De unde vine creșterea exponențială a erorii? După ``n`` pași, cu regula de derivare a funcțiilor compuse,
```math
\lvert f^n(x_0 + \varepsilon) - f^n(x_0)\rvert \approx \lvert (f^n)'(x_0)\rvert\, \varepsilon
= \lvert f'(x_0)\, f'(x_1) \cdots f'(x_{n-1})\rvert\, \varepsilon = e^{\lambda n}\, \varepsilon,
```
unde
```math
\lambda = \frac{1}{n} \sum_{k=0}^{n-1} \ln \lvert f'(x_k)\rvert
```
este media logaritmilor factorilor de amplificare de-a lungul traiectoriei: **exponentul Liapunov**.
* ``\lambda < 0``: erorile se sting ⟹ traiectoria este atrasă de un echilibru sau de un ciclu stabil;
* ``\lambda > 0``: erorile cresc exponențial ⟹ **haos**.

Pentru ``r = 4`` se poate arăta că ``\lambda = \ln 2``: eroarea se dublează la fiecare generație, deci o eroare inițială de ``10^{-6} \approx 2^{-20}`` devine de ordinul populației însăși după circa 20 de generații.

### 5.4. Ce înseamnă „haos”?
Haosul **nu** înseamnă dezordine întîmplătoare. Un sistem haotic este:
1. **determinist**: nu conține nimic aleator, ``x_{n+1}`` este complet determinat de ``x_n``;
2. **mărginit**: populația nu explodează și nu dispare;
3. **neregulat**: traiectoria nu se repetă niciodată periodic;
4. **sensibil la condițiile inițiale** (*efectul fluture*): ``\lambda > 0``.

> 🌍 **Consecințe ecologice.**
> * Fluctuațiile mari și neregulate ale unei populații **nu sînt neapărat** cauzate de factori externi întîmplători (vreme, boli); ele pot fi produse de dinamica internă a populației. Acesta a fost mesajul principal al articolului lui May din 1976.
> * Haosul a fost observat experimental: în anii '90, Costantino și colaboratorii au împins populații de gîndaci de făină (*Tribolium*) în regim haotic, modificînd în laborator rata de mortalitate a adulților.
> * Pentru o populație haotică, predicțiile pe termen lung sînt imposibile în practică, chiar dacă modelul este cunoscut perfect.
"""

# ╔═╡ 7782ae4f-db38-4ab1-896a-155c7115628a
md"""
## 🐟 6. Alte modele: Ricker și Beverton–Holt

Modelul logistic discret are un defect: dacă ``x_n > 1``, atunci ``x_{n+1} < 0``, ceea ce nu are sens biologic. În gestionarea pescuitului se folosesc două modele care evită această problemă (cu capacitatea de suport normalizată la ``K = 1``):

1. **Verhulst (logistic):** ``x_{n+1} = r\, x_n (1 - x_n)``.
2. **Ricker** (1954, populații de somon): ``x_{n+1} = x_n\, e^{r(1 - x_n)}``. Cînd populația este mare, urmașii sînt decimați (canibalism, boli), dar populația nu devine niciodată negativă.
3. **Beverton–Holt** (1957): ``x_{n+1} = \dfrac{r\, x_n}{1 + (r-1)\, x_n}``. Numărul de urmași crește cu populația, dar se saturează (**compensare**, fără supracompensare).

Toate trei au echilibrele ``x^* = 0`` și ``x^* = 1`` (pentru logistic, ``x^* = 1 - 1/r``). Criteriul derivatei (Exercițiile 1 și 2) dă:

| Model | ``f'(x^*)`` în echilibrul nenul | Echilibrul nenul este stabil pentru | Haos? |
| :--- | :---: | :---: | :---: |
| Verhulst | ``2 - r`` | ``1 < r < 3`` | da, ``r > 3.57`` |
| Ricker | ``1 - r`` | ``0 < r < 2`` | da, ``r > 2.69`` |
| Beverton–Holt | ``1/r`` | orice ``r > 1`` | nu |

Comparați cele trei modele, folosind același ``r`` și aceeași populație inițială:

Parametrul ``r``: $(@bind r_cmp PlutoUI.Slider(0.1:0.1:4.0, show_value=true, default=2.3))

Populația inițială ``x_0``: $(@bind x0_cmp PlutoUI.Slider(0.05:0.05:1.0, show_value=true, default=0.2))
"""

# ╔═╡ 68ba3408-3ef6-452e-a011-2e17b8305fc0
let
	n_pasi = 60
	modele = [
		("Verhulst (logistic)", x -> f_logistic(x, r_cmp), :red, 1.0),
		("Ricker", x -> f_ricker(x, r_cmp), :green, max(3.0, 1.1 * exp(r_cmp - 1) / r_cmp)),
		("Beverton–Holt", x -> f_beverton(x, r_cmp), :blue, r_cmp >= 1 ? 2.0 : 1.0),
	]

	fig = Figure(size = (1400, 800))
	for (j, (nume, f, culoare, xmax)) in enumerate(modele)
		# rîndul 1: pînza de păianjen
		ax = Axis(fig[1, j], title = nume, xlabel = "xₙ", ylabel = "xₙ₊₁")
		panza!(ax, f, x0_cmp, 30, 0, xmax, culoare = culoare)
		ylims!(ax, 0, xmax)

		# rîndul 2: evoluția în timp
		ax2 = Axis(fig[2, j], title = "Simulare $nume  (r = $(r_cmp))",
			xlabel = "Generația (n)", ylabel = "Populație")
		sim = itereaza(f, x0_cmp, n_pasi)
		hlines!(ax2, [0], color = :gray)
		scatterlines!(ax2, 0:n_pasi, sim, color = culoare, markersize = 6)
	end
	rowgap!(fig.layout, 30)
	colgap!(fig.layout, 25)
	fig
end

# ╔═╡ c817bb66-b11f-4b0b-9dce-b48e00a04bc8
md"""
**Ce observăm:**
* Pentru ``r = 2.3`` modelul logistic oscilează amortizat, modelul Ricker are deja un ciclu de perioadă 2 (de ce? comparați cu tabelul), iar Beverton–Holt urcă liniștit către ``1``.
* Curba Beverton–Holt este **crescătoare**: mai mulți părinți înseamnă întotdeauna mai mulți urmași, deci nu poate apărea supracompensarea. Curbele Verhulst și Ricker au o „cocoașă”: dincolo de un anumit nivel, mai mulți părinți înseamnă **mai puțini** urmași.
* Modelul logistic are sens doar pentru ``0 \le x \le 1`` (altfel produce populații negative), pe cînd modelele Ricker și Beverton–Holt rămîn pozitive pentru orice ``x_0 > 0``.

> 💡 Forma funcției de reproducere ``f`` (cu sau fără cocoașă) contează mai mult decît formula exactă: **supracompensarea** este cea care produce oscilații, cicluri și haos.
"""

# ╔═╡ 6f8fcc5a-f312-4879-9aeb-2c50556dc200
md"""
## 📝 7. Recapitulare și exerciții

### 🧠 Recapitulare
- ⏱️ Speciile cu **generații care nu se suprapun** se modelează prin ecuații cu diferențe ``x_{n+1} = f(x_n)``.
- 📈 **Creștere geometrică:** ``x_{n+1} = \lambda x_n \implies x_n = \lambda^n x_0``; legătura cu timpul continuu este ``\lambda = e^r``.
- ⚖️ **Echilibre:** ``x^* = f(x^*)``, intersecția graficului lui ``f`` cu bisectoarea. **Stabilitate:** ``\lvert f'(x^*)\rvert < 1``; dacă ``f'(x^*) < 0``, apropierea se face cu oscilații.
- 🕸️ **Diagrama pînză de păianjen:** vertical pînă la grafic, orizontal pînă la bisectoare.
- 🔄 **Ciclu de perioadă 2:** echilibru al lui ``f \circ f``, stabil dacă ``\lvert f'(p)\, f'(q)\rvert < 1``.
- 🌪️ **Modelul logistic:** dublări de perioadă la ``r = 3,\ 3.449,\ 3.544, \dots``, haos după ``r_\infty \approx 3.57``; haosul este **determinist**, dar **imprevizibil** (exponent Liapunov ``\lambda > 0``).
- 🐟 **Ricker** are haos, **Beverton–Holt** nu: oscilațiile vin din **supracompensare**.

---
### ✍️ Exerciții
**Exercițiul 1 (modelul Ricker).** Fie ``x_{n+1} = x_n e^{r(1 - x_n)}``, cu ``r > 0``.

a) Găsiți punctele de echilibru (rezolvați ``x = f(x)``).

b) Calculați ``f'(x)`` și arătați că echilibrul nenul este stabil pentru ``0 < r < 2``. Cînd se face apropierea monoton și cînd cu oscilații?

c) Folosiți funcțiile `itereaza` și `panza!` pentru a desena pînza de păianjen pentru ``r = 1.5`` și ``r = 2.3``.
"""

# ╔═╡ 16d50a4b-9b4d-4454-8d75-2e6db9c61e45
hint(md"a) ``x = x e^{r(1-x)}`` dă ``x^* = 0`` sau ``e^{r(1-x)} = 1``, adică ``x^* = 1``. b) ``f'(x) = e^{r(1-x)}(1 - rx)``, deci ``f'(0) = e^r > 1`` și ``f'(1) = 1 - r``. Condiția ``\lvert 1 - r\rvert < 1`` înseamnă ``0 < r < 2``; monoton pentru ``r < 1``, cu oscilații pentru ``1 < r < 2``.")

# ╔═╡ ae98bc88-1aa3-4fb7-9755-20a146a8723f
# EXERCIȚIUL 1 c): completați și rulați
begin
fig_ex1 = Figure(size = (900, 450))
ax_ex1 = Axis(fig_ex1[1, 1], xlabel = "xₙ", ylabel = "xₙ₊₁", title = "Ricker, r = 1.5")
panza!(ax_ex1, x -> f_ricker(x, 1.5), 0.2, 30, 0, 3)
ax_ex2 = Axis(fig_ex1[1, 2], xlabel = "xₙ", ylabel = "xₙ₊₁", title = "Ricker, r = 1.5")
panza!(ax_ex2, x -> f_ricker(x, 2.3), 0.2, 30, 0, 3)
fig_ex1
end

# ╔═╡ 5317987e-98c6-4a3c-877a-8e9cc60f2e51
md"""
**Exercițiul 2 (modelul Beverton–Holt).** Fie ``x_{n+1} = \dfrac{r x_n}{1 + (r-1) x_n}``, cu ``r > 1``.

a) Arătați că ``f'(1) = 1/r`` și deduceți că echilibrul ``x^* = 1`` este întotdeauna stabil, iar apropierea de el este monotonă.

b) Cu substituția ``y_n = 1/x_n``, arătați că ``y_{n+1} = \frac{1}{r}\, y_n + \frac{r-1}{r}`` și rezolvați această recurență liniară.

c) Deduceți formula ``x_n = \dfrac{1}{1 + \left(\dfrac{1}{x_0} - 1\right) r^{-n}}`` și comparați-o cu soluția ecuației logistice continue din Cursul 3. Ce rol joacă ``r^{-n}``?
"""

# ╔═╡ 2804940d-cf50-40e0-964c-6b2a6ad14690
hint(md"a) ``f'(x) = r / (1 + (r-1)x)^2``, deci ``f'(1) = r/r^2 = 1/r \in (0, 1)``. b) Abaterea ``y_n - 1`` se înmulțește cu ``1/r`` la fiecare pas, deci ``y_n = 1 + (y_0 - 1)\, r^{-n}``. c) Soluția continuă era ``K / (1 + (K/x_0 - 1)e^{-rt})``; cu ``K = 1``, formula discretă este aceeași, cu ``e^{-rt}`` înlocuit de ``r^{-n}``.")

# ╔═╡ 7425776c-1b87-4446-94e8-4f55c1c7c912
md"""
**Exercițiul 3 (metoda Euler poate produce haos).** Aplicăm metoda Euler cu pasul ``\Delta t`` ecuației logistice continue ``x' = \rho\, x(1 - x/K)``.

a) Arătați că, notînd ``R = 1 + \rho \Delta t`` și ``y_n = \dfrac{\rho\Delta t}{K R}\, x_n``, schema Euler devine exact modelul logistic discret ``y_{n+1} = R\, y_n(1 - y_n)``.

b) Pentru ``\rho = 0.5``, ce pas ``\Delta t`` face ca soluția numerică să oscileze? Dar să devină haotică? Ce concluzie trageți despre alegerea pasului în simulări?
"""

# ╔═╡ 9fad8262-63bd-4681-9336-f098bbb5d30f
hint(md"b) Oscilații apar pentru ``R > 2``, adică ``\Delta t > 2``; ciclul de perioadă 2 pentru ``R > 3``, adică ``\Delta t > 4``; haos pentru ``R > 3.57``, adică ``\Delta t > 5.14``. Soluția exactă a ecuației continue este însă întotdeauna monotonă: oscilațiile și haosul sînt artefacte numerice ale unui pas prea mare.")

# ╔═╡ 1da460a7-91a8-4fb7-b704-99ea528e991f
md"""
**Exercițiul 4 (discuție).** 🎣 Sînteți responsabili cu gestionarea populației de pești dintr-un lac. Datele arată că populația urmează un model logistic discret cu ``r = 3.8``.

a) Pe baza diagramei de bifurcație și a exponentului Liapunov, ce i-ați spune directorului despre predictibilitatea populației pe termen lung? De ce este periculos să presupuneți că populația va fi aceeași de la un an la altul?

b) Pescuitul reduce rata efectivă de reproducere ``r``. Ce s-ar putea întîmpla cu dinamica populației dacă se pescuiește moderat? Dar dacă se pescuiește excesiv?
"""

# ╔═╡ 3f2caf5d-f427-4320-aed4-bca3377d25a5
hint(md"a) Pentru ``r = 3.8`` sistemul este haotic (``\lambda > 0``): orice eroare de numărare crește exponențial, deci previziunile au sens doar pe cîțiva ani; cote fixe calculate pe baza unei populații „medii” pot fi catastrofale într-un an slab. b) Scăzînd ``r`` coborîm în diagrama de bifurcație spre stînga: un pescuit moderat poate **stabiliza** populația (cicluri, apoi echilibru), dar dacă ``r`` scade sub ``1``, populația dispare.")

# ╔═╡ Cell order:
# ╟─e7e9e376-adb0-11f1-91f0-f7cc7081a2ec
# ╟─ce44eb4e-4b93-43fa-9105-1a10e9015ccd
# ╟─f612f2c6-5756-406d-ae1d-f4d55b03e87f
# ╟─5a49c8fe-6d17-4ce3-8865-3c0ba68af9c8
# ╟─e07e94e7-938f-4546-b9cc-71f297ba1234
# ╟─9f47681d-543a-40f2-a663-cc685f395057
# ╟─1e558223-6779-4a42-b75c-39e7c1cfb046
# ╟─aa081a08-4f61-48a9-865b-f3e487b5d231
# ╟─09fd950e-89f8-4627-adfd-ae595cacc545
# ╟─182c62ba-84ab-4233-b7d1-12a534582095
# ╟─2c8f0921-8e71-49a9-b910-9940abcd25bd
# ╟─dabc61af-b537-48ac-9e20-caca44da4b65
# ╟─2a3b4bbb-7891-455c-99e0-d820d1b4ffe3
# ╟─b46b9ded-0cd2-4fec-a51c-17baaaf52db9
# ╟─4ce78424-8bc8-4cc4-986b-3ad350e627e1
# ╟─01fa65ac-bed0-41de-b98c-1c1ca5ea5622
# ╟─7128b3da-1907-4946-a117-41cbfa3381e5
# ╟─37806866-5fb1-4f6e-8d76-0133ee350ebb
# ╟─ef07c739-7314-48ac-bad0-0206cf48f2fe
# ╟─c1b2337e-73a6-4842-8e62-7bfc97e8b519
# ╟─bc222929-5b58-48b3-8b1a-00aa502ed3e5
# ╟─3b9199fc-11d0-4908-beaa-9da23c27ffe8
# ╟─a815fcaa-7fc1-4d0b-997e-6e6fe63158dd
# ╟─4f5b0ffe-0db5-485e-94fc-6949b6416e62
# ╟─6eb206ba-143d-423a-a1e4-0f1151ea024b
# ╟─d1bdcb58-19f6-4759-8c16-49eb24fcff8b
# ╟─62ff91b1-fadf-435f-8a4e-0f2a146cf252
# ╟─05a9b88d-8511-4a57-9665-a56b88b64a8b
# ╟─188dfdbb-4e1c-4dec-9741-108142504561
# ╟─36c2d7fd-4356-4af2-a394-ae46993a7a3d
# ╟─28c40066-0753-4f87-aed7-9e5a164fba7a
# ╟─a2198010-22d4-4e79-82f9-18f5a8ed9dbe
# ╟─83cd0f6c-d829-47f4-b679-d5ac5076ac81
# ╟─7782ae4f-db38-4ab1-896a-155c7115628a
# ╟─68ba3408-3ef6-452e-a011-2e17b8305fc0
# ╟─c817bb66-b11f-4b0b-9dce-b48e00a04bc8
# ╟─6f8fcc5a-f312-4879-9aeb-2c50556dc200
# ╟─16d50a4b-9b4d-4454-8d75-2e6db9c61e45
# ╠═ae98bc88-1aa3-4fb7-9755-20a146a8723f
# ╟─5317987e-98c6-4a3c-877a-8e9cc60f2e51
# ╟─2804940d-cf50-40e0-964c-6b2a6ad14690
# ╟─7425776c-1b87-4446-94e8-4f55c1c7c912
# ╟─9fad8262-63bd-4681-9336-f098bbb5d30f
# ╟─1da460a7-91a8-4fb7-b704-99ea528e991f
# ╟─3f2caf5d-f427-4320-aed4-bca3377d25a5
