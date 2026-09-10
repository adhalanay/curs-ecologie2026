### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ a1b2587e-f150-11f0-0b8e-3baef1756e57
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	
	using PlutoTeachingTools
	
	# using LoopVectorization
	TableOfContents()
end

# ╔═╡ cb55a251-0a77-4326-b9e9-d1b9e5705ff8
begin
	using Distributions,Plots
	n = 10.0
	p = 0.6
	d = Binomial(n, p)
	x_values = support(d)
	pmf_values = pdf.(d, x_values)
	k=0:11
	bar(k, pmf_values,title="Distribuția binomială cu n=$n și p=$p", legend=false)
end

# ╔═╡ 213ad65f-ef04-4e20-8702-474493326ecc
WidthOverDocs()

# ╔═╡ 59ba7db8-6a87-45da-9454-a46a2747c4c9
md"""
# Variabile aleatoare

## Variabile aleatoare discrete

Adesea în experimentele probabilistice sîntem mai interesați de anumite proprietăți numerice ale rezultatelor decît în 
rezultatele înseși. De exemplu aruncăm mai multe monede și sîntem interesați cîte „capete” apar, fără să ne inereseze cum 
obținem aceasta. 

O **variabilă aleatoare discretă** este o funcție $X:\Omega \to M$, unde $M$ este o submulțime discretă a lui $\mathbb{R}.$
Dacă aruncăm două monede atunci $\Omega=\left\{HH,HT,TH,TT\right\}.$ Fie $X$ variabila care numără de cîte ori apare $T$.
$X:\Omega \to \left\{0,1,2\right\}$, $X(HH)=0$, $X(HT)=X(TH)=1$ și $X(TT)=2$.

Fiecărei variabile aleatoare îi asociem mai multe funcții. Prima este **densitatea de probabilitate**: $p_i=P(X=i).$ Densitatea
de probabilitate se reprezintă de obicei prin histograme ale căror înălțime este aleasă astfel încît aria este proporțională cu
$p_i$. Avem întotdeauna că $p_i \geq 0$ pentru orice $i$ și mai ales
```math
\sum_i p_i=1.
```
Deci aria totală a histogramei este $1$.

O altă mărime importantă asociată unei variabile aleatoare $X$ este **densitatea cumulativă de probabilitate** a lui $X:$ 
```math
F_i=P(X\leq i).
```

Pentru variabile aleatoare discrete avem ca și în cazul statisticii media și varianța. **Media** unei variabile discrete $X$
este
```math
E[X]=\sum_i ip_i,
```
unde $p_i$ este densitatea de probabilitate. Tradițional, media se mai numește și valoarea așteptată a lui $X$ sau centrul de 
masă.

**Varianța** lui $X$ măsoară cît de mult se poate îndepărta $X$ de medie:
```math
Var[X]=\sum_i\left(i-E[X]\right)^2p_i.
```

**Deviația standard** a lui $X$ este radicalul lui $Var[X].$

Două exemple de variabile aleatoare des întîlnite sînt variabilele de tip Bernoulli și cele de tip binomial. 

Variabilele Bernoulli modelează evenimentele care pot avea doar două valori. O **variabilă Bernoulli** cu pondere $p$ este 
o variabilă cu $P(X=1)=p$ și $P(X=0)=1-p$.

Variabilele binomiale generalizează variabilele Bernoulli la mai multe experimente independente. O **variabilă binomială** cu
parametrii $p$ și $n$ ia valori în $\left\{0,1,\dots,n\right\}$ și are densitatea de probabilitate
```math
P(X=i)=C_n^ip^i(1-p)^{n-i}.
```

La oameni probabilitatea ca un nou născut să fie fată este aproximativ 48%. Considerăm o familie cu 4 copii. Care este 
probabilitatea ca 3 din copii să fie fete, ținînd cont de faptul că sexele copiilor sînt independente?

Sexul fiecărui copil este modelat de o distrbuție Bernoulli cu $n=4$ și $p=0.48$.
```math
P(X=i)=C_4^i 0.48^i(1-0.48)^{4-i}=\frac{4!}{i!(4-i)!}0.48^i 0.52^{4-i}.
```
Răspunsul este dat de calculul pentru $i=3:$
```math
P(X=3)=\frac{4!}{3!1!}0.48^30.52=4\cdot0.48^30.52 \approx 0.23.
```

"""

# ╔═╡ e7e0cd69-46de-46bb-a62a-1e0652b6fc8b
md"""
Alte tipuri de distribuții pentru variabile discrete sînt:
  - distribuția **multinomială**: corespunde unui experiment în care avem mai mule rezultate posibile;
  - distribuția **geometrică**: considerăm un experiment cu două rezultate posibile și ne oprim atunci cînd obținem un succes: ``P(X=k)=(1-p)^kp``. ``X=k`` semnifică faptul că primele ``k-1`` experimente au eșuat și al ``k-``ulea a avut succes;
  - distrbuția **Poisson** este una dintre cele mai importante distribuții in biologie modelează de exemplu substituțiile amino-acizilor în proteine, probabilitatea de evitare a parazitismului de către gazdă, Avem 
```math 
P(X=k) = e^-\lambda\frac{\lambda^k}{k!},
```
   unde ``\lambda > 0`` este un parametru real.
"""

# ╔═╡ 85836ef2-e7cc-42ac-8b2f-59535ba95f13
begin
	# n = 10.0
	p1 = 0.3
	d1 = Geometric(p1)
	max_trials = quantile(d1, 0.999) + 1
	trials = 1:max_trials
	pmf_trials = [pdf(d1, k-1) for k in trials]
	bar(trials, pmf_trials,title="Distribuția geometrică cu p=$p1", legend=false)
end

# ╔═╡ c6ae7f67-126f-448d-b5a6-6e92f52b58e5
@bind λ PlutoUI.Slider(0.1:0.1:15,show_value=true)

# ╔═╡ 7fff8acf-e940-4d43-8d2a-dd1c4cd97dbc
begin
	dist = Poisson(λ)
    k2 = 0:min(50, ceil(Int, λ + 4√λ))  # Reasonable range
    pmf = pdf.(dist, k2)
    
    bar(k2, pmf,
        xlabel = "k",
        ylabel = "P(X = k)",
        title = "Distribuția Poisson cu λ = $λ",
        legend = false,
        color = :viridis,
        fillalpha = 0.7,
        ylims = (0, 0.5),
        size = (800, 400))
    
    # Add mean and mode
    vline!([λ], color=:red, linestyle=:dash, label="Media = $λ")
    
    # Mode is floor(λ) for Poisson
    mode_val = floor(Int, λ)
    Plots.scatter!([mode_val], [pdf(dist, mode_val)], 
             color=:black, markersize=8, label="Modulul = $mode_val")
end

# ╔═╡ 8a4985ab-583d-4ace-9c1b-bbe58bed81cd
md"""
## Variabile aleatoare continue

O variabilă aleatoare continuă este o funcție $X:\Omega \to I$, unde $I$ este un interval din $\mathbb{R}$. Putem să ne gîndim
la o variabilă aleatoare continuă ca provenind dintr-una discretă, dar cu valorile cutiilor din ce în ce mai mici.

Ca și în cazul variabilelor discrete putem asocia unei variabile aleatoare continue o densitate de probabilitate, adică o 
funcție $f(x)$ pentru care 
```math
P(a \leq X \leq b)=\int_a^b f(x) dx.
```

Analog cu cazul discret $f(x) \geq 0$ și
```math
\int_{-\infty}^{\infty} f(x) dx = 1.
```
Adică aria de sub graficul lui $f$ este egală cu $1$. De asemenea se poate defini și densitatea cumulativă de probabilitate
prin
```math
F(x)=P(X\leq x).
```

Avem noțiunile analoage de medie și varianță pentru variabile aleatoare continue. 
```math
E[X]=\int_{-\infty}^{\infty}xf(x) dx.
```
Ca și înainte media mai este numită și valoare așteptată.

Varianța are formula
```math
Var[X]=\int_{-\infty}^{\infty}(x-E[X])^2f(x)dx
```
sau dacă efectuăm calculele
```math
Var[X]=\int_{-\infty}^{\infty}x^2f(x)dx-E[X]^2.
```
Multe fenomene au probabilitatea descendentă exponențial. O variabilă aleatoare $T$ cu parametru $c$ ia valori în $[0,\infty)$
și are densitatea de probabilitate $f(x)=ce^{-cx}.$ Media lui $T$ este $\frac{1}{c}$.

Alte fenomene sînt modelate de variabile aleatoare **normale**. O variabilă aleatoare normală are distribuția
```math
f(x)=\frac{1}{\sigma\sqrt{2\pi}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}.
```
Media unei astfel de variabile este $\mu$, iar varianța este $\sigma^2$.

Adesea este util să considerăm variabile normale cu media $0$ și varianța $1$.
"""

# ╔═╡ 81864038-3ac2-4cba-b443-4146ec751462
begin
	n1 = 40.0
	p2 = 0.6
	m = n1*p2
	s = sqrt(n1*p2*(1-p2))
	d2 = Binomial(n1, p2)
	x_values1 = support(d2)
	pmf_values1 = pdf.(d2, x_values1)
	k1=0:40
	bar(k1, pmf_values1,title="Distribuția binomială cu n=$n1 și p=$p2 și distribuția \n normală corespunzătoare", legend=false)
	f(x) = 1/(s*sqrt(2*pi))*exp(-1/2*((x-m)/s)^2)
	xs = LinRange(0,40,400)
	ys = [f(x) for x ∈ xs]
	Plots.plot!(xs,ys)
end

# ╔═╡ 1d099793-df63-4e40-9b5c-fa9dfe347c86
md"""
## Alte distrbuții continue

- distribuția **uniformă** este de forma
```math
f(x)=\left\{\begin{array}{ccc} \frac{1}{b-a} \text{ dacă } x \in (a,b) \\
0 \text{ în rest. }
\end{array}\right.
```
Media este ``\frac{a+b}{2}``. iar varianța este ``\frac{(b-a)^2}{12}.``
 - distribuția **exponențială** are forma 
```math
f(x)=\left\{\begin{array}{ccc} \lambda e^{-\lambda x } \text{ dacă } x \geq 0 \\
0 \text{ dacă } x < 0.
\end{array}\right.
```
 - distrbuția Poisson : 
```math
f_{n,\lambda}(x)=\frac{\lambda^nx^{n-1}}{(n-1)!}e^{-\lambda x} \text{ pentru } x \geq 0.
```  
 - distribuția Weibull: 
```math
f_{\lambda,k}(x)=\left\{\begin{array}{ccc} \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}\exp\left(\left(-\frac{x}{\lambda}\right)^k\right) \text{ dacă } x \geq 0 \\
0 \text{ dacă } x < 0.
\end{array}\right.
```
"""

# ╔═╡ 7a82278d-cfa0-4f00-b928-b30c88d3afd4
begin
	xs0 = LinRange(0,20,100)
	a=5
	b=10
	λ0 = 0.3
	λp = 1.9
	n0 = 20
	function fu(x) 
		if x > a && x < b 
			1/(b-a)
		else
			0
		end
	end
	function fe(x)
		if x ≥ 0 
			λ0*exp(-λ0*x) 
		else 
			0
		end
	end
	function fp(x)
		if x ≥ 0 
			λp^n0*x^(n0-1)/factorial(n0-1)*exp(-λp*x) 
		else 
			0
		end
	end
	fw(x) = x/λp* (x/λp)^(n0-1)*exp(-(x/λp)^n0)
	ys_u = [fu(x) for x ∈ xs0]
	ys_e = [fe(x) for x ∈ xs0]
	ys_p = [fp(x) for x ∈ xs0]
	ys_w = [fw(x) for x ∈xs0]
	Plots.plot(xs0,ys_u,title="Distribuția uniformă cu a=$a și b=$b",legend=false)
end

# ╔═╡ 8a820ac3-853c-4491-9cba-f8087653d6a8
Plots.plot(xs0,ys_e,title="Distribuția exponențială cu λ=$λ0",legend=false)

# ╔═╡ 9f87ab7d-d2c4-4f0d-8ef9-284cdcb8a2e8
Plots.plot(xs0,ys_p, title="Distribuția Poisson cu n=$n0 și λ=$λp",legend=false)

# ╔═╡ da3e8bd0-b23e-4f1d-8de7-b775ff31d6b3
md"""
k = $(@bind kw PlutoUI.Slider(0.1:0.1:5.0,show_value=true))

λ = $(@bind λw PlutoUI.Slider(0.1:0.1:2.0,show_value=true))
"""

# ╔═╡ fe3343e3-f3c8-4cdb-b91b-f01ffb311262
let
function fw(x)
		if x ≥ 0 
			kw/λw*(x/λw)^(kw-1)*exp(-(x/λw)^kw)
		else 
			0
		end
	end
	ys_w = [fw(x) for x ∈ xs0]
	Plots.plot(xs0,ys_w, title="Distribuția Weibull pentru k=$kw și λ=$λw",legend=false)
end

# ╔═╡ ee127f8e-325f-4004-93f9-7f6c3064d2c7
md"""
# Legea numerelor mari

Nivelul antigenului specifice PSA este folosit pentru a detecta cancerul de prostată. Nivelul peste 0.5 ng/ml arată că au rămas celule canceroase după rezecția prostatei. Într-un studiu au fost urmăriți 429 de pacienți. După 5 ani 8% au avut nivelul crescut. Un alt studiu pe 30 de pacienți a constat creștere la 3 pacienți, adică
în 10% din cazuri. Care dintre rezultate este mai credibil? Evident primul. Aceasta este o manifestare a **legii numerelor mari** care spune că o estimare este mai bună cu cît eșantionul este mai mare.

Riguros, să presupunem că avem un șir de variabile aleatoare independente $X_1,\dots,X_n$ toate cu aceeași distribuție. Presupunem că $E(X_i)=\mu$ și $Var(X_i)=\sigma^2$
pentru orice $i$. Ne interesează variabila medie aritmetică:
```math
\overline{X}_n=\frac{1}{n}\sum_{i=1}^n X_i.
```

Legea (slabă) a numerelor mari spune că șirul mediilor $\left\{\overline{X}_n\right\}_{n \in \mathbb{N}}$ tinde în probabilitate la $\mu.$ Adică 
```math
\lim_{n \to \infty} Prob(|\overline{X}_n -\mu| > \epsilon)=0
```
pentru orice $\epsilon > 0$.

# Teorema limită centrală
cesta este unul din cele mai importante rezultate din teoria probabilităților. Acesta ne spune că pentru un șir de variabile i.i.d. avînd media și varianța finite, sumele formează un șir de variabile aleatoare aproximativ normal distribuite. Concret, dacă notăm 
```math
S_n =\sum_{i=1}^n X_i,
```
atunci 
```math
Prob\left(\frac{S_n - n\mu}{\sqrt{n\sigma^2}} \leq x\right) \to F(x),
```
unde 
```math
F(x)=\frac{1}{\sqrt{2\pi}}e^\frac{-x^2}{2}.
```
"""

# ╔═╡ Cell order:
# ╠═a1b2587e-f150-11f0-0b8e-3baef1756e57
# ╠═213ad65f-ef04-4e20-8702-474493326ecc
# ╟─59ba7db8-6a87-45da-9454-a46a2747c4c9
# ╠═cb55a251-0a77-4326-b9e9-d1b9e5705ff8
# ╠═e7e0cd69-46de-46bb-a62a-1e0652b6fc8b
# ╠═85836ef2-e7cc-42ac-8b2f-59535ba95f13
# ╠═c6ae7f67-126f-448d-b5a6-6e92f52b58e5
# ╠═7fff8acf-e940-4d43-8d2a-dd1c4cd97dbc
# ╟─8a4985ab-583d-4ace-9c1b-bbe58bed81cd
# ╠═81864038-3ac2-4cba-b443-4146ec751462
# ╟─1d099793-df63-4e40-9b5c-fa9dfe347c86
# ╠═7a82278d-cfa0-4f00-b928-b30c88d3afd4
# ╠═8a820ac3-853c-4491-9cba-f8087653d6a8
# ╠═9f87ab7d-d2c4-4f0d-8ef9-284cdcb8a2e8
# ╠═da3e8bd0-b23e-4f1d-8de7-b775ff31d6b3
# ╠═fe3343e3-f3c8-4cdb-b91b-f01ffb311262
# ╠═ee127f8e-325f-4004-93f9-7f6c3064d2c7
