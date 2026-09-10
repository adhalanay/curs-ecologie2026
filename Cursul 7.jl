### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 1db058d6-c4b0-11f0-31d2-8fd829c839fa
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using DifferentialEquations
	using PlutoTeachingTools
	using DynamicalSystems
	# using LoopVectorization
	TableOfContents()
end

# ╔═╡ 07efda97-46ed-4d2c-b347-d33e55561166
begin
	let
using LinearAlgebra

A=[1 2; 4 3]
p=[1,0.5]
println("Valorile proprii sînt: ",eigvals(A))
println("Vectorii proprii sînt: ",eigvecs(A))
println("Inversa lui $(A) este:", inv(A))
	end
end

# ╔═╡ deef08bd-5eec-4e82-9940-6b53b2be99bd
WidthOverDocs()

# ╔═╡ 8c4c4fd7-0566-42be-ab06-5f448055b34e
md"""
# Algebră Liniară

Pînă acum nu am făcut decît să observăm comportamentul sistemelor dinamice, dar avem nevoie de metode riguroase pentru a studia teoria calitativă a acestora. Folosirea simulărilor este doar un prim pas pentru a intui comportamentul, dar trebuie să putem face calcule riguroase. Primul instrument teoretic de care avem nevoie este algebra liniară.

## Aplicații liniare și sisteme dinamice
Vom studia aplicații liniare ``f:\mathbb{R}^n \to \mathbb{R}^n``. Vom studia în paralel două aplicații: 
 - la sisteme discrete ``\overrightarrow{X_{N+1}}=f\left(\overrightarrow{X_N}\right)``;
 - la ecuații diferențiale ``\overrightarrow{X'}=f\left(\overrightarrow{X}\right)``.

Vom nota elementele lui ``\mathbb{R}^n`` vor fi notatate cu ``\overrightarrow{x}`` și
```math
\overrightarrow{x}=\left(\begin{array}{c}x_1 \\ \vdots \\ x_n\end{array}\right).
```

Vectorii vor fi scriși atît pe verticală cît și pe orizontală. 

Am văzut la începutul cursului că spațiul stărilor pentu un sistem dinamic este ``\mathbb{R}^n`` și la fel este și spațiul tangent în care trăiesc variațiile. Tot atunci am văzut că pentru vectori avem două operații:
  - adunarea vectorilor:
```math
\overrightarrow{x} + \overrightarrow{y} = \left(\begin{array}{ccc} x_1+y_1 \\ \vdots \\ x_n+y_n\end{array}\right)
```
  - înmulțirea cu scalari:
    ```math
    a\cdot\overrightarrow{x}=\left(\begin{array}{c} ax_1 \\ \vdots \\ ax_n \end{array}\right).
    ```
    ## Combinații liniare și baze
În ``mathbb{R}^n`` avem o mulțime specială de vectori, anume 
```math
\overrightarrow{e_1}=(1,0,\dots,0),\overrightarrow{e_2}=(0,1,\dots,0),\dots,\overrightarrow{e_n}=(0,0,\dots,1).   
```
Orice alt vector ``\overrightarrow{x}=(x_1,\dots,x_n)`` se scrie în mod unic ca ``\overrightarrow{x}=x_1\overrightarrow{e_1}+\dots +x_n\overrightarrow{e_n}.``

Spunem că mulțimea ``\left\{\overrightarrow{e_1},\dots,\overrightarrow{e_n}\right\}`` este o bază a lui ``\mathbb{R}^n``, numită *baza canonică*. Vom 
vedea mai tîrziu și alte exemple de baze.

Un vector care se obține dintr-o mulțime de alți vectori folosind doar înmulțirea cu scalari și adunarea se numește **o combinație liniară** a 
respectivilor vectori.

O bază este un sistem de vectori ``\left\{b_1,\dots,b_n\right\}`` astfel încît orice vector $x \in \mathbb{R}^n$ se scrie **în mod unic**
ca o combinație liniară ``x=x_1b_1+\dots x_nb_n``. ``x_1,\dots,x_n`` se numesc *coordonatele* lui ``x`` în raport cu baza respectivă. 

Concret dacă ``b_1=\left(b_{11},b_{21},\dots,b_{n1}\right),\dots,b_n=\left(b_{1n},\dots,b_{nn}\right)`` definim matricea ``\Lambda`` cu coloanele ``b_1,\dots,b_n``. Atunci coordonatele lui ``x`` în raport baza nouă
vor fi date de formula 
```math
\Lambda^{-1}\left(\begin{array}{c} x_1 \\ \vdots \\ x_n\end{array}\right)
```

## Aplicații liniare

Am văzut tot la început că o aplicație $f:\mathbb{R} \to \mathbb{R}$ este liniară dacă ``f(x+y)=f(x)+f(y)`` și ``f(ax)=af(x).``
În cazul ``n-``dimensional definiția este similară, doar că produsul este înlocuit de înmulțirea cu scalari:

**Definiție.** O aplicație $f:\mathbb{R}^n \to \mathbb{R}^m$ se numește liniară dacă:
```math
f(\overrightarrow{x}+\overrightarrow{y})=f(\overrightarrow{x})+f(\overrightarrow{y}); \\
f(a\cdot\overrightarrow{x})=a\cdot f(\overrightarrow{x}).
```

Din definiție nu rezultă cum arată o aplicație liniară. Dar folosind definiția putem deduce acest lucru.

Dacă ``f:\mathbb{R}^1 \to \mathbb{R}^1``, putem să ne gîndim că ``\mathbb{R}^1`` are baza formată dintr-un singur vector
``\overrightarrow{e}=1``. Atunci orice vector ``\overrightarrow{x}=x\cdot\overrightarrow{e}``, deci din liniaritate rezultă că
``f(\overrightarrow{x})=x\cdot f(1)``. Deci de îndată ce știm ``f(1)``, știm toate valorile pe care le ia $f$. Deci o aplicație liniară
este în acest caz de forma ``f(x)=x\cdot a``, unde ``a=f(1)``.

Trecem acum la cazul ``f:\mathbb{R}^2 \to \mathbb{R}``. Un vector este de forma
```math
\overrightarrow{x}=x_1\overrightarrow{e_1}+x_2\overrightarrow{e_2}.
```
Deci 
```math
f(\overrightarrow{x})=x_1 f(\overrightarrow{e_1})+x_2
f(\overrightarrow{e_2}).
```

Dacă ``f(\overrightarrow{e_1})=a`` și ``f(\overrightarrow{e_2})=b``$, atunci ``f(\overrightarrow{x})=x_1\cdot a+x_2\cdot b.``

Pentru ``f:\mathbb{R}^n \to \mathbb{R}`` situația este similară cu cazul precedent: un vector va fi
```math
\overrightarrow{x}=x_1\overrightarrow{e_1}+\dots+x_n\overrightarrow{e_n}.
```

Presupunem că
```math
f(\overrightarrow{e_1})=a_1,\dots,f(\overrightarrow{e_n})=a_n.
```
Atunci
```math
f(\overrightarrow{x})=a_1x_1+\dots+a_nx_n.
```

Deci orice funcție liniară ``\mathbb{R}^n \to \mathbb{R}`` poate fi scrisă ca o combinație liniară de componentele argumentului său.

Studiem acum aplicațiile cu valori într-un spațiu de dimensiune mai mare decît ``1``. Începem cu aplicațiile ``f:\mathbb{R}^2 \to
\mathbb{R}^2``. Ca și în cazurile precedente aplicația este determinată de valorile lui ``f`` pe elementele bazei. Presupunem
``f(\overrightarrow{e_1})=\overrightarrow{a_1}`` și ``f(\overrightarrow{e_2})=\overrightarrow{a_2}.`` Dacă 
``\overrightarrow{x}=x_1\overrightarrow{e_1}+x_2\overrightarrow{e_2}$ atunci $f(\overrightarrow{x})=x_1\overrightarrow{a_1}+x_2\overrightarrow{a_2}.`` Mai departe fiecare din
vectorii 
```math
\overrightarrow{a_1}=\left(\begin{array}{c} a_{11} \\ a_{21} \end{array}\right) \text{ și } 
\overrightarrow{a_2}=\left(\begin{array}{c} a_{12} \\ a_{22} \end{array}\right).
```
Prin urmare 
```math
f(\overrightarrow{x})=\left(\begin{array}{c} a_{11}x_1+a_{12}x_2 \\ a_{21}x_1 + a_{22}x_2 \end{array}\right).
```
Deci ``f`` este unic determinată de tabloul 
```math
\left(\begin{array}{cc} a_{11} & a_{12} \\ a_{21} & a_{22} \end{array}\right),
```
numit matricea lui $f$ și 
```math
f(\overrightarrow{x})=\left(\begin{array}{cc} a_{11} & a_{12} \\ a_{21} & a_{22} \end{array}\right)
\left(\begin{array}{c} x_{1} \\ x_{2} \end{array}\right).
```

Atunci cînd considerăm o altă bază cu matricea de trecere $\Lambda$, atunci matricea lui ``f`` în raport cu noua bază va fi ``\Lambda^{-1} A \Lambda``.

## Un model matricial de populație

Revenim la modele privind ursul american (Ursus Americanus)
"""

# ╔═╡ 957700f1-603d-4d07-83ee-2cf788e0eac4
PlutoUI.LocalResource("./Black_Bear.jpg", :width => 1000)

# ╔═╡ 2f5ba151-7278-4947-9003-04620344aa76
md"""
este răspîndit pe întreg teritoriul SUA. Femelele ating maturitatea sexuală în jurul vîrstei de 3-4 ani și
trăiesc în sălbăticie intre 15-20 de ani. Fiecare femelă va da naștere o dată la 2 ani și în general va avea 2 pui. Vom împărți populația de urși
în două clase: juvenili (pui și animale care nu au atins maturitatea sexuală) și adulți, pe care le vom organiza ca un 2-vector
```math
\left(\begin{array}{c} J \\ A \end{array}\right).
```
În medie o femelă va naște un pui pe an. Ca să simplificăm discuția vom studia doar femelele și deci putem presupune că o femelă va naște $0,5$ pui
pe an. în fiecare an 10% din juvenili mor și 25% se maturizează la adulți, deci 65% vor rămîne juvenili. Dacă notăm cu ``J_{N}`` numărul de
juvenili la momentul ``N`` și cu ``A_N`` numărul de adulți la momentul ``N``, vom avea că ``J_{N+1}=0.65 J_N + 0.5 A_N.`` Speranța de viață a unui adult
este de aproximativ 14 ani și este adult pentru aproximativ 10 ani. Deci rata de deces per capita este ``1/10``, prin urmare 90% din adulți rămîn
adulți și în următorul an, de asemenea după cum am văzut 25% din juvenili devin adulți. Obținem în anul ``N+1`` o populație de adulți
``A_{N+1}=0.25 J_N + 0.9 A_N.`` Modelul de evoluție al populației de urși negri este dat de o funcție liniară
```math
\left(\begin{array}{c} J_{N+1} \\ A_{N+1} \end{array}\right)= A \cdot \left(\begin{array}{c} J_{N} \\ A_{N} \end{array}\right),
```
cu
```math
A=\left(\begin{array}{cc} 0.65 & 0.5  \\ 0.25 & 0.9 \end{array}\right).
```
Deci dacă știm că într-un an sînt 100 de juvenili și 50 de adulți, atunci în anul următor vom avea 90 de juvenili și 70 de adulți.

## Compunerea funcțiilor. Produsul matricelor

Funcțiile pot fi aplicate în mod înlănțuit. Am văzut că putem face acest lucru pentru funcții $\mathbb{R} \to \mathbb{R}.$ Aceeași construcție
funcționează și în dimensiuni mai mari. Dacă $f:\mathbb{R}^n \to \mathbb{R}^k$ și $g:\mathbb{R}^k \to \mathbb{R}^p$, atunci putem să construim
$g\circ f (\overrightarrow{x})=g(f(\overrightarrow{x})).$ Ne vom concentra mai ales pe cazul în care toate spațiile sînt $\mathbb{R}^n$.
Compunerea a două aplicații liniare $f,g :\mathbb{R}^n \to \mathbb{R}^n$ este tot o aplicație liniară $h:\mathbb{R}^n \to \mathbb{R}^n.$ Dacă $f$
și $g$ au matricele $A$, respectiv $B$ în raport cu baza canonică, vrem să determinăm matricea $C$ a lui $h$. Coloanele lui $C$ sînt valorile lui
$h$ pentru elementele bazei canonice. Prima coloană este $h(\overrightarrow{e_1})$. Avem $h(\overrightarrow{e_1})=f(g(\overrightarrow{e_1}))$,
$g(\overrightarrow{e_1})=b_{11}\overrightarrow{e_1}+b_{21}\overrightarrow{e_2}+\dots+b_{n1}\overrightarrow{e_n}$. Aplicăm acum funcția $f$,
ținem seama de faptul că aceasta este liniară. Deci prima coloană a lui $C$ se obține înmulțind $g(\overrightarrow{e_1})$ cu $A$ și la fel pentru
restul coloanelor. Deci $C=A\cdot B,$ produsul matricelor.

Revenim la studiul urșilor. Presupunem că este un an prost și natalitatea scade la 40%, mortalitatea juvenililor crește la 40% și doar 10% devin
adulți și mortalitatea adulților este 20%, deci
```math
\left(\begin{array}{c} J_{N+1} \\ A_{N+1} \end{array}\right)= A \cdot \left(\begin{array}{c} J_{N} \\ A_{N} \end{array}\right),
```
cu
```math
A=\left(\begin{array}{cc} 0.5 & 0.4  \\ 0.1 & 0.8 \end{array}\right).
```
Deci dacă un an prost vine după unul bun populația de urși va fi
```math
\left(\begin{array}{cc} 0.5 & 0.4  \\ 0.1 & 0.8 \end{array}\right)
\left(\begin{array}{cc} 0.65 & 0.5  \\ 0.25 & 0.9 \end{array}\right) \left(\begin{array}{c} 100 \\ 50 \end{array}\right)\approx
\left(\begin{array}{c} 70 \\ 62 \end{array}\right).
```
"""

# ╔═╡ 25da8b84-53c0-4337-9169-27f4b63efaf3
begin
A=[0.65 0.5; 0.25 0.9]
B=[0.5 0.4; 0.1 0.8]
init=[100,50]
C=A*B
C*init
end

# ╔═╡ 2091aec9-7c11-4df4-a02b-a302753b34d1
md"""

## Comportarea pe termen lung a modelelor matriciale
Modelele de tipul celui prezentat se numesc modele matriciale. În particular aceste tipuri de modele pentru populații segregate pe tipuri de
vîrstă se numesc modele *Leslie*. Pornind cu o condiție inițială $\overrightarrow{x}$, atunci
```math
\overrightarrow{x}_N=A^N\overrightarrow{x},
```
cu
```math
A=\left(\begin{array}{cc} 0.65 & 0.5  \\ 0.25 & 0.9 \end{array}\right).
```
Pe diagonala lui $A$ se află proporția din populația care rămîne în aceeași grupă de vîrstă, iar valorile care nu sînt pe diagonală ne dau
rata de tranziție și rata de deces. Să vedem comportamentul pe termen lung.
"""

# ╔═╡ d4e88ef9-2b2a-4b4a-bac9-c1008ac2070a
begin
	x=[50,10]
v_j=[]
v_a=[]
n_points=16
count=range(1,n_points,length=n_points)
for n ∈ count
    push!(v_j,x[1])
    push!(v_a,x[2])
    x=A*x
end
fig1=Figure()
    ax = Axis(fig1[1,1])
    scatter!(count,v_j, label="juvenili")
    scatter!(count,v_a, label="adulți")
    axislegend(ax, position = :rt)
fig1
end

# ╔═╡ 1a52a77b-ed00-47de-8b9e-9e61f07d4705
md"""
Se observă că după cîteva iterații numărul de juvenili devine egal cu cel de adulți și ambele cresc. Pentru anii proști avem
```math
A=\left(\begin{array}{cc} 0.5 & 0.4  \\ 0.1 & 0.8 \end{array}\right).
```
"""

# ╔═╡ a2d39229-dce8-4a2e-8690-c6c628fea06e
begin
	let
	x=[10,50]
A=[0.5 0.4; 0.1 0.8]
v_j=[]
v_a=[]
n_points=16
count=range(1,n_points,length=n_points)
for n ∈ count
    push!(v_j,x[1])
    push!(v_a,x[2])
    x=A*x
end
fig2=Figure()
    ax = Axis(fig2[1,1])
    scatter!(count,v_j,label="juvenili")
    scatter!(count,v_a,label="adulți")
    axislegend(ax, position = :rt)
fig2
	end
end

# ╔═╡ beb127e3-7cbd-420c-8948-d922de7113c0
md"""
De asemenea valorile pentru juvenili și pentru adulți devin egale, dar ambele scad. La început numărul de juvenili crește, înainte de a scădea,
Un astfel de comportament se numește *tranzient*. Să luăm acum un alt model Leslie
```math
A=\left(\begin{array}{cc} 0.1 & 1.4  \\ 0.4 & 0.2 \end{array}\right).
```

"""

# ╔═╡ 598082cd-f79a-4ac4-bb62-e11e20460dad
begin
	let
		x=[10,50]
A=[0.1 1.4; 0.4 0.2]
v_j=[]
v_a=[]
n_points=30
count=range(1,n_points,length=n_points)
for n ∈ count
    push!(v_j,x[1])
    push!(v_a,x[2])
    x=A*x
end
fig3=Figure()
    ax = Axis(fig3[1,1])
    scatter!(count,v_j,label="juvenili")
    scatter!(count,v_a,label="adulți")
    axislegend(ax, position = :rt)
fig3
	end
end

# ╔═╡ 00e29c5f-8466-46df-b56f-b11d066fe3ec
md"""
## Vectori și valori proprii
Pînă acum am văzut că un sistem liniar discret poate avea diferite comportamente, dar am observat acest lucru empiric, iterînd matricea de un număr mare de ori. Vom vedea cum putem determina riguros comportamentul. Pentru aceasta avem nevoie să studiem valorile și vectorii proprii.

Să luăm pentru început o aplicație de două variabile formată din două funcții de o singură variabilă:
```math
\left(\begin{array}{cc} u \\ v\end{array}\right) =f \left(\begin{array}{cc} x \\ y\end{array}\right) = \left(\begin{array}{cc} a_{11}x+a_{12}y \\ a_{21}x+a_{22}y\end{array}\right)=A\cdot\left(\begin{array}{cc} x \\ y\end{array}\right),
```
unde 
```math
A=\left(\begin{array}{cc} a_{11} & a_{12} \\ a_{21} & a_{22} \end{array}\right). 
```
Un caz foarte simplu este atunci cînd 
```math
\left(\begin{array}{cc} a & 0 \\ 0 & b\end{array}\right).
```
Ideal ar fi să putem aduce orice aplicație liniară la forma aceasta, pentru a putea studia mai ușor acțiunea sa. Pornim la drum cu o aplicație liniară ``f:\mathbb{R}^2 \to \mathbb{R}^2`` și căutăm doi vectori nenuli pentru care ``f(v_1) = av_1`` și respectiv
``f(v_2)=bv_2``. Dacă reușim să facem aceasta spunem că ``v_1`` și ``v_2`` sînt **vectori proprii**, iar numerele ``a`` și ``b`` se numesc **valori proprii**. Se poate vedea că pentru a obține astfel de vectori trebuie ca ``a`` și ``b`` trebuie să satisfacă 
**ecuația caracteristică** ``\lambda^2-(x_{11}+x_{22})\lambda+(x_{11}x_{22}-x_{12}x_{21}).`` Dacă această ecuație are două soluții distincte atunci pentru fiecare valoare proprie vom avea cîte un vector propriu, iar cei doi vectori vor fi liniar-independenți,
adică vor forma o bază.

Atunci dacă un vector are coordonatele ``(a_1,a_2)`` în raport cu această bază imaginea sa prin ``f`` va avea componentele ``(aa_1,ba_2).``

Revenind la exemplul cu urșii americani, considerăm sistemul discret 
```math
\left(\begin{array}{c} x_{n+1} \\ y_{n+1} \end{array}\right) = M \cdot \left(\begin{array}{c} x_{n} \\ y_{n} \end{array}\right) = M^{n+1} \cdot \left(\begin{array}{c} x_{0} \\ y_{0} \end{array}\right)
```
cu valorile inițiale ``(x_0,y_0).`` Presupunem că ``M`` este diagonalizabilă (adică are două valori proprii) și ``\Lambda`` este matricea de trecere la baza de vectori proprii. Atunci ``M^n=\Lambda D^n\Lambda^{-1}``, unde
```math
D=\left(\begin{array}{c} a & 0 \\ 0 & b \end{array}\right)
```
și deci
```math
D^n=\left(\begin{array}{c} a^n & 0 \\ 0 & b^n \end{array}\right).
```
"""

# ╔═╡ c1f23cef-817d-4cde-a641-12e76e199619
danger(md" Nu orice matrice este diagonalizabilă.")

# ╔═╡ adbd10e4-c7fc-4fdc-9f9d-bf1451c7cfb4
md"""
Aceasta se poate întîmpla din două motive. Pe de o parte s-ar putea ca ca ecuația caracteristică să nu aibă rădăcini reale sau să aibă o singură rădăcină. Cînd are o singură rădăcină atunci se poate întîmpla să aibă doar un vector propriu.

## Ce se întîmplă în dimensiune mai mare

Pentru matrice de ordin mai mare ecuația caracteristică este ``\det(\lambda I_n -M)=0`` unde ``I_n`` este matricea identică (are ``1`` pe diagonală și în rest ``0.``).

Pentru a calcula valorile proprii și vectorii proprii corespunzători folosim funcțiile julia analoage cu cele de mai jos:

"""

# ╔═╡ Cell order:
# ╟─1db058d6-c4b0-11f0-31d2-8fd829c839fa
# ╠═deef08bd-5eec-4e82-9940-6b53b2be99bd
# ╟─8c4c4fd7-0566-42be-ab06-5f448055b34e
# ╟─957700f1-603d-4d07-83ee-2cf788e0eac4
# ╟─2f5ba151-7278-4947-9003-04620344aa76
# ╟─25da8b84-53c0-4337-9169-27f4b63efaf3
# ╟─2091aec9-7c11-4df4-a02b-a302753b34d1
# ╟─d4e88ef9-2b2a-4b4a-bac9-c1008ac2070a
# ╟─1a52a77b-ed00-47de-8b9e-9e61f07d4705
# ╟─a2d39229-dce8-4a2e-8690-c6c628fea06e
# ╟─beb127e3-7cbd-420c-8948-d922de7113c0
# ╟─598082cd-f79a-4ac4-bb62-e11e20460dad
# ╟─00e29c5f-8466-46df-b56f-b11d066fe3ec
# ╟─c1f23cef-817d-4cde-a641-12e76e199619
# ╟─adbd10e4-c7fc-4fdc-9f9d-bf1451c7cfb4
# ╠═07efda97-46ed-4d2c-b347-d33e55561166
