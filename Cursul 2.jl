### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 2effecea-a41d-11f0-030d-0dea9a2cf25b
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	using DifferentialEquations
	TableOfContents()
end

# ╔═╡ f828c688-8ac5-45a4-8183-b89c16f51194
begin
	using FastDifferentiation
	
	@variables x1

	
	f3 = exp(sin(3*x1^2))
    g1 = make_function(derivative([f3],x1),[x1])
    derivative([f3],x1)
	f4=make_function([f3],[x1])
	t3 = range(-3,3,length=120)
    y3 = [f4(t)[1] for t ∈ t3]
    y4 = [g1(t)[1] for t ∈ t3]
	fig6 = Figure()
    ax6 = Axis(fig6[1,1])
    lines!(t3,y3)
    lines!(t3,y4)
    fig6
end

# ╔═╡ 1038b0f7-4372-4f1e-a498-f0ea0ebcae4f
html"""
<style>
	@media screen {
		main {
			margin: 0 auto;
			max-width: 2000px;
    		padding-left: max(283px, 10%);
    		padding-right: max(383px, 10%); 
            # 383px to accomodate TableOfContents(aside=true)
		}
	}
</style>
"""

# ╔═╡ e20606c5-ffa3-4d6b-81b8-b73df73e0305
md"""
# Cursul II

## Schimbarea în spațiul stărilor

După cum am văzut în cursul precedent sîntem interesați în principal de evoluția în timp a variabilelor de stare. Pentru aceasta avem nevoie să formulăm modele matematice. 

Scopul oricărui model este să înțelegem cauzele schimbării unui sistem în timp și să prezicem evoluția lui în viitor. Dacă $x$ este o variabilă de stare, vom nota cu $x'$ sau cu $\Delta x$ variația sa. De exemplu dacă $x$ este valoarea deținută într-un cont atunci $x'=d+r$, unde $d$ este suma depusă, iar $r$ este dobînda amîndouă măsurate în RON/lună.

Să presupunem acum că $x$ măsoară cantitatea de apă dintr-un bazin care se golește pe la fund. Atunci cantitatea de apă care se scurge într-un minut, să zicem, va depinde de cantitatea de apă din bazin. Să presupunem că viteza de curgere este proporțională cu cantitatea de apă. Avem 
```math
x'=-kx,
```
pentru o constantă *pozitivă* $k$. $x'$ se măsoară în $\frac{l}{min}$, $x$ în $l$, deci $k$ trebuie să fie măsurată în $\frac{1}{min}$.

În modelele pe care le vom studia, schimbarea variabilei de stare va depinde atît de valoarea sa, dar și de *parametrii*, constante care nu se schimbă în timp. 

Dacă vrem să modelăm cum se răcește o cană de cafea într-o cameră mai rece, începem prin a nota cu $T$ temperatura sa în $^\circ$K. Conform legii lui Newton rata de schimbare este proporțională cu diferența între $T$ și $r$, unde $r$ este temperatura camerei:
```math
T'=c\cdot(T-r),
```
unde $c$ este o constantă. Sistemul nostru va depinde de doi parametrii $r$ și $c$. Vrem să studiem mai atent pe $c$. Cum cafeaua este mai caldă decît camera avem $T-r> 0$. În timp cafeaua se va răci deci $T' < 0$ prin urmare $c <0$. La fel dacă este vorba de o cafea cu gheață care se încălzește, atunci $T-r < 0$, $T' > 0$ și iarăși $c<0$. Ecuația finală va fi
```math
T'=k(r-T),
```
unde $k=-c$.

## Cîteva modele biologice

### Un model simplificat de populație

Notăm cu $x$ numărul de animale dintr-o specie trăind pe un anumit teritoriu. Avem $x'=$rata nașterilor-rata deceselor.Trebuie acum să reprezentăm matematic cele două mărimi. Pentru început să presupunem că animalele nu mor (de exemplu dacă intervalul de timp este mult mai mic decît durata de viață). Facem o nouă presupunere simplificatoare: animalele nu au sex (oricare este capabil să dea naștere), toate animalele au aceeași fertilitate de-a lungul vieții și toate animalele au aceeași șansă să nască. Să presupunem că rata de naștere este $b=0.5$, adică un animal naște o dată la 2 ani. Deci
```math
x'=0.5x
```
($b$ este rata de nașteri per capita). Unitatea de măsură pentru $b$ este animal/animal/an=$\frac{1}{an}$. Să luăm un alt model în care animalele doar mor. Atunci $x'=-dx$, unde $d$ este rata de morți. Am presupus de asemenea că toate animalele trăiesc același număr de ani, rata de decese nu depinde de numărul de animale, rata de decese nu depinde de timp. Avem deci $x'=-dx$. 

Cele două modele simplificate pot fi unificate prin $x'=bx-dx=(b-d)x=rx$, unde $r$ este rata de creștere a populației. 

### Un model cu aglomerare
Evident modelul $x'=bx$ este complet nerealist. Luăm în considerare faptul că resursele sînt limitate, deci există competiție pentru ele. Prin urmare $x'=b\cdot f\cdot x$, unde $f<1$ este factorul de aglomerare. Să presupunem că mediul are o capacitate de încărcare $k$ (numărul maxim de animale pe care o poate susține). Atunci $\frac{x}{k}$ reprezintă resursele utilizate de populația curentă și deci $\left(1-\frac{x}{k}\right)$ sînt resursele rămase libere. Deci 
```math
x'=bx\left(1-\frac{x}{k}\right).
```
Aceasta este **ecuația logistică**, foarte importantă și cu care ne vom mai întîlni adesea.

### Arcuri
Să luăm acum un sistem format dintr-un arc și o greutate care se mișcă la capătul său. În această situație avem nevoie de 2 variabile, anume poziția (centrului de greutate) greutății $x$, dar și viteza sa $v$. Ecuațiile de schimbare vor fi
```math
\begin{align}
x' & = & f(x,v) \\
v' & = & g(x,v).
\end{align}
```
Prima ecuație este foarte simplă, dată de definiția vitezei: $x'=v$. Pentru $v'$ trebuie să aducem aminte de legea lui Newton $F=ma$, unde $a=v'$. Deci $v'=\frac{F}{m}$. Valoarea lui $F$ provine din legea lui Hooke $F=-kx$.Deci sistemul de ecuații este:
```math
\begin{align}
x' & = & v \\
v' & = & -\frac{k}{m} x.
\end{align}
```
Dacă ne alegem sistemul de unități astfel încît $\frac{k}{m}=1$, atunci sistemul ia forma mai simplă
```math
\begin{align}
x' & = & v \\
v' & = & -x.
\end{align}
```

### Rîși și iepuri

Notăm cu $x$ numărul de rîși și cu $y$ numărul de iepuri. Vrem să scriem $x'=$ și $y'=$. Ce schimbă numărul de rîși?
În primul rînd rîșii se nasc și mor. Presupunem că rata de deces $d$ este constantă, iar rata de nașteri este proporțională cu cantitatea de pradă. Notăm cu $m$ rata de proporționalitate. De cîte ori un rîs întîlnește un iepure probabilitatea să-l prindă este $\beta$. În concluzie $x'=m\beta xy - dx.$ Un raționament analog ne va da $y'=bx-\beta xy$. Sistemul obținut:
```math
\begin{align}
x' & = & m\beta xy & - & dx \\
y' & = & by & - & \beta xy
\end{align}
```
se numește sistemul *Lotka-Volterra*.

Pentru a studia calitativ aceste ecuații presupunem că toți parametrii sînt egali cu $1$:

```math
\begin{align}
x' & = & xy & - & x \\
y' & = & y & - & xy
\end{align}
```
care produc comportamentul oscilatoriu observat.

### Infecția cu HIV într-un individ

Sindromul imuno-deficienței dobîndite (SIDA) este produs de virusul HIV. HIV infectează o anumită clasă de limfocite T ($CD4^+$). După infectarea unei celule, fie aceasta începe imediat să producă viruși și moare în cîteva zile sau preia materialul genetic al virusului, rămîne aparent sănătoasă, dar se poate activa mai tîrziu. Cînd o persoană este infectată concentrația virusului crește foarte mult, dar după cîteva luni scade la un nivel mult mai mic. La început s-a crezut că aceasta se datorează răspunsului imunitar, dar un model din 1996 a arătat că aceasta este posibil și în absența oricărui răspuns imunitar.

Variabilele noastre vor fi $V$, cantitatea de viruși, $R$ numărul de celule neinfectate și $E$ numărul celor infectate. Virușii sînt produși de celule, apoi infectează noi celule sau mor. Numărul celor care infectează este mic, deci poate fi ignorat. Presupunem că fiecare celulă infectată produce $100$ de viruși pe zi, și rata de deces este $2$ (adică un virus trăiește în medie 12 ore), deci $$V'=100E-2V.$$ Celulele neinfectate sînt produse de corp, mor natural și pot deveni infectate. În medie sînt produse $0,272$ celule (per mm$^3$), rata de deces este $0,00136$. Pentru a fi infectată o celulă trebuie să întîlnească un virus, ca și în cazul rîșilor și iepurilor aceasta este proporțional cu produsul $RV$. Pe scurt 
```math
R'=0,272-0,00136R-0,00027RV.
```
Pentru celulele infectate avem
```math
E'=0,00027RV-0,33E.
```
Rata de mortalitate a celulelor infectate este $0,33$. Grupat sistemul (introdus de Andrew Phillips pornind de la lucrările Angelei McLean) este:
```math
\begin{eqnarray}
V'&=& 100E-2V \\
R'&=& 0,272-0,00136R-0,00027RV \\
E'&=& 0,00027RV-0,33E.
\end{eqnarray}
```
Acest model prezice scăderea numărului de viruși, dar nu prevede declanșarea bolii. Există versiuni îmbunătățite, dar mai dificile, care prezic apariția bolii.

## Reacții chimice

Reacțiile chimice se scriu de forma
```math
A+B \xrightarrow{k} C
```
Cantitățile de subtanțe se modifică în timp și vrem să putem scrie un sistem de ecuații pentru evoluția acestora. 

Ca să vedem cum putem deduce un astfel de sistem. Pornim cu o reacție mai simplă 
```math
A \xrightarrow{k} B
```

Aici $k$ este rata cu care $A$ se transformă în $B$. Deci concentrația lui $A$, notată tot cu $A$, variază după regula
```math
A'=-kA.
```
Luăm acum reacția mai complicată 
```math
2A \xrightarrow{k} B \text{ ,adică } A+A \xrightarrow{k} B.
```
Ca și în exemplul precedent concentrația lui A poate doar să scadă. Lege acțiunii masă ne spune că
```math
A'=-2kA^2.
```

Revenind la reacția
```math
A+B \xrightarrow{k} C
```
ovservăm că $A$ și $B$ scad întotdeauna, atunci cînd se întîlnesc. Deci
```math
A'=-kAB.
```
## Interpretarea geometrică a schimbărilor

Spațiul stărilor este mulțimea tuturor valorilor posibile pentru variabilele de stare $x$. Dar care este noțiunea similară pentru $x'$? Nu poate fi spațiul stărilor pentru că $x'$ se măsoară cu alte unități decît $x$. De asemenea de multe ori $x$ ia doar valori pozitive, dar $x'$ poate fi și negativ.

Spațiul unde iau valori variațiile $x'$ este spațiul tangent. Schimbarea poate fi gîndită ca o mișcare prin spațiul stărilor, deci și noi vom reprezenta elementele din spațiul tangent prin săgeți care are vîrful în direcția schimbării și lungimea indică mărimea schimbării. 

Prin urmare modelul este o ecuație diferențială care asociază fiecărui punct din spațiul stărilor cîte un vector tangent. Deci ecuația dinamică este o funcție, numită **cîmp de vectori.** În general vom desena cîmpurile de vectori peste reprezentarea spațiului stărilor. Cel mai interesant din punct de vedere grafic este cazul $2$ dimensional. 

Ne putem imagina că un punct din spațiul stărilor se deplasează, descriind o curbă. În fiecare punct viteza instantanee coincide cu valoarea în punctul respectiv a cîmpului de vectori. O astfel de curbă se numește **curbă integrală** sau **traiectorie**.

Pentru modelul arcului cîmpul de vectori asociat împreună cu traiectoriile asociate este
"""

# ╔═╡ 17aa5b94-b3be-45d2-b694-b86a290104ab
begin
	f1(x) = Point2f(-x[2],x[1])
	fig1 = Figure()
	ax1 = Axis(fig1[1,1])
	xs1 = LinRange(-2, 2, 10)
    ys1 = LinRange(-2, 2, 10)
	us1 = [-y for x ∈ xs1, y ∈ ys1]
	vs1 = [x for x ∈ xs1, y ∈ ys1]
	strength = vec(sqrt.(us1 .^ 2 .+ vs1 .^ 2))
	streamplot!(ax1,f1,xs1,ys1,colormap = :magma,linewidth = 1)
	fig1
end

# ╔═╡ 91a4c1d4-94b7-4950-815d-9a93b9b7479a
begin
	fig3 = Figure()
	ax3 = Axis(fig3[1,1])
	arrows!(ax3,xs1,ys1,f1,lengthscale = 0.3,color=strength)
	fig3
end

# ╔═╡ 52134e56-2475-43c8-b28e-a434a6989249
md"""
Pentru sistemul care modelează interacțiunea pradă-prădător (Lotka-Volterra)
"""

# ╔═╡ 4c6ab08d-3adb-403f-b45e-fa7f50896a4d
begin
	f2(x) = Point2f(x[1]-x[1]*x[2],x[1]*x[2]-x[2])
	fig2 = Figure()
	ax2 = Axis(fig2[1,1])
	xs2 = LinRange(0, 2, 15)
    ys2 = LinRange(0, 2, 15)
	us2 = [x-x*y for x ∈ xs2, y ∈ ys2]
	vs2 = [x*y-y for x ∈ xs2, y ∈ ys2]
	strength2 = vec(sqrt.(us2 .^ 2 .+ vs2 .^ 2))
	streamplot!(ax2,f2,xs2,ys2,colormap = :magma,linewidth = 1)
	fig2
end

# ╔═╡ e7f47501-4ac7-437b-aa02-e2e586972a31
begin
	fig4 = Figure()
	ax4 = Axis(fig4[1,1])
	arrows!(ax4,xs2,ys2,f2,lengthscale = 0.3,color=strength2)
	fig4
end

# ╔═╡ 601f61a9-7b09-4817-b7d0-89f93d5f7bc5
md"""
O reprezentare alternativă este aceea de a desena serii temporale (grafice) pentru fiecare componentă în parte. 

Pentru ecuațiile diferențiale cu care ne vom întîlni avem că pentru orice valoare inițială există și este unică o traiectorie. O problemă naturală este cum desenăm traiectoria pentru o anumită valoare inițială. O idee ar fi să urmăm vectorul pînă la capăt, dar acesta presupune că vom urma o unitate întreagă de timp să zicem un an. Dar ce se întîmplă după 6 luni? Deci următorul moment trebuie să fie la mai puțin de un an, mai puțin de 6 luni, etc. Cum între orice două numere reale se află o infinitate de numere practic nu este posibil să găsim un moment "următor".

Oricît am urmări vectorul de schimbare am greși. Deci trebuie să schimbăm direcția infinit de repede. Metoda prin care putem desena traiectoria este

## Metoda Euler
Să presupunem că timpul în care urmăm săgeata este $\Delta t$ (acest timp este foarte scurt). Ideea este să presupunem că $\Delta t$ devine din ce în ce mai mic. Vom vedea în cursurile următoare cum vom putea face riguroasă această construcție.

Metoda Euler presupune că $\Delta t$ este mic, dar nu $0$. Presupunem că $x'=f(x)$ este o ecuație diferențială cu o condiție inițială $x_0$. Schimbarea lui $x_0$ este $f(x_0)$, dar urmăm vectorul de schimbare doar $\Delta t$ și deci schimbarea adevărată este $\Delta t\cdot f(x_0)$. Obținem deci noua valoare
```math
x_1=x_0+\Delta t \cdot f(x_0).
```
Următorul pas este acum
```math
x_2 = x_1+\Delta t \cdot f(x_1).
```
Spre exemplu dacă avem ecuația $x'=0,5x,$ $x_0=100$ și $\Delta t=0,001$, atunci $x_1=100+0.001\cdot 50 = 100.5$, $x_2=100.5+0,001\cdot 50,25=100,55025$, etc.

În cazul multi-dimensional procedăm oarecum analog. Să luăm spre exemplu sistemul Lotka-Volterra
```math
\begin{eqnarray}
x' & = & xy - x \\
y' & = & -xy + y
\end{eqnarray}
```
cu valorile inițiale $(x_0,y_0)=(2,3)$ și fie $\Delta t=0,1.$ Atunci $x_1=x_0+\Delta t\cdot (6-2)=2,4$ și $y_1=y_0+\Delta t\cdot(-6+3)=2,7$.

# Derivate 

Am văzut că atunci cînd aplicăm metoda lui Euler, momentul următor este $t+\Delta t$ și noua valoare a variabilei de stare este $X(t+\Delta t)$, aveam $X(t+\Delta t)=X(t) + \Delta t \cdot X'.$ O ilustrare este mai jos:


"""

# ╔═╡ 3e0bc92d-b1d9-4d5b-b0ce-a066f13efcd6
begin
	t=range(0,10,length=15)
    y=cos.(t)
    t2 = range(0,10,length=150)
    y2 = cos.(t2)
    fig5 = Figure()
	ax5 = Axis(fig5[1,1])
    scatter!(t,y)
	lines!(t,y)
	lines!(t2,y2)
	fig5
end

# ╔═╡ 277dd939-0ca4-42c5-b36a-6024d8d9f350
md"""
De fapt $X(t+\Delta t) \approx X(t) + \Delta t\cdot X'$ pentru că nu avem egalitate ci doar aproximare. Prin urmare

```math
X' \approx \frac {X(t+\Delta t)-X(t)}{\Delta t}.
```

În cuvinte variația variabilei este aproximată de diferența între valorile acesteia la două momente apropiate raportat la intervalul de timp. Putem să ne gîndim la următoarea analogie:

Dacă ne deplasăm din punctul $A$ în punctul $B$ și parcurgem o distanță de 10 km în 30 de minute (1/2 ore) atunci viteze medie va fi $$\frac{10}{0.5}=20 \text{ km/h.}$$ Să notăm acum cu $x(t)$ distanța pe care am parcurs-o în timpul $t$. Atunci viteza medie pe intervalul $(t_1,t_2)$ este 
```math
\frac{x(t_2)-x(t_1)}{t_2-t_1}=:\frac{\Delta X}{\Delta t},
```
dacă alegem un interval mai mic $(t_1,t_3)$ atît $\Delta t$ cît și $\Delta x$ vor fi mai mici. Viteza medie este evident că poate fi calculată pe un interval oricît de mic. 

Ce facem însă pentru a calcula viteza instantanee la un moment dat? Dacă vrem să folosim viteza medie atunci am avea că 
```math
v = \frac{0}{0}
```
ceea ce este imposibil. Rezolvarea este să ne uităm la toate vitezele medii pentru intervale de tipul $t_0+\Delta t$ cînd $\Delta t$ este din ce în ce mai mic, adică să facem $\Delta t \to 0$. Dacă în felul acesta obținem un număr, atunci acesta este viteza instantanee.
```math
v_{t_0}=\lim_{\Delta t \to 0} \frac{\Delta x}{\Delta t}=\lim_{\Delta t \to 0}\frac{x(t_0+\Delta t)-x(t_0)}{\Delta t}.
```
Putem acum să definim rata de variație a lui $x$ în $t$ ca "viteza instantanee" a lui $x$ în $t$:
```math
x'(t)=\lim_{\Delta t \to 0}\frac{x(t+\Delta t)-x(t)}{\Delta t},
```
numită **derivata** lui $x$ în $t$.

Galileo a descoperit legea căderii corpurilor: Dacă notăm cu $H(t)$ distanța corpului pînă la sol, atunci 
```math
H(t)=H(0)-4,9t^2.
```
Să calculăm viteza de cădere la momentul $t=1,5$s după ce a fost aruncat de la $100$ m. Calculăm mai întîi vitezele medii pentru valori ale lui $\Delta t$ din ce în ce mai mici:
"""

# ╔═╡ 47c525f6-ca11-4d3c-a287-5aed70ad7549
begin
H0=100
H(t)=H0-4.9*t^2
dt=0.1
v1=(H(1.6)-H(1.5))/0.1
println(v1)
v2 = (H(1.51)-H(1.5))/0.01
println(v2)
v3 = (H(1.501)-H(1.5))/0.001
println(v3)
println(H(1.5))
4.9*(1.5^2)
end

# ╔═╡ 1394f8f4-d425-4522-aa0f-ae55964d30ce
md"""
Se pare că ne apropiem de $-14,7$, dar să vedem acest lucru riguros. Mai întîi calculăm $H(1,5)=88.975$, apoi 
```math
\begin{aligned}
H(1,5+\Delta t)=100-4,9(1,5+\Delta t)^2=100-(11,025+14,7\Delta t+(\Delta t^2))=\\
=88,975-14,7\Delta t-\Delta t^2.
\end{aligned}
```
Derivata este
```math
H'(1,5)=\lim_{\Delta t \to 0} \frac{88,975-14,7\Delta t-\Delta t^2-88,975}{\Delta t}=-14,7. 
```
Pentru fiecare valoare a lui $t$ procedăm similar și avem
```math
\begin{align}
H'(t) = \lim_{\Delta t \to 0} \frac{H(t+\Delta t)-H(t)}{\Delta t}=\frac{100-4,9(t^2+2t\Delta t+\Delta^2)-100+4,9t^2}{\Delta t}= \\ 
=\lim_{\Delta t \to 0} (-9,8t-\Delta t)=-9,8t.
\end{align}
```
**Notație:** Derivata lui $x$ în $t_0$ se va nota $x'(t_0)$ sau 
```math
\left.\frac{dx}{dt}\right|_{t_0}.
```
Notația a doua provine din definiția echivalentă a derivatei:
```math
\lim_{t \to t_0}\frac{x(t)-x(t_0)}{t-t_0}=\lim_{t \to t_0} \frac{\Delta x}{\Delta t}.
```

## Interpretarea geometrică

Putem studia derivata și într-o manieră geometrică. Să luăm o funcție $Y(X)$. Atunci rata medie de schimbare este
```math
\left. \frac{\Delta Y}{\Delta X}\right|_{X_0}=\frac{Y_2-Y_1}{X_2-X_1}.
```
Dacă ne uităm la graficul funcției atunci $\Delta Y$ este modificarea verticală, iar $\Delta X$ este modificarea orizontală. Dacă desenăm secanta între punctele $(Y_1,X_1)$ și $(Y_2,X_2)$, panta acestei drepte este chiar $\frac{\Delta Y}{\Delta X}$. Deci rata medie de schimbare este totuna cu panta secantei. Cînd $X_2$ se tot apropie de $X_1$ atunci secanta se apropie din ce în ce mai mult de curbă, astfel că la final va atinge curba doar într-un punct. Această dreaptă se numește **tangenta** la curbă. Deci **panta tangentei** în $X_1$ este **derivata** $\left.\frac{dY}{dX}\right|_{X_1}$.

Ecuația generală a unei drepte este $Y=mX+b$. Vrem să aflăm $m$ și $b$ pentru tangenta în $X_1$. Pentru $b$, ținem cont că $Y_1=mX_1+b$, deci $b=Y_1-mX_1$. Întroducem înapoi în ecuație și obținem că $Y-Y_1=m(X-X_1)$. Ne aducem aminte cine este $m$ și obținem în final
```math
Y-Y_1=\left.\frac{dY}{dX}\right|_{X_1}(X-X_1),
```
adică $Y-Y_1$ depinde liniar de $X-X_1$.

O functie $f$ pentru care 
```math
\begin{aligned}
f(X_1+X_2) & = & f(X_1)+f(X_2) \\
f(aX) & = & a f(X)
\end{aligned}
```
se numește **liniară**. Rescriind $Y=f(X)$ vedem că derivata produce aproximarea funcției cu o funcție liniară.

**Atenție: Nu orice funcție este derivabilă**.

## Derivata unei funcții

Am văzut că dacă o funcție $f$ este derivabilă într-un punct, atunci îi putem atașa un număr, derivata sa în acel punct.Variind punctul obținem o nouă funcție pe care o vom nota
```math
\frac{df}{dx}.
```
Explicit derivata funcției este definită astfel:
```math
\frac{df}{dx}(x)=\lim_{\Delta x \to 0} \frac{f(x+\Delta x)-f(x)}{\Delta x}.
```
Dacă de exemplu $f(x)=x^2-x$, atunci 
```math
\begin{aligned}
\frac{df}{dx}=\lim_{\Delta x \to 0}\frac{(x+\Delta x)^2-(x+\Delta x)-x^2+x}{\Delta x}= \\ 
\lim_{\Delta x \to 0} \frac{2x\Delta x+\Delta x^2-\Delta x}{\Delta x}=\lim_{\Delta x \to 0} 2x+\Delta x-1=2x-1.
\end{aligned}
```

Cum derivata este o funcție putem repeta construcția și obținem derivatele de ordin superior. Notațiile sînt:
```math
\frac{d^nf}{dx^n}
```
sau $f'',f''',$ etc. De exemplu dacă $H(t)$ este înălțimea, atunci $H'(t)$ este viteza, $H''(t)$ este accelerația, $H'''(t)$ este supraaceleratia (jerk).

Pentru cele mai uzuale funcții derivatele sînt:

 - $\frac{d}{dx} c = 0,$ unde $c$ este o constantă;
 - $\frac{d}{dx} x^n = nx^{n-1},$ unde $n$ poate fi orice număr real;
 - $\frac{d}{dx} e^{kx}=ke^{kx},$ unde $e$ este baza logaritmilor naturali și $k$ este orice număr;
 - $\frac{d}{dx} ln(x)=\frac{1}{x};$
 - $\frac{d}{dx} \sin(x) = \cos(x);$
 - $\frac{d}{dx} \cos(x) = -\sin(x).$
 
 Funcțiile pe care le întîlnim în practică sînt obținute din cele uzuale folosind cîteva operații. Este important să vedem cum se comportă derivata referitor la aceste operații:
 
 - $\frac{d (c\cdot f)}{dx}=c\cdot\frac{df}{dx};$
 - $\frac{d(f+g)}{dx}=\frac{df}{dx}+\frac{dg}{dx};$
 - $\frac{d(f\cdot g)}{dx}=\frac{df}{dx}\cdot g+f\cdot\frac{dg}{dx};$
 - $\frac{d\frac{f}{g}}{dx}=\frac{1}{g^2}\left(\frac{df}{dx}\cdot g - f\cdot\frac{dg}{dx}\right);$
 - $\frac{d (f\circ g)}{dx}=\frac{df}{dg}\cdot\frac{dg}{dx}.$
 
 La ultima relație $\frac{df}{dg}$ semnifică faptul că ne referim la $f$ ca o funcție doar de $g$, altfel spus punem $y=g(x)$ și considerăm $f(y)$.
"""

# ╔═╡ 657e2e23-4f22-4822-855d-e247e6c7e98f
md"""
## Derivata implicită

Adesea avem nevoie să calculăm derivata unei funcții, dar nu avem o expresie explicită a acestei funcții, ci doar o ecuație pe care o satisface funcția. Spre exemplu modelul Kermack-McKendrick de propagare a unei maladii infecțioase conduce la expresia

```math
\rho e^{-qA} = 1-A,
```
unde $A$ este partea (procentul) de populație infectată, $q$ măsoară transmisibilitatea și $\rho$ este procentul de populație susceptibilă la infecție. Studiem cum se modifică $A$ atunci cînd susceptibilitatea crește. Pentru aceasta trebuie să determinăm derivata lui $A$, chiar dacă nu avem o formulă 
explicită a sa. Putem totuși deriva în ambii membri:

```math
\frac{d}{dq}\left(\rho e^{-qA}\right)=\frac{d}{dq}(1-A)
```

și obținem

```math
-\frac{dA}{dq}=\rho e^{-qA} \left(-q\frac{dA}{dq}-A\right).
```

Rezolvăm pentru $\frac{dA}{dq}$ și obținem 

```math
\frac{dA}{dq}=\frac{A\rho e^{-qA}}{1-pqe^{-qA}}=\frac{A\rho}{e^{qA}-pq}.
```
"""

# ╔═╡ Cell order:
# ╟─2effecea-a41d-11f0-030d-0dea9a2cf25b
# ╠═1038b0f7-4372-4f1e-a498-f0ea0ebcae4f
# ╟─e20606c5-ffa3-4d6b-81b8-b73df73e0305
# ╟─17aa5b94-b3be-45d2-b694-b86a290104ab
# ╟─91a4c1d4-94b7-4950-815d-9a93b9b7479a
# ╟─52134e56-2475-43c8-b28e-a434a6989249
# ╟─4c6ab08d-3adb-403f-b45e-fa7f50896a4d
# ╟─e7f47501-4ac7-437b-aa02-e2e586972a31
# ╟─601f61a9-7b09-4817-b7d0-89f93d5f7bc5
# ╟─3e0bc92d-b1d9-4d5b-b0ce-a066f13efcd6
# ╟─277dd939-0ca4-42c5-b36a-6024d8d9f350
# ╟─47c525f6-ca11-4d3c-a287-5aed70ad7549
# ╟─1394f8f4-d425-4522-aa0f-ae55964d30ce
# ╟─f828c688-8ac5-45a4-8183-b89c16f51194
# ╟─657e2e23-4f22-4822-855d-e247e6c7e98f
