### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 9e1db0e6-d997-11f0-348a-d9f585410f74
begin
	using Pkg
	Pkg.activate(".")
	using Makie, CairoMakie
	using PlutoUI
	
	using PlutoTeachingTools
	
	# using LoopVectorization
	TableOfContents()
end

# ╔═╡ f7a6e388-7d3e-4eb6-8788-bb13dad71898
begin
using StatsBase
lst=[5,0,0,54,5,12,27,24,36,5,56,43,15,42,12,62,36,34,58,23]
	println("Datele sînt: ",lst)
println("Media este: ", mean(lst))
println("Mediana: ", median(lst))
println("Modulul este: ", mode(lst))
println("Deviația standard este: ",std(lst))
println("Prima quantilă: ",percentile(lst,25))
println("A treia quantilă: ",percentile(lst,75))
println("IQR = ",iqr(lst))
end

# ╔═╡ 51f4e264-f290-46f1-838d-db3a4dd74ba6
begin
	using StatsPlots
data = [5,0,0,54,5,12,27,24,36,5,56,43,15,42,12,62,36,34,58,23]
data_cuc = [19 19 20 22 22 22 23 24 25 35]
# data = reshape(lst, length(lst),1)
StatsPlots.boxplot([" "],data,leg=false,fillalpha=0.75;orientation = :v,ylim=(minimum(data)-1,maximum(data)+1))
dotplot!([" "],data,leg=false,fillalpha=0.75;ylim=(minimum(data)-1,maximum(data)+1)) 
# StatsPlots.plot(bp,dp)
end

# ╔═╡ dfb42990-8e1a-4b3e-ad1d-11054d7b0c1d
begin
	using Plots
locatii = ["col uterin", "gură", "organe genitale"]
poz = [1,2,3]
counts = [500000, 10000, 50000]
bar(locatii,counts,leg=false, title="Localizarea cancerelor induse de HPV în SUA")
end

# ╔═╡ 9f0d0ed7-1ed4-469d-8bcd-80f89965afb7
WidthOverDocs()

# ╔═╡ bc79a3dd-ccc3-43be-b90d-048cf6059b4f
md"""
# Elemente de statistică descriptivă

Statistica descriptivă este o colecție de tehnici care ne permit să extragem informații despre datele numerice pe care le 
obținem în urma unei serii de experimente sau unei serii de măsurători. Obiectele care sînt măsurate se numesc **indivizi**, iar 
proprietățile care sînt măsurate se numesc **variabile**. 

Variabilele sînt de mai multe tipuri:
 - **categoriale**: descriu apartenența la diferite categorii. De exemplu notele, genotipul, prezența sau absența unei boli. 
   Dacă valorile au o ordine naturală, atunci variabilele se numesc **ordinale**. Altfel se numesc **nominale**;
 - **numerice**: variabilele iau valori numerice care pot fi ordonate și asupra cărora are sens să facem operații aritmetice. 
   Variabilelele numerice pot fi **discrete** sau **continue**, în funcție de valorile posibile. Variabilele discrete iau valori
   într-o submulțime discretă (de exemplu $\mathbb{N}$ sau $\mathbb{Z}$) iar cele continue într-un interval de numere reale.
   
Un exemplu de variabile categoriale apare în experimentele lui Gregor Mendel cu plantele de mazăre. Una dintre trăsăturile
urmărite era poziția florilor: axială sau terminală. Cum nu avem nici o posibilitate de ordonare avem variabile nominale.

Atunci cînd studiem variabile numerice avem multe valori posibile și multe structuri care pot apărea. Ca să studiem astfel de 
variabile o primă mărime care ne interesează este "mijlocul" mulțimii de valori. Avem mai multe metode de identificare a 
acestuia:
  - Fie $x_1,\dots,x_n$ o mulțime de valori. **Media** lor este
    ```math
    \overline{x}=\frac{x_1+\dots+x_n}{n}.
    ```
  - O altă măsură a centralității este **mediana**. Dacă numărul de variabile este impar este egală cu valoarea din mijloc, iar
    dacă numărul este par este media aritmetică a valorilor din mijloc.
  - O a treia măsură a tendinței centrale este **modulul**. Acesta este valoarea variabilei care apare cel mai des.

Spre exemplu luăm ca mulțime a valorilor cantitatea de paraziți pe un număr de libelule $\left\{5,0,0,54,5,12,27,24,36,5,56,43,
15,42,12,62,36,34,58,23\right\}$. Atunci media este $27,45$, mediana este $25,5$ și modulul este $5$.

Mărimile "centrale" dau valorile spre care se concentrează variabilele. Un alt aspect important este **variabilitatea**, 
adică felul în care valorile se îndepărtează de medie.  Spre exemplu ambele mulțimi de valori:
```math
\begin{aligned}
\left\{50, 58, 78, 81, 93\right\} \\
\left\{72, 71, 72, 72, 73\right\}
\end{aligned}
```
au media $72$, dar în mod clar prima este mult mai dispersată decît a doua. Măsurile dispersiei sînt:
 - **deviația standard**: Dacă $x_1,\dots,x_n$ este o mulțime de valori cu media $\overline{x}$, atunci
   ```math
   d.s. = \sqrt{\frac{1}{n}\left(\left(\overline{x}-x_1\right)^2+\dots+\left(\overline{x}-x_n\right)^2\right)}.   
   ```
 - **varianța** este pătratul deviației standard.
 
Un alt indicator al variabilității este dat de valorile minimă și maximă. De asemenea putem lua mediana primei jumătăți a 
mulțimii numită **primul sfert** ($Q_1$), precum și mijlocul celei de-a doua jumătăți, numită **al treilea sfert** ($Q_3$).

Avem sumarul cu $5-$numere a variabilelor: min, $Q_1$, mediana, $Q_3$, max. De la ele putem calcula cîteva măsuri ale întinderii
datelor: întinderea datelor este $r=max-min$. De asemenea avem întinderea între sferturi (IQR) $IQR=Q_3-Q_1$. 

Adesea datele colectate în științele vieții conțin date **aberante** (outliers), adică date foarte îndepărtate de mediană. Acestea pot
proveni din erori de măsurare sau să fie datorate variabilității naturale. Se admite că o valoare este aberantă dacă este la cu mai mult decît $1.5\cdot IQR$ deasupra lui $Q_3$ sau mai mult decît $1.5\cdot IQR$ sub $Q_1.$
"""

# ╔═╡ 13dfb805-fdb6-4e4a-93c1-154303f4b331
md"""
 Datele pot fi reprezentate grafic cu ajutorul unui box plot precum mai jos, unde vedem cele cinci numere care caracterizează 
 setul nostru: 
"""

# ╔═╡ faefd827-a505-4f0d-b338-d443860a0c17
md"""
Pentru datele categoriale avem de asemenea mai multe posibilități de reprezentare. Primul este cu ajutorul unui **grafic de 
bare**. Avem cîte o bară pentru fiecare categorie, iar lungimea sa este proporțională cu numărul de indivizi din fiecare 
categorie. Un exemplu pentru localizarea cancerelor HPV este:
"""

# ╔═╡ dbebc2ae-a819-447c-9ce8-a0d63d5eb673
md"""
O altă metodă de vizualizare este **graficul proporților (pie chart)**. Aceasta este un cerc împărțit în sectoare, cîte unul pentru fiecare categorie. Unghiul central și deci aria fiecărui sector este proporțional cu numărul de indivizi din fiecare 
categorie. De exemplu datele pentru mortalitatea SARS în perioada 2002-2003 pot fi prezentate astfel:
"""

# ╔═╡ f163a685-b7f1-49b2-bf25-6edb98904eb4
begin
	sizes = [43, 349, 299, 37,33]
labels = ["Canada","China","Hong Kong","Taiwan","Singapore"]
#fig1 = Figure()
#ax = Axis(fig1[1,1])
Plots.pie(labels,sizes,l=0.5)
end

# ╔═╡ 7001202f-936f-4982-beb5-31627d20ba68
md"""
Ca să vizualizăm variabilele numerice, mai întîi transformăm datele în variabile categorice (mai precis ordinare) plasînd valorile între anumite limite și numărînd numărul de indivizi din fiecare categorie. Vizualizare însăși se face cu ajutorul **histogramelor**. O histogramă este formată din dreptunghiuri 
cu baze egale, numite **cutii**. Înălțimea fiecărui dreptunghi este proporțională cu numărul de indivizi din fiecare cutie. 
Valoarea concretă a înălțimii este determinată în mai multe situații:
  - pentru **histograme de frecvență** înălțimea este egală cu numărul de inidivizi în fiecare cutie;
  - pentru **histograme de frecvență relativă** înălțimea este proporțională cu numărul de indivizi din cutie;
  - pentru **histograme de densitate** înălțimea este aleasă astfel încît aria dreptunghiului este proporțională cu numărul de indivizi din fiecare cutie.

Indiferent de convenția folosită, histogramele au aceeași formă. În general nu există nici un criteriu clar pentru a calcula 
mărimea unei cutii, dar există anumite reguli. Printre ele regula lui Sturges: calculăm numărul natural $k$, cel mai aproape de 
$1+\log_2n$, unde $n$ este numărul de puncte. Atunci dacă $R$ este întinderea datelor, luăm numărul de cutii $R/k$.

Histograma de densitate este adesea cea mai convenabilă, astfel mediana este localizată la valoarea care împarte aria totală
a histogramei în jumătate.
"""

# ╔═╡ 931374e3-8c34-4d8e-b820-b066479d7c07
begin
	data1=[5,0,0,54,5,12,27,24,36,5,56,43,15,42,12,62,36,34,58,23]
	n_bins = range(minimum(data),maximum(data),10)
	fig1 = Figure()
	ax1 = Axis(fig1[1,1])
		hist!(ax1,data1,bins=n_bins,strokewidth = 1, strokecolor = :black,normalization= :probability)
	fig1
end

# ╔═╡ ed7533c9-f235-4cce-8d10-793e115f6fdc
md"""
Dacă variem numărul de cutii obținem alte forme:
"""

# ╔═╡ 4fd1a3dd-5a23-4854-a56e-49be07d731e9
begin
	data2=[5,0,0,54,5,12,27,24,36,5,56,43,15,42,12,62,36,34,58,23]
	n_bins2 = range(minimum(data2),maximum(data2),9)
	fig2 = Figure()
	ax2 = Axis(fig2[1,1],title="9 cutii")
		hist!(ax2,data2,bins=n_bins2,strokewidth = 1, strokecolor = :black,normalization= :probability)
	fig2
end

# ╔═╡ e8d5b769-79fc-4249-8d7e-d15157235a77
begin
	n_bins3 = range(minimum(data),maximum(data),11)
	fig3 = Figure()
	ax3 = Axis(fig3[1,1],title="11 cutii")
		hist!(ax3,data,bins=n_bins3,strokewidth = 1, strokecolor = :black,normalization= :probability)
	fig3
end

# ╔═╡ 3b4efbfb-6e1f-4b48-ba3d-aef6dc3ba97f
md"""
Histogramele ne ajută să vedem dacă datele sînt simetrice în jurul medianei sau se concetrează într-o parte. În acest caz spunem că datele sînt **asimetrice** (skewed). Dacă avem multe valori mai mici decît media spunem că valorile sînt concentrate la stînga, în celălalt caz că sînt concentrate la dreapta.
 

Adesea este util să aproximăm datele discrete precum o histogramă printr-o funcție continuă. Pentru foarte multe date, chiar o cutie îngustă va conține multe puncte. Histograma de mai jos reprezintă valorile pulsului unor pacienți internați la terapie intensivă. Funcția care aproximează histograma este funcția normală:
```math
f(x)=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}.
```
Această funcție are un singur punct de maxim în $0$. Ca să adaptăm funcția la situația noastră introducem doi parametri $\mu$ și
$\sigma$, numite **centru**, respectiv **dispersie**:
```math
f(x)=\frac{1}{\sigma \sqrt{2\pi}}e^{-\frac{\left(x-\mu\right)^2}{2\sigma^2}}.
```
"""

# ╔═╡ ba743f67-c441-407e-9501-b6165d101eba
begin
f(x;σ,μ) = 1/(σ*sqrt(2*π))*exp(-(x-μ)^2/(2*σ^2))
data3=[52,68,69,70,78,78,79,81,83,88,88,88,88,89,89,89,92,92,95,95,95,98,98,98,99,99,101,101,103,106,108,108,109,109,113,115,115,
119,128,139]
μ = mean(data3)
σ = std(data3)
n_bins4 = range(minimum(data3),maximum(data3),8)
# println(μ)
# println(σ)
histogram(data3,bins=n_bins4,leg=false,normalize =:pdf)
g(x)=f(x;σ,μ)
Plots.plot!(g,lw=5,color= :green)
end

# ╔═╡ 73c46e44-72bb-4a9e-9d77-ede852eaf80f
md"""
Funcția normală are un singur punct de maxim pentru $x=\mu$ cu valoarea $$\frac{1}{\sigma\sqrt{2\pi}}.$$ Graficul său este simetric față de dreapta
$x=\mu$. Ariile de sub grafic sînt astfel
 - 68.2% este între $\mu-\sigma$ și $\mu+\sigma$;
 - 95.5% este între $\mu-2\sigma$ și $\mu+2\sigma$;
 - 99.7% este între $\mu-3\sigma$ și $\mu+3\sigma$.

Mai jos avem graficele funcției normale pentru diferite valori ale lui $\sigma$ și $\mu$.
"""

# ╔═╡ 4cbd0e4a-b734-4c49-9fa8-5dd344bc2719
begin
	f(x,μ,σ) = 1/(σ*sqrt(2*π))*exp(-(x-μ)^2/(2*σ^2))
	xs = LinRange(-3,3,400)
	fig5=Figure(size=(1000,500))
	ax5=Axis(fig5[1,1])
	for (μ,σ) ∈ ((0,1.0),(0,0.4),(0,2.0),(-2,0.6))
		ys = [f(x,μ,σ) for x ∈ xs]
		Makie.lines!(ax5,xs,ys,linewidth=5,label="Funcția normală cu σ=$σ și μ=$μ")
	end
	axislegend(ax5)
	fig5
end

# ╔═╡ 6f95003a-23a7-493f-a279-11e42c6988b9
md"""
# Relații între variabile
## Relații între variabile categoriale
Dacă avem două variabile categoriale care caracterizează un set de date putem construi **tabela de contingență.** Pentru a studia efectele antiviralului Rimantadină au
fost luați 76 de copii împărțiți în două grupuri: unuia i s-a administrat  medicamentul, celorlalți un placebo. Acestei populații i se asociază două variabile categoriale: dacă au primit medicament sau placebo și dacă au fost sau nu bolnavi:

|                  | Infectați  | Neinfectați | Total pe linie |  
|------------------|------------|-------------|----------------|
| Placebo          | 20         | 21          | 41             |  
| Rimantadină      | 1          | 34          | 35             | 
| Total pe coloană | 21         | 55          | 76             |


Tabelul sugerează că faptul de a fi infectat se corelează cu faptul de a fi luat medicamentul. Putem folosi un desen cu bare pentru a vedea că în grupul care a primit rimatadină aproape nimeni nu s-a infectat, în vreme ce în grupul placebo aproximativ jumătate s-a infectat.
"""

# ╔═╡ 88999291-b9de-44bf-8d40-8eda791f8060
begin
	ctg = repeat(["Infectați","Neinfectați"],inner = 2)
nume = repeat(["Placebo","Rimatadină"],outer = 2)
placebo=[20,21]
rima = [1,34]
vals = [20,1,21,34]
groupedbar(nume, vals, group = ctg)
end

# ╔═╡ a577b98e-559f-46d4-8504-3bea2af59882
md"""
## Date numerice și categoriale

Dacă avem de comparat o variabilă numerică și una categorială, transformăm variabila numerică într-una ordinală ca mai sus și comparăm cele două variabile.

## Două variabile numerice

Pentru început este cel mai bine să desenăm un **scatter plot**. Anume luăm una dintre variabile pe axa $x$ și cealaltă pe $y$. Un studiu a fost efectuat în Marea Britanie pentru a vedea legătura între fumat și mortalitate. Un eșantion de bărbați a fost împărțit în 25 de categorii bazate pe ocupație și în fiecare grup s-a măsurat
consumul mediu de țigări și probabilitatea să moară de cancer de plămîni.
"""

# ╔═╡ f1cf2e4a-331f-494c-9fa0-84becaadbf7d
begin
x = [77,102,115,137,91,105,117,104,87,94,107,91,116,112,100,102,113,76,111,110,66,93,125,88,133]
y = [84,88,128,116,104,115,123,129,79,128,86,85,155,96,120,101,144,60,118,139,51,113,113,104,146]
fig4=Figure()
ax4 = Axis(fig4[1,1])
Makie.scatter!(ax4,x,y,color=:red)
fig4
end

# ╔═╡ 035ec8de-9760-4f4b-9635-9594ed73a793
md"""
Se vede că persoanele care fumează mai mult au o probabilitate mai mare de a muri de cancer pulmonar. Adesea însă corelația nu este foarte clară doar din desen. În cazul nostru pare să existe o relație liniară între cele două variabile. Ca să vedem în mod riguros acest lucru trebuie să desenăm o dreaptă de care punctele noastre 
se apropie cel mai mult. 

Metoda cea mai folosită este **metoda celor mai mici pătrate**, adică trebuie găsită dreapta pentru care suma distanțelor față de puncte este minimă. Începem prin a calcula distanța verticală de la puncte la dreapta $y=ax+b$ ($a$ și $b$ trebuie găsite)
```math
\Delta y_i = y_i-(ax_i+b)
```
le ridicăm la pătrat și le adunăm. Avem deci $S=[y_1-y(x_1)]^2+\dots+[y_n-y(x_n)]^2.$ $S$ este o funcție de două variabile $S=S(a,b)$ pe care trebuie să o minimizăm.
Vom obține 
```math
\begin{eqnarray*}
a=\frac{\overline{xy}-\overline{x}\overline{y}}{\overline{x^2}-\overline{x}^2} \\
b=\overline{y}-a\overline{x},
\end{eqnarray*}
```
unde cu $\overline{}$ am notat media. Revenim la exemplul nostru
"""

# ╔═╡ 23efaf54-7395-4d9a-b4ea-2c11efac17be
begin
	x1 = [77,102,115,137,91,105,117,104,87,94,107,91,116,112,100,102,113,76,111,110,66,93,125,88,133]
xsq = x.^2
y1 = [84,88,128,116,104,115,123,129,79,128,86,85,155,96,120,101,144,60,118,139,51,113,113,104,146]
ysq = y1.^2
xpy = x1.*y
mx = mean(x)
my = mean(y)
mxy = mean(xpy)
mx2 = mean(xsq)
my2 = mean(ysq)
a = (mxy-mx*my)/(mx2-mx^2)
b = my-a*mx
y0(x) = a*x+b
Plots.scatter(x,y,color=:red,leg=false,xlims=(minimum(x)-1,maximum(x)+1),ylims=(minimum(y)-1,maximum(y)+1))
Plots.plot!(y0,xlims=(minimum(x)-1,maximum(x)+1),ylims=(minimum(y)-1,maximum(y)+1),lw=5)
end

# ╔═╡ 5a02d5c3-c9e9-48ba-9e24-9d7886675a53
md"""
# Recapitulare combinatorică
O **permutare** a unei mulțimi de obiecte este o aranjare ordonată a elementelor. Dacă avem $n$ obiecte atunci numărul de permutări
se poate calcula astfel: avem $n$ posibiltăți pentru alegerea obiectului de pe locul întîi, apoi $n-1$ posibilități pentru a 
doua poziție, $n-2$ pentru a treia etc. Deci  numărul de permutări de $n$ obiecte este $n\cdot(n-1)\cdot(n-2)\cdot\dots1=n!.$

Dacă luăm cîte $k$ elemente din cele $n$, numărul permutărilor posibile este 
```math
A_n^k=\frac{n!}{(n-k)!}.
```

Spre exemplu dacă într-un vas se află $4$ bile colorate (de culori diferite). Dacă luăm cîte două bile o dată, atunci 
avem $$\frac{4!}{2!}=\frac{24}{2}=12$$ 
posibilități.

O **combinație** este submulțime *neordonată* a unei mulțimi. Deci într-o combinație nu contează ordinea elementelor. Numărul
de combinații de $n$ elemente luate cîte $k$ este
```math
C_n^k=\frac{A_n^k}{k!}=\frac{n!}{k!(n-k)!}.
```

Dacă într-un vas sînt $6$ bile colorate și luăm cîte $4$, atunci avem 
```math
C_6^4=\frac{6!}{2!4!}=\frac{5\cdot6}{2}=15.
```

Un alt exemplu este dat de comportamentul HIV. După infectarea unei persoane cu HIV acesta suferă numeroase mutații rezultînd
numeroase genotipuri. La transmitere însă, doar cîteva genotipuri se transmit. Să presupunem că o persoană infectată posedă
150 de genotipuri, dar transmite doar 10 genotipuri. Cîte posibilități pentru transmiterea genotipurilor sînt?

Răspunsul este 
```math
C_{150}^{10}=\frac{150!}{10!140!}=1.169.554.298.222.310
```
"""

# ╔═╡ 4f68b47b-4c2b-42da-a53b-7934639bcd17
begin
N1=factorial(big(150))
println(N1,"\n")
N2=factorial(big(140))
println(N2,"\n")
N3=factorial(10)
println(N3)
N=N1/(N2*N3)
println(N)
binomial(150,10)
end

# ╔═╡ b919f79b-fc64-4356-9403-7b8333a400fe
md"""
# Elemente de teoria probabilităților

## Probabilități

Teoria Probabilităților măsoară incertitudinea. Avem nevoie de unele noțiuni pentru a aplica riguros această măsură. O noțiune 
centrală este aceea de **experiment**: este o procedură care poate fi repetată de un număr indefinit de ori în condiții 
identice. Fiecare aplicare a unui experiment se numește **încercare**. Mulțimea tuturor rezultatelor posibile se numește
**spațiul de selecție**.

Spre exemplu pentru aruncarea unei monezi, spațiu de selecție este $\Omega=\left\{H,T\right\}$. Pentru aruncarea unui zar 
$\Omega=\left\{1,2,3,4,5,6\right\}$. Orice submulțime a spațiului de selecție se numește **eveniment**. Dacă evenimentul
constă dintr-un singur rezultat acesta se numește simplu, iar dacă nu se numește compus.

Pentru un spațiu de selecție $\Omega$ și o mulțime de evenimente o **probabilitate** este o funcție $P$ care asociază 
oricărui eveniment un număr $P(E)$ care reprezintă șansa ca evenimentul $E$ să aibă loc pentru un experiment. Intuitiv este
numărul de apariții ale lui $E$ pentru un număr suficient de mare de experimente.

Dacă de exemplu experimentul constă în aruncarea unui zar și evenimentul este să iasă $3$, atunci $P(E)=\frac{1}{6}.$ Dacă 
evenimentul este să iasă o față impară adică $1,3$ sau $5$ atunci probabilitatea să iasă oricare din fețe este egală cu 
$\frac{1}{6}$ și deci pentru evenimentul nostru $P(E)=\frac{1}{2}.$ Orice funcție de probabilitate trebuie să satisfacă două
proprietăți fundamentale:
 1. $0 \leq P(E) \leq 1;$
 2. $P(\Omega)=1.$

Două evenimente $E$ și $F$ se numesc exclusive dacă $E \cap F = \emptyset$. Pentru astfel de evenimente probabilitatea trebuie
să satisfacă: 

   3. $P(E \cup F)=P(E)+P(F).$

Cele trei condiții 1,2,3 se numesc *axiomele probabilităților*. O primă consecință importantă a axiomelor este un fapt care
este extrem de util în multe situații: 
```math
P(\Omega \setminus E)=1-P(E).
```

Dacă avem două evenimente $E$ și $F$ nu neapărat exclusive atunci probabilitatea reuniunii lor este
  
  4. $P(E \cup F)=P(E)+P(F)-P(E\cap F).$

Dacă toate evenimentele din $\Omega$ au aceeași probabilitate, atunci avem 

```math
P(E)=\frac{n(E)}{n(\Omega)},
```
unde $n(\cdot)$ reprezintă numărul de elemente.

## Probabilități condiționate și independență
Unul din experimentele lui Gregor Mendel se referea la culoarea florilor de mazăre. Aceasta este determinată de o singură genă cu două alele $C$ și $c$. Fiecare plantă are două copii ale genei, deci trei genotipuri posibile $CC$, $Cc$ sau $cc$. Primele două determină flori roz, iar ultimul flori albe. Presupunem că încrucișăm doi părinți cu genotipul $Cc$ și obținem o plantă cu flori roz. Care este probabilitatea ca aceasta să aibă genotipul $CC$? Cum florile sînt roz genotipurile posibile vor fi
$CC$, $Cc$ sau $cC$. Deci probabilitatea este $\frac{1}{3}$. Acest tip de probabilitate calculat avînd condiții suplimentare se numește **probabilitate condționată**.

Noțiune de probabilitate condiționată este fundamentală în teoria probabilităților. Aceasta se calculează impunînd niște 
condiții suplimentare în spațiul de selecție. De exemplu dorim să calculăm probabilitatea să iasă 1 dacă rezultatul aruncării
este impar. Evident în acest caz spațiul de selecție este format doar din $1,3$ și $5$. Deci probabilitatea căutată este 
$\frac{1}{3}$.

Riguros să presupunem că $E$ și $F$ sînt evenimente cu $P(F) > 0.$ Atunci 
```math
P(E|F)=\frac{P(E\cap F)}{P(F)}.
```
$P(E|F)$ se numește **probabilitatea lui E, condiționată de F**. 

Boala lui Huntington este o afecțiune genetică produsă de o mutație într-o singură genă. Oamenii pot avea $0$, $1$ sau $2$ copii
ale mutației. Cazurile corespunzătoare se numesc homozigotic normal, heterozigot, homozigotic mutant. Indivizii heterozigoți
și cei homozigoți mutanți fac boala. Care este probabilitatea ca un copil să fie homozigot, dacă știm că are boala?

Fie $D$ evenimentul „bolnav” și $H$ evenimentul „homozigot”. Vedem că $P(D)=\frac{3}{4}$ și că $P(H \cap D)=\frac{1}{4}$, deci
```math
P(H|D)=\frac{P(H \cap D)}{P(D)}=\frac{\frac{1}{4}}{\frac{3}{4}}=\frac{1}{3}.
```

Formula pentru probabilitatea condiționată se poate rescrie $P(E\cap F)=P(E|F)P(F)$. Dacă fenomenul $F$ nu influențează $E$, atunci $P(E|F)=P(E)$. Spunem că $E$ și $F$ sînt **independente** dacă $$P(E \cap F)=P(E)P(F).$$ 

Să arătăm că evenimentele  $E:$„să iasă 10” și $F:$ „să fie caro” sînt indepente. Cum avem $52$ de cărți si patru dintre ele 
sînt 10 avem 
```math
P(E)=\frac{4}{52}=\frac{1}{13}.
```
Pentru că sînt cîte $13$ cărți de aceeași culoare avem 
```math
P(F)=\frac{13}{52}=\frac{1}{4}.
```
Probabilitatea să tragem 10 de caro este
```math
P(E \cap F)=\frac{1}{52}=\frac{1}{4}\cdot\frac{1}{13}=P(E)P(F).
```
Deci $E$ și $F$ sînt independente.

Putem să folosim probabilități condiționate ca să calculăm probabiltățile necondiționate, dacă reușim să partiționăm
spațiul de selecție în părți în care calculele sînt mai ușoare. Concret să presupunem că 
```math
\Omega = F_1 \cup F_2 \cup \dots \cup F_k
```
cu toate $F_j$ disjuncte două cîte două. Atunci pentru orice eveniment $E$
```math
P(E)=P(E|F_1)P(F_1)+P(E|F_2)P(F_2)+\dots+P(E|F_k)P(F_k).
```
Această relație se numește **legea probabilității totale.**
Să presupunem că un bărbat cu boala Huntington și o femeie fără boală decid să aibă un copil. Care este probabilitatea ca acesta
să fie bolnav?

Notăm cu $A$ mutația care dă boala și cu $a$ mutația care nu o dă. Atunci genotipul mamei este $aa$, iar cel al tatălui poate fi
$AA$ sau $Aa$. Notăm cu $F_1$ primul eveniment și cu $F_2$ cel de-al doilea. Presupunem că $P(F_1)=0.05$ și $P(F_2)=0.95.$ 
Notăm cu $E$ evenimentul „copilul are boala”. În primul caz copilul va primi o genă $A$ de la tată deci cu siguranță se va 
îmbolnăvi prin urmare $P(E|F_1)=1$. În cazul al doilea copilul primește de la tată fie gena $A$ fie $a$. Deci probabilitatea
$P(E|F_2)=\frac{1}{2}.$ Folosind principiul aditivității avem
```math
P(E)=P(E|F_1)P(F_1)+P(E|F_2)P(F_2)=1\cdot 0.05+\frac{1}{2}0.95=0.525.
```
"""

# ╔═╡ Cell order:
# ╟─9e1db0e6-d997-11f0-348a-d9f585410f74
# ╟─9f0d0ed7-1ed4-469d-8bcd-80f89965afb7
# ╟─bc79a3dd-ccc3-43be-b90d-048cf6059b4f
# ╟─f7a6e388-7d3e-4eb6-8788-bb13dad71898
# ╟─13dfb805-fdb6-4e4a-93c1-154303f4b331
# ╟─51f4e264-f290-46f1-838d-db3a4dd74ba6
# ╟─faefd827-a505-4f0d-b338-d443860a0c17
# ╟─dfb42990-8e1a-4b3e-ad1d-11054d7b0c1d
# ╟─dbebc2ae-a819-447c-9ce8-a0d63d5eb673
# ╠═f163a685-b7f1-49b2-bf25-6edb98904eb4
# ╟─7001202f-936f-4982-beb5-31627d20ba68
# ╠═931374e3-8c34-4d8e-b820-b066479d7c07
# ╠═ed7533c9-f235-4cce-8d10-793e115f6fdc
# ╠═4fd1a3dd-5a23-4854-a56e-49be07d731e9
# ╠═e8d5b769-79fc-4249-8d7e-d15157235a77
# ╟─3b4efbfb-6e1f-4b48-ba3d-aef6dc3ba97f
# ╟─ba743f67-c441-407e-9501-b6165d101eba
# ╟─73c46e44-72bb-4a9e-9d77-ede852eaf80f
# ╟─4cbd0e4a-b734-4c49-9fa8-5dd344bc2719
# ╠═6f95003a-23a7-493f-a279-11e42c6988b9
# ╟─88999291-b9de-44bf-8d40-8eda791f8060
# ╟─a577b98e-559f-46d4-8504-3bea2af59882
# ╟─f1cf2e4a-331f-494c-9fa0-84becaadbf7d
# ╟─035ec8de-9760-4f4b-9635-9594ed73a793
# ╟─23efaf54-7395-4d9a-b4ea-2c11efac17be
# ╟─5a02d5c3-c9e9-48ba-9e24-9d7886675a53
# ╟─4f68b47b-4c2b-42da-a53b-7934639bcd17
# ╟─b919f79b-fc64-4356-9403-7b8333a400fe
