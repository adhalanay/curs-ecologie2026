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

[Notebook-uri](https://codeberg.org/adhalanay/curs_ecologie2025) 

# Cum se va face notarea
   - 40 de puncte pentru activitate (prezențe la curs și laborator, răspunsuri, teme etc.);
   - 30 de puncte un proiect (poate fi făcut în echipe de maxim 2 persoane);
   - 60 de puncte examenul final.

Pentru cei care doresc va fi un examen parțial în săptămîna a 8-a.

Condiții de trecere: cel puțin 50 de puncte în total și cel puțin 30 de puncte la examen.
"""

# ╔═╡ 11efb0cf-73b4-4d19-b0f7-612a5abec2d3
md"""
# Sisteme cu feedback

Modelarea matematică a unor procese biologice a început să fie folosită începînd cu mijlocul secolului XIX. Printre numeroasele succese ale se pot menționa premiile Nobel pentru

    - Ronald Ross în 1902 pentru înțelegerea procesului de transmisie al malariei;
    - Alan Lloyd Hodgkin și Andrew Fielding Huxley în 1963 pentru înțelegerea transmisiei impulsului nervos;
    - Alan Cormack și Geoffrey Housfield în 1979 pentru dezvoltarea metodologiei din spatele tehnologiei CT.

În zilele noastre matematica are aplicații în numeroase domenii ale biologiei: este indispensabilă în înțelegerea unor aspecte fundamentale precum funcționarea sistemului imunitar și apariție bolilor autoimune sau pentru înțelegerea mecanismelor rezistenței la anumite medicamente. De asemenea designul modern al medicamentelor sau al unor dispozitive medicale folosește în mod fundamental metode matematice. Procesele de mișcare a organismelor vii luate atît individual cît și colectiv, spre exemplu dinamica populațiilor, de asemenea pentru înțelegerea efectelor modificării habitatelor și a recoltării.

Putem considera ca unul din primele studii științifice de ecologie este reprezentat de studiul fluctuațiilor a două populații de animale: Rîsul canadian (Lynx canadensis):

"""

# ╔═╡ 9f282b4f-da9c-49b1-94ef-03cbe39ca1f9
PlutoUI.LocalResource("Canadalynx.jpg",:width => 300)

# ╔═╡ 4d08da5d-e63f-49d9-a25a-05c9ff0c2cd4
md"""
care se hrănește aproape exclusiv iepurele de zăpadă (Lepus americanus):

"""

# ╔═╡ 9968fa47-4085-4417-8100-7c8385175aca
PlutoUI.LocalResource("SnowshoeHare.jpg",:width=>300)

# ╔═╡ ccd56781-7fbe-48ed-9617-6c88829384c0
md"""
Numărul de piei de rîși și de iepuri este reprezentat de graficul:
"""

# ╔═╡ d42b9ecf-086a-47fe-b2e8-742a3b0ed214
PlutoUI.LocalResource("fourrures.jpg")

# ╔═╡ 474aa1f8-42b2-44d7-a7ef-441d0ddffa36
md"""
sau dacă ne limităm la o perioadă mai scurtă:
"""

# ╔═╡ a6e282d9-dd41-440d-a8c6-bc068f23720a
PlutoUI.LocalResource("rîși.png")

# ╔═╡ 391dd0e7-66e6-42cb-a967-0345457f85f7
md"""
(reprodus după A.Garfinkel et. al. 2017) 

Numărul de rîși este de fapt numărul de piei de rîs colectate de vînătorii companiei Hudson Bay Company în aproape 100 de ani (numărul de iepuri a fost doar estimat). 

Se observă că populațiile oscilează, dar aceste oscilații nu sînt la întîmplare ci au aproximativ o perioadă de 10 ani. De asemenea se observă că populația de predători crește sau scade cu o oarecare întîrziere față de cea a prăzii.

Dacă vrem să analizăm care este cauza acestor oscilații trebuie să **modelăm** interacțiunea între cele două specii. Chiar dacă nu este foarte riguros, observăm că populația de pradă influențează în mod pozitiv populația de prădători, iar dimpotrivă populația de prădători influențează în mod negativ populația prăzii. Dacă știm numărul de rîși și numărul de iepuri la un moment dat și vrem să estimăm cum vor evolua cele două populații, pentru început pare clar că numărul de iepuri va scădea, iar cel de rîși va crește. Dar pe măsură ce numărul de rîși crește, iar cel de iepuri scade, populația prădătorilor va deveni nesustenabilă și va începe să scadă. Este greu de spus cu precizie ce se întîmplă mai departe. De aceea avem nevoie de un model matematic care să modeleze evoluția populației. Următorul aspect este deosebit de important

**Fiecare populație determină valoarea celeilalte populații**.

Spunem că avem feedback pozitiv, dacă o valoare pozitivă a unei variabile determină creșterea acesteia. Spre exemplu

  - dacă o persoană are mulți bani poate să-i investească pentru a obține mai mulți;
  - animalele fac pui ceea ce duce la creșterea populației deci la creșterea numărului de pui, prin urmare populația     va crește atîta timp cît resursele o permit;
  - creșterea nivelului de $CO_2$ duce la creșterea temperaturii, ceea ce duce la creșterea metabolismului microorganismelor din sol care descompun materia organică și deci produc mai mult $CO_2$.
  
Avem feedback negativ dacă o valoare pozitivă a unei variabile duce la scăderea acesteia, iar o valoare negativă duce la creșterea acesteia. Un exemplu de astfel de sistem este cel al unui sistem de aer condiționat: cînd temperatura crește, sistemul intră în funcțiune și o scade. Putem presupune că sistemul poate și încălzi aerul și atunci cînd temperatura scade sub un anumit nivel începe să încălzească. Mai precis dacă notăm cu $T_0$ temperatura setată, și cu $C$ temperatura curentă definim variabila temperatură prin $T=C-T_0$. Deci dacă $T > 0$ sistemul răcește, iar dacă $T<0$ sistemul încălzește.

Un alt exemplu clasic de sistem cu feedback negativ este dat de concentrația de glucoză/insulină în sînge. Cînd concentrația de glucoză crește, atunci crește și cea de insulină ceea ce duce la scăderea glicemiei. Pentru o persoană care are o pompă care eliberează in mod continu insulină concentrația de insulină/glucoză variază astfel:
"""

# ╔═╡ 9dd94e6d-04e5-4669-adff-414b25cb1dd3
PlutoUI.LocalResource("insulina.png")

# ╔═╡ 291381ae-569e-4335-9110-b3af9bee2c9f
md"""

Sturis, J., Polonsky, K. S., Blackman, J. D., Knudsen, C., Mosekilde, E., & Van Cauter E.
(1991a). Aspects of oscillatory insulin secretion. In E. Mosekilde & L. Mosekilde (Eds.),
Complexity, Chaos, and Biological Evolution (vol. 270, pp. 75–93). New York: Plenum Press.
Sturis, J., Polonsky, K. S., Mosekilde, E., & Van Cauter, E. (1991b). Computer model for mech-
anisms underlying ultradian oscillations of insulin and glucose. American Journal of Physiology-
Endocrinology and Metabolism, 260 (5), E801–E809. 

foarte asemănător cu sistemul rîși/iepuri.

În epidemiologie contactele între persoanele infectate și cele susceptibile duce la creșterea numărului de infectați și scăderea numărului de susceptibili, deci numărul de infecții va scădea. Modelele epidemiologice sînt folosite pentru a elabora strategiile de intervenție în cazul unei epidemii. 

Sistemele cu feedback au adesea comportări ne intuitive. Să presupunem că vrem să reducem numărul de prădători prin scoaterea lor din mediu. Atunci cantitatea de pradă va crește ducînd la creșterea numărului de prădători peste valoarea inițială. Acesta este fenomenul de **rebound**.

Rezultatul unei intervenții într-un astfel de sistem depinde de faza ciclului cînd are loc intervenția. 

Sisteme aparent simple au comportamente foarte ne intuitive, deci formularea și studiul modelelor matematice este indispensabilă.

Primul concept de care vom avea nevoie este acela de funcție.

# Funcții

Un concept fundamental de care avem nevoie pentru a descrie un fenomen este cel de **funcție**. Am văzut în graficele prezentate că fiecărui moment îi corespunde exact o valoare. Mai precis o funcție este o relație între o mulțime de date de intrare și una de date de ieșire astfel încît fiecărui input îi este asociat **exact** un output, nu mai multe și nu nici unul. 

Funcțiile pot fi definite cu ajutorul unui tabel cu două coloane de lungime egală, dar cel mai adesea vor fi definite prin formule. Spre exemplu funcția $f(x)=x^2$ asociază fiecărui număr $x$ valoarea $x^2$. Putem scrie formula $y=x^2$. 

Două aspecte sînt extrem de importante:
  - Nu orice formulă definește o funcție: spre exmplu formula $y^2=x^2$ nu definește o funcție deoarece pentru $x=2$ avem două valori pentru $y$, anume $y=2$ și $y=-2;$
  - Nu orice serie temporală (grafic) poate fi descrisă cu o formulă explicită. De fapt nici unul din graficele prezentate nu este graficul unei funcții explicit formulate.
  
Mulțimea valorilor de intrare ale unei funcții se numește **domeniul** acesteia. De obicei domeniul funcțiilor va fi mulțimea numerelor reale $\mathbb{R}$, sau o submulțime a acesteia. Putem alege ce domeniu dorim pentru funcția noastră, atîta timp cît ea are sens. Spre exemplu pentru funcția $f(x)=\frac{1}{x}$ domeniul nu trebuie să-l conțină pe $0$, dar poate fi altfel arbitrar. 

Mulțimea valorilor de ieșire se numește **codomeniul** acesteia. Spre exemplu codomeniul funcției $f(x)=x^2$ este mulțimea $\mathbb{R}_+$ a numerelor reale pozitive. Pentru o funcție cu domeniu $X$ și codomeniu $Y$ notăm $f:X \to Y$.

Un exemplu interesant provine din biologia moleculară. După cum se știe ADN-ul este compus dintr-un șir de patru baze $A,C,G,T.$ Cînd ADN-ul este transpus în ARN mesager baza $T$ este înlocuită de o altă bază $U$. Avem deci o funcție 
```math
transcriere: \left\{A,C,G,T\right\} \to \left\{A,C,G,U\right\}.
```

De fapt transcrierea ia o tripletă de baze ADN (codonul) și o duce într-un codon ARN. Un codon ARN determină ca un anumit amino-acid să fie adăugat unei proteine. Cum fiecare codon determină un singur amino-acid, avem o nouă funcție, numită translație 
```math
translație : \text{codoni ARN} \to \text{amino-acizi}.
```
Avem un proces numit expresie genetică 
```math
expr-gen: \text{codoni ADN} \to \text{amino-acizi. }
```
Funcția se obține aplicînd succesiv cele două funcții. Acest proces se numește compunere de funcții. 

Un alt exemplu: dacă $f(x)=2x^2+1$ și $g(x)=\sqrt{x}$, atunci $(g\circ f)(x)=\sqrt{2x^2+1}$ și $(f\circ g)(x)=
2x+1.$

**Exercițiu:** Presupunem că există viață pe Marte similară cu o structură genetică foarte asemănătoare cu a noastră, dar codonii ARN se comportă diferit. Spre exemplu codonul $AUC$ produce în 60% din cazuri izoleucină și în 40% treonină. Mai este expresia genetică o funcție?

O funcție poate fi reprezentată în patru forme:
  - verbal (prin cuvinte);
  - numeric (printr-un tabel de valori);
  - visual (printr-un grafic);
  - algebric (printr-o formulă).

# Reprezentarea grafică
Graficul obișnuit al unei funcții folosește un sistem de două axe perpendiculare:
 - axa orizontală $Ox$ reprezintă valorile variabilei independente $x$;
 - axa verticală $Oy$ reprezintă valorile funcției $y=f(x).$

Spre exemplu dacă $f:\mathbb{R} \to \mathbb{R}$, $f(x)=x^5$ atunci graficul său este

"""

# ╔═╡ e4d39336-7ddd-4619-9288-6169aadd6a0e
begin
	f(x) = x^5
	xs = range(0,5,length=100)
	ys = [f(x) for x ∈ xs]
	fig = Figure()
	ax = Axis(fig[1,1])
	lines!(xs,ys)
	fig
end

# ╔═╡ 839ddca1-d3ba-4d99-b15f-a795b040c161
md"""
Uneori valorile lui $y$ sînt foarte mari față de valorile lui $x$. Atunci este mai folositor să folosim **scara logaritmică**. Astfel de scări apar atunci cînd considerăm pH-ul, scara Richter, scara decibelică. În aceste cazuri valorile corespund exponenților lui $10$, adică $1$ corespunde lui $10$, $2$ lui $100=10^2$, 3 lui $1000=10^3$ etc. 

Adesea în biologie se folosesc graficele **semi-logaritmice**, anume reprezentăm perechile $(x,\log y)$, unde logaritmul este în baza 10. Pentru funcția noastră, acest grafic este:
"""

# ╔═╡ 019adb7d-8418-4149-a81b-ec99341efe50
begin
	fig1 = Figure()
	ax1 = Axis(fig1[1,1],yscale=log10)
	lines!(xs,ys)
	fig1
end

# ╔═╡ 5477506e-4926-41d0-a31c-7bf486040402
md"""
Dacă într-un astfel de grafic obținem o dreaptă atunci spunem că datele noastre satisfac o lege exponențială.
 
Dacă luăm scările logaritmice atît pentru $x$ cît și pentru $y$ atunci obținem o reprezentare **logaritmică**.
"""

# ╔═╡ 803b1d6b-e26b-44eb-917d-d0265c228b14
begin
	fig2 = Figure()
	ax2 = Axis(fig2[1,1],xscale=log10,yscale=log10)
	lines!(xs,ys)
	fig2
end

# ╔═╡ 4f5a7238-8855-44b1-8afa-b774958caaf0
md"""
# Spațiul stărilor

Variabilele de stare sînt instrumentele de bază pentru a descrie cantitativ sistemele biologice la un moment dat. Spre exemplu starea unei populații este mărimea ei. Este extrem de important să alegem ce variabile să studiem în funcție de problema care ne interesează. În cazul populației putem fi interesați de reportul între sexe, distribuția
pe clase de vîrstă, distribuția teritorială etc. Cea mai importantă și dificilă parte a creerii unui model este alegerea unor variabile de stare relevante precum și a unor unități de măsură potrivite. Această alegere depinde atît de structura sistemului (de exemplu animalele sînt distribuite omogen pe teritoriu sau există subpopulații distincte) cît și de ceea ce urmărim cu modelul nostru.

Cum orice variabilă de stare poate lua o singură valoare la un moment dat:**Variabilele de stare sînt funcții de timp**. Scopul modelării matematice este să înțelegem evoluția în timp a sistemului, deci este important să considerăm toate valorile posibile pentru variabilele de stare. Mulțimea tuturor acestor valori posibile se numește
**spațiul stărilor.** Vom presupune că spațiul stărilor este o muțime de numere, și de asemenea vom presupune că pot lua orice valori între anumite limite (ipoteza de continuitate).

Dacă avem o singură variabilă de stare atunci spațiul stărilor va fi mulțimea numerelor reale sau un interval.

# Sisteme de dimensiune 1

În funcție de mărimea pe care o considerăm spațiul stărilor poate fi:
  - ``[-\infty,\infty]`` pentru voltaj;
  - ``[0,\infty]`` pentru mărimea populației, temperatură (în $^\circ K$);
  - ``[0,1]`` pentru procentul unui grup într-o populație.

**Exercițiu:** Care va fi spațiul stărilor pentru temperatura măsurată în 
``^\circ C``?

Cum orice variabilă de stare este funcție de timp, putem să o desenăm ca un grafic 
cu timpul pe axa orizontală și cu valoarea pe axa verticală. Acest tip de reprezentare se numește **serie temporală**. Astfel de serii sînt numărul de rîși și cel de iepuri.

# Sisteme de dimensiune mai mare

Dacă avem mai multe variabile de stare, de exemplu pentru rîși și iepuri, atunci spațiul stărilor va fi format din perechi de numere $(R,I).$ Chiar dacă spațiul stărilor nu mai este format din numere reale putem să redefinim două operații
   - adunarea lor: $(R_1,I_1)+(R_2,I_2)=(R_1+R_2,I_1+I_2)$ (adunarea pe componente);
   - înmulțirea cu un număr: $a(R,I)=(aR,aI)$ (înmulțirea cu scalari).
 
Dacă spațiul stărilor este format din perechi de numere atunci putem desena cîte o axă pentru fiecare componentă:
"""

# ╔═╡ 7d5f4320-933e-4358-8322-dcf613a8e702
begin
	fig3 = Figure()
	ax3 = Axis(fig3[1,1])
	x = [2,1,3,5]
    y = [3,2,4,1]
    scatter!([(2,3),(1,2),(3,4),(5,1)])
	fig3
end

# ╔═╡ 920b2efd-f83f-4300-b34f-b9695cf1079b
md"""
Vom numi în continuare o pereche de numere un **2-vector**, iar fiecare număr care îl constituie va fi numit o **componentă**. Un vector poate fi reprezentat ca o pereche ca mai sus sau ca o săgeată:
"""

# ╔═╡ 5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
begin
	fig4 = Figure()
	ax4 = Axis(fig4[1,1])
	arrows!([0],[0],[1.985],[2.985])
	scatter!(1.985,2.985,markersize=7)
	fig4
end

# ╔═╡ 5b1a6088-d669-4892-9ab4-1068dc68ec3f
md"""
Această reprezentare ne permite o interpretare geometrică a înmulțirii cu scalari: dacă modulul scalarului este mai mic decît $1$, atunci vectorul se scurtează, dacă modulul este mai mare decît $1$, acesta se lungește.

Adesea sîntem nevoiți să considerăm mai multe variabile. Atunci pentru a descrie spațiul fazelor vom avea nevoie de mai multe axe, cîte una pentru fiecare variabilă. Numărul de axe necesar se numește **dimensiunea** spațiului stărilor. Evident că nu putem vizualiza astfel de vectori, dar va trebui să lucrăm cu vectori cu orice număr de
componente. Ca și în cazul $2-$vectorilor putem aduna doi vectori:
  - Dacă $\mathbf{a}=(a_1,\dots,a_n)$ și $\mathbf{b}=(b_1,\dots,b_n)$ atunci 
$\mathbf{a}+\mathbf{b}=(a_1+b_1,\dots,a_n+b_n).$ 

**Ca să putem aduna doi vectori, trebuie ca aceștia să aibă aceeași dimensiune.**

  - Dacă $\mathbf{a}=(a_1\dots,a_n)$, atunci $\alpha\cdot\mathbf{a}=(\alpha a_1,\dots,\alpha a_n).$

Mulțimea $n-$vectorilor se notează cu $\mathbb{R}^n$. 


"""

# ╔═╡ Cell order:
# ╠═3c523909-557c-4473-9736-4e271c3e958a
# ╟─a8e79b63-4624-4f3c-a985-5bf6a1f51d43
# ╟─9d43ec2a-5d67-4219-829e-fc747d8d4057
# ╟─215a9db0-9667-11f0-3e8b-db501bdbd2d5
# ╟─11efb0cf-73b4-4d19-b0f7-612a5abec2d3
# ╟─9f282b4f-da9c-49b1-94ef-03cbe39ca1f9
# ╟─4d08da5d-e63f-49d9-a25a-05c9ff0c2cd4
# ╟─9968fa47-4085-4417-8100-7c8385175aca
# ╟─ccd56781-7fbe-48ed-9617-6c88829384c0
# ╟─d42b9ecf-086a-47fe-b2e8-742a3b0ed214
# ╟─474aa1f8-42b2-44d7-a7ef-441d0ddffa36
# ╟─a6e282d9-dd41-440d-a8c6-bc068f23720a
# ╟─391dd0e7-66e6-42cb-a967-0345457f85f7
# ╟─9dd94e6d-04e5-4669-adff-414b25cb1dd3
# ╟─291381ae-569e-4335-9110-b3af9bee2c9f
# ╟─e4d39336-7ddd-4619-9288-6169aadd6a0e
# ╟─839ddca1-d3ba-4d99-b15f-a795b040c161
# ╟─019adb7d-8418-4149-a81b-ec99341efe50
# ╟─5477506e-4926-41d0-a31c-7bf486040402
# ╟─803b1d6b-e26b-44eb-917d-d0265c228b14
# ╟─4f5a7238-8855-44b1-8afa-b774958caaf0
# ╟─7d5f4320-933e-4358-8322-dcf613a8e702
# ╟─920b2efd-f83f-4300-b34f-b9695cf1079b
# ╟─5db5c824-059a-49b1-af9f-eb8f5ed5dbf8
# ╟─5b1a6088-d669-4892-9ab4-1068dc68ec3f
