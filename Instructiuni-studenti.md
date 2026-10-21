# Lucrul pe serverul JupyterHub al cursului

## Unde se păstrează ce lucrați

Tot ce se află în directorul home (`/home/jovyan`) se păstrează de la o sesiune la alta:

- `~/course` — notebook-urile cursului și modificările voastre;
- orice fișier sau director creat de voi în home.

Ce salvați în afara directorului home se pierde la oprirea serverului.

La fiecare pornire, materialele noi ale cursului sînt copiate în `~/course`. Fișierele existente
**nu sînt suprascrise**, deci modificările voastre rămîn. Pentru a primi versiunea nouă a unui
fișier, redenumiți-l sau ștergeți-l și reporniți serverul
(*File → Hub Control Panel → Stop My Server → Start My Server*).

## Pachetele cursului

Pachetele folosite la curs (CairoMakie, DifferentialEquations, PlutoUI etc.) sînt deja instalate
și precompilate; nu trebuie instalat nimic.

### În Pluto

Notebook-urile cursului activează singure mediul cursului. Pentru un notebook nou, începeți cu
celula următoare:

```julia
begin
    using Pkg
    Pkg.activate(joinpath(homedir(), "course"))
    using CairoMakie, PlutoUI, DifferentialEquations
end
```

Fără `Pkg.activate`, Pluto creează un mediu separat pentru notebook și descarcă și compilează
din nou pachetele la fiecare sesiune, ceea ce poate dura multe minute.

**Salvați imediat notebook-ul nou în `~/course`.** Pluto salvează notebook-urile noi implicit
într-un director din afara home-ului, care se pierde la oprirea serverului. În cîmpul cu numele
fișierului, din partea de sus a notebook-ului, scrieți calea completă, de exemplu
`/home/jovyan/course/tema1.jl`, și apăsați *Enter*.

### În Jupyter

Alegeți kernel-ul **Ecologie 2026**; `using CairoMakie` etc. funcționează imediat.

## Pachete suplimentare

Pachetele instalate de voi **nu se păstrează** între sesiuni: la oprirea serverului fișierele
lor se șterg. Se păstrează doar lista lor, dacă le instalați într-un mediu propriu din directorul home:

```julia
using Pkg
Pkg.activate(joinpath(homedir(), "mywork"))
Pkg.add("NumePachet")
```

În fiecare sesiune nouă activați același mediu și reinstalați pachetele din listă
(durează cîteva minute, în funcție de pachete):

```julia
using Pkg
Pkg.activate(joinpath(homedir(), "mywork"))
Pkg.instantiate()
```

- Nu folosiți `Pkg.add` direct din kernel-ul Jupyter, fără `Pkg.activate`: mediul acelui kernel
  nu se păstrează, iar după repornire nici pachetul, nici lista nu mai există.
- Nu adăugați pachete în mediul din `~/course`: se pot schimba versiunile pachetelor cursului,
  iar recompilarea durează mult.
