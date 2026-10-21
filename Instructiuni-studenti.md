# Lucrul pe serverul JupyterHub al cursului

## Unde se păstrează lucrul vostru

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
celula următoare și salvați notebook-ul în `~/course`:

```julia
begin
    using Pkg
    Pkg.activate(joinpath(homedir(), "course"))
    using CairoMakie, PlutoUI, DifferentialEquations
end
```

Fără `Pkg.activate`, Pluto creează un mediu separat pentru notebook și poate descărca și compila
din nou pachetele, ceea ce poate dura multe minute.

### În Jupyter

Alegeți kernel-ul **Ecologie 2026 1.13**; `using CairoMakie` etc. funcționează imediat.

## Pachete suplimentare

Instalați-le într-un mediu propriu, în directorul home:

```julia
using Pkg
Pkg.activate(joinpath(homedir(), "mywork"))
Pkg.add("NumePachet")
```

În sesiunile următoare activați din nou același mediu (`Pkg.activate(joinpath(homedir(), "mywork"))`)
și pachetul este disponibil.

- Nu folosiți `Pkg.add` direct din kernel-ul Jupyter, fără `Pkg.activate`: mediul acelui kernel
  nu se păstrează, iar după repornire pachetul nu mai apare.
- Nu adăugați pachete în mediul din `~/course`: se pot schimba versiunile pachetelor cursului,
  iar recompilarea durează mult.
