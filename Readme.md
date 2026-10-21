Aveți nevoie de:

1. [Git](https://git-scm.com/) pentru a clona codul de față;
2. Limbajul [Julia](https://julialang.org/);

După instalare, din directorul proiectului rulați (după fiecare linie apăsați **Enter**)
 
 ````
    julia
    ]
    activate .
    instantiate
````
:arrow_left: (Backspace)

````
    using Pluto
    Pluto.run() 
````

## Varianta cu Docker

Dacă aveți [Docker](https://www.docker.com/) instalat, nu mai este nevoie de Julia. Din directorul proiectului:

````
    docker compose build
    docker compose up
````

Apoi deschideți în browser `http://localhost:8888` (JupyterLab); Pluto se pornește din Launcher, iar notebook-urile sînt în directorul `course`.
Ele sînt aduse de pe GitHub la fiecare pornire (e nevoie de internet), deci în container apare ultima versiune publicată, nu copia locală.

## Serverul JupyterHub al cursului

Pentru lucrul pe serverul JupyterHub (unde se păstrează lucrul vostru, cum folosiți pachetele cursului
în Pluto și Jupyter, cum instalați pachete suplimentare), citiți [Instructiuni-studenti.md](Instructiuni-studenti.md).

Un notebook poate fi rulat și fără browser, pentru a verifica dacă are erori
(se folosește versiunea din directorul proiectului, nu copia din imagine):

````
    docker compose run --rm -v .:/work -w /work curs xvfb-run -a julia --project=. Cursul5.jl
````
