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


Un notebook poate fi rulat și fără browser, pentru a verifica dacă are erori:

````
    docker compose run --rm curs run Cursul5.jl
````
