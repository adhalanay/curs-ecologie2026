# Imagine single-user pentru JupyterHub: JupyterLab + kernel Julia + Pluto (prin jupyter-server-proxy).
#
#   docker build -t curs-ecologie2026 .
#   docker run --rm -p 8888:8888 curs-ecologie2026      # test local, fără JupyterHub
#
# În JupyterHub: c.DockerSpawner.image = "curs-ecologie2026" (sau KubeSpawner.image).
# Notebook-urile sînt copiate în ~/course la pornire (fără să suprascrie ce a modificat studentul).

FROM quay.io/jupyter/base-notebook:python-3.11.10

USER root

# Xvfb + Mesa: GLMakie (Test1.jl, Cursul10.jl) are nevoie de un context OpenGL.
# CairoMakie își aduce bibliotecile prin JLL-uri, nu are nevoie de pachete -dev.
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates git \
        xvfb xauth libgl1 libglx-mesa0 libgl1-mesa-dri \
        fonts-liberation fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

ARG JULIA_VERSION=1.13.2
ARG JULIA_SHA256=2ecdfc9890ab0a565fdd92404c915ba2a04a933dbd5d292cc224d13819d17d07

RUN curl -fsSL -o /tmp/julia.tar.gz \
        "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
    && echo "${JULIA_SHA256}  /tmp/julia.tar.gz" | sha256sum -c - \
    && mkdir /opt/julia \
    && tar -xzf /tmp/julia.tar.gz -C /opt/julia --strip-components=1 \
    && rm /tmp/julia.tar.gz \
    && ln -s /opt/julia/bin/julia /usr/local/bin/julia

RUN pip install --no-cache-dir jupyter-server-proxy \
    && fix-permissions "${CONDA_DIR}"

# Depozitul Julia e în afara lui $HOME, fiindcă JupyterHub montează un volum peste /home/jovyan.
ENV JULIA_PROJECT=/opt/course/julia \
    JULIA_DEPOT_PATH=/opt/course/depot \
    JUPYTER=${CONDA_DIR}/bin/jupyter \
    IJULIA_NODEFAULTKERNEL=1 \
    JULIA_CPU_TARGET="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"
# JULIA_CPU_TARGET: codul precompilat merge și pe procesorul serverului JupyterHub, nu doar pe cel de build.

# Pachetele sînt descărcate și precompilate o singură dată, la build (durează 15–30 min).
# fix-permissions în același RUN, altfel stratul cu depozitul (cîțiva GB) ar fi duplicat.
COPY Project.toml Manifest.toml /opt/course/julia/
RUN JUPYTER_DATA_DIR="${CONDA_DIR}/share/jupyter" xvfb-run -a julia -e ' \
        using Pkg; Pkg.instantiate(); Pkg.precompile(); \
        using IJulia; IJulia.installkernel("Ecologie 2026")' \
    && fix-permissions /opt/course "${CONDA_DIR}/share/jupyter"

# Lansatorul Pluto din JupyterLab (jupyter-server-proxy). Trebuie să fie un fișier
# jupyter_server_config.json; JSON-urile din jupyter_server_config.d/ nu sînt citite pentru ServerProxy.
COPY pluto.json /usr/local/etc/jupyter/jupyter_server_config.json

# Materialele cursului; vezi .dockerignore pentru ce este exclus (partial.jl, PDF-uri, .xcf).
COPY Project.toml Manifest.toml *.jl *.png *.PNG *.jpg /opt/course/notebooks/
COPY course-materials.sh /usr/local/bin/before-notebook.d/course-materials.sh
RUN chmod a+rx /usr/local/bin/before-notebook.d/course-materials.sh \
    && fix-permissions /opt/course/notebooks

USER ${NB_UID}
WORKDIR /home/${NB_USER}
