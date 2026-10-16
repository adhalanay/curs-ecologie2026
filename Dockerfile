# Imagine Docker pentru notebook-urile Pluto ale cursului.
#
#   docker compose build
#   docker compose up                                 # Pluto pe http://localhost:1234
#   docker compose run --rm curs run Cursul5.jl       # rulează un notebook fără browser

FROM debian:trixie-slim

ARG JULIA_VERSION=1.13.2
ARG JULIA_SHA256=2ecdfc9890ab0a565fdd92404c915ba2a04a933dbd5d292cc224d13819d17d07

# Xvfb + Mesa: GLMakie (folosit în Test1.jl) are nevoie de un context OpenGL.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl xvfb xauth libgl1 libglx-mesa0 libgl1-mesa-dri \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL -o /tmp/julia.tar.gz \
        "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
    && echo "${JULIA_SHA256}  /tmp/julia.tar.gz" | sha256sum -c - \
    && mkdir /opt/julia \
    && tar -xzf /tmp/julia.tar.gz -C /opt/julia --strip-components=1 \
    && rm /tmp/julia.tar.gz \
    && ln -s /opt/julia/bin/julia /usr/local/bin/julia

ENV JULIA_DEPOT_PATH=/opt/julia-depot \
    HOME=/tmp

# Mediul Julia al cursului; pachetele sînt descărcate și precompilate o singură dată, la build.
WORKDIR /opt/curs
COPY Project.toml Manifest.toml ./
RUN xvfb-run -a julia --project=/opt/curs -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' \
    && chmod -R a+rwX /opt/julia-depot

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Notebook-urile se montează aici (vezi compose.yaml), ca modificările să ajungă pe disc.
WORKDIR /notebooks
EXPOSE 1234
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["pluto"]
