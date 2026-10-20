#!/bin/bash
# Hook before-notebook.d: copiază materialele cursului în ~/course, fără să suprascrie fișierele existente.
# Este *sourced* de start.sh, deci fără `set -e` / `exit`.

SRC="/opt/course/notebooks"
if [ "$(id -u)" = 0 ]; then
    DEST="/home/${NB_USER}/course"
else
    DEST="${HOME}/course"
fi

if [ -d "$SRC" ]; then
    mkdir -p "$DEST"
    cp -rn "$SRC"/. "$DEST"/ 2>/dev/null || true
    # Containerul pornit ca root: fișierele trebuie să fie ale utilizatorului.
    if [ "$(id -u)" = 0 ]; then
        chown -R "${NB_UID}:${NB_GID}" "$DEST"
    fi
fi
