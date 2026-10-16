#!/bin/sh
# pluto [notebook.jl]  — pornește serverul Pluto (implicit)
# run <notebook.jl>    — rulează notebook-ul ca script, pentru a verifica erorile
# orice altă comandă   — este executată ca atare (ex.: julia, bash)
set -e

# Ecran virtual pentru GLMakie (xvfb-run se blochează cînd rulează ca PID 1).
Xvfb :99 -screen 0 1280x1024x24 -nolisten tcp >/dev/null 2>&1 &
export DISPLAY=:99

case "$1" in
    pluto)
        shift
        nb=""
        [ -n "$1" ] && nb=", notebook=\"$1\""
        exec julia --project=/opt/curs -e "
            using Pluto
            Pluto.run(; host=\"0.0.0.0\", port=1234, launch_browser=false$nb)"
        ;;
    run)
        shift
        exec julia --project=/notebooks "$@"
        ;;
    *)
        exec "$@"
        ;;
esac
