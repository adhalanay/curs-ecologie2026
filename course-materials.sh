#!/bin/bash
# Hook before-notebook.d: sincronizează ~/course cu depozitul cursului (course-sync), în fundal,
# ca serverul să pornească imediat. Jurnalul este în /tmp/course-sync.log.
# Este *sourced* de start.sh, deci fără `set -e` / `exit`.

if [ "$(id -u)" = 0 ]; then
    # Containerul pornit ca root: git trebuie rulat ca utilizatorul, altfel fișierele ar fi ale lui root.
    sudo -E -H -u "${NB_USER}" env PATH="${PATH}" course-sync >/tmp/course-sync.log 2>&1 &
else
    course-sync >/tmp/course-sync.log 2>&1 &
fi
