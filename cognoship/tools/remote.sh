#!/bin/bash
# Run tools/models.py experiments on a GPU server over ssh (no rsync or node needed there).
#   tools/remote.sh push [output files...]   copy the Cognoship tree (tracked + untracked, not ignored)
#                                          and the named output/models files (e.g. a sweep that
#                                          --like reads) to the server, plus aquarium233.json (node's
#                                          reading of aquarium233.js, for a server without node)
#   tools/remote.sh run <name> <cmd...>      run a command there in the background (nohup), log in
#                                          output/models/<name>.log; the venv's python is $PY
#   tools/remote.sh log <name>               the log's last lines
#   tools/remote.sh pull <path...>           copy output/models paths back here (e.g. 5N7KKM/steer-limus-b.json)
# Server: $HOST (default r4), folder $DIR (default ~/claude/cognoship-runs), venv $DIR/.venv.
set -e
HOST=${HOST:-r4}
DIR=${DIR:-claude/cognoship-runs}
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
case "$1" in
  push)
    shift
    mkdir -p output/models
    node -e 'process.stdout.write(JSON.stringify(require(process.argv[1])))' "$HERE/aquarium233.js" > output/models/aquarium233.json
    { git ls-files -co --exclude-standard; echo output/models/aquarium233.json; for f in "$@"; do echo "output/models/$f"; done; } |
      tar -cf - -T - | ssh "$HOST" "mkdir -p ~/$DIR/Cognoship && tar -xf - -C ~/$DIR/Cognoship"
    echo "pushed to $HOST:~/$DIR/Cognoship" ;;
  run)
    name=$2; shift 2
    printf '%s\n' "PY=~/$DIR/.venv/bin/python" "$*" | ssh "$HOST" "cd ~/$DIR/Cognoship && mkdir -p output/models && cat > output/models/$name.sh"
    ssh -n "$HOST" "cd ~/$DIR/Cognoship && setsid -f bash output/models/$name.sh > output/models/$name.log 2>&1 < /dev/null; echo started $name" ;;
  log)
    ssh "$HOST" "tail -${3:-5} ~/$DIR/Cognoship/output/models/$2.log" ;;
  pull)
    shift
    ssh "$HOST" "cd ~/$DIR/Cognoship/output/models && tar -cf - $*" | tar -xf - -C output/models
    echo "pulled $*" ;;
  *) sed -n '2,12p' "$0"; exit 1 ;;
esac
