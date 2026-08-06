#!/usr/bin/env bash

set -euo pipefail



L_QDRANT=16333

L_NEO4J_HTTP=17474

L_NEO4J_BOLT=17687

PIDFILE="${HOME}/.tunnel_stack.pid"



usage() { echo "Usage: $0 start <JOBID> | stop | status"; }



get_node() {

  local jobid="$1"

  local node

  node="$(squeue -j "$jobid" -h -o "%N" | head -n1)"

  [[ -n "$node" ]] || { echo "ERROR: JOBID $jobid not running?"; exit 1; }

  echo "$node"

}



is_running() {

  [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null

}



start() {

  local jobid="$1"

  local node

  node="$(get_node "$jobid")"



  if is_running; then

    echo "Tunnel already running (PID $(cat "$PIDFILE"))."

    exit 0

  fi



  echo "Tunneling to NODE=$node"

  echo "Qdrant:  http://localhost:${L_QDRANT}  (-> ${node}:6333)"

  echo "Neo4j:   http://localhost:${L_NEO4J_HTTP} (-> ${node}:7474)"

  echo "Bolt:    bolt://localhost:${L_NEO4J_BOLT} (-> ${node}:7687)"



  ssh -4 -N \
    -o ExitOnForwardFailure=yes \
    -o ServerAliveInterval=60 \
    -o ServerAliveCountMax=3 \
    -o TCPKeepAlive=yes \
    -L ${L_QDRANT}:localhost:6333 \
    -L ${L_NEO4J_HTTP}:localhost:7474 \
    -L ${L_NEO4J_BOLT}:localhost:7687 \
    fnatali@${node} &

  echo $! > "$PIDFILE"

  echo "Tunnel started (PID $(cat "$PIDFILE"))."

}



stop() {

  if ! [[ -f "$PIDFILE" ]]; then

    echo "No tunnel PIDFILE found."

    exit 0

  fi

  local pid

  pid="$(cat "$PIDFILE")"

  kill "$pid" 2>/dev/null || true

  rm -f "$PIDFILE"

  echo "Tunnel stopped."

}



status() {

  if is_running; then

    echo "Tunnel running (PID $(cat "$PIDFILE"))."

    lsof -iTCP:${L_QDRANT} -sTCP:LISTEN -nP 2>/dev/null || true

    lsof -iTCP:${L_NEO4J_HTTP} -sTCP:LISTEN -nP 2>/dev/null || true

    lsof -iTCP:${L_NEO4J_BOLT} -sTCP:LISTEN -nP 2>/dev/null || true

  else

    echo "Tunnel NOT running."

  fi

}



cmd="${1:-}"

case "$cmd" in

  start) [[ $# -eq 2 ]] || { usage; exit 1; }; start "$2" ;;

  stop) stop ;;

  status) status ;;

  *) usage; exit 1 ;;

esac

