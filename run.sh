#!/usr/bin/env bash
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

# Free port 5050 if in use
PORT=${PORT:-5050}
PID=$(lsof -ti :$PORT 2>/dev/null || true)
if [ -n "$PID" ]; then
  echo "Freeing port $PORT (PID $PID)..."
  kill -9 $PID 2>/dev/null || true
  sleep 1
fi

export PORT=$PORT
echo "Starting SeqLit web server on http://127.0.0.1:$PORT..."
exec python3 app.py
