#!/usr/bin/env bash
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

PORT=${PORT:-5050}
PID=$(lsof -ti :$PORT 2>/dev/null || true)
if [ -n "$PID" ]; then
  echo "Freeing port $PORT (PID $PID)..."
  kill -9 $PID 2>/dev/null || true
  sleep 1
fi

export PORT=$PORT
# Configure Python SSL cert path
export SSL_CERT_FILE=$(python3 -c "import certifi; print(certifi.where())" 2>/dev/null || true)
export REQUESTS_CA_BUNDLE="$SSL_CERT_FILE"

echo "Using SSL Certificate Bundle: $SSL_CERT_FILE"
echo "Starting SeqLit web server on http://127.0.0.1:$PORT..."
exec python3 app.py
