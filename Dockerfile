FROM python:3.11-slim

WORKDIR /app

# Install minimal build tools for C-extensions if needed
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source code
COPY . .

# Create writable data directory for cache and matplotlib
RUN mkdir -p data/matplotlib && chmod -R 777 data

# Default port (7860 for Hugging Face Spaces; dynamically overridden by $PORT on Render/Cloud Run)
ENV PORT=7860
EXPOSE 7860

# Universal entrypoint compatible with any host's $PORT
CMD ["sh", "-c", "gunicorn app:app --bind 0.0.0.0:${PORT} --timeout 120 --workers 2 --threads 4"]
