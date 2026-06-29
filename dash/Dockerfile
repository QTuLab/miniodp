# syntax=docker/dockerfile:1
FROM python:3.11-slim

ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

LABEL org.opencontainers.image.revision=$VCS_REF \
      org.opencontainers.image.created=$BUILD_DATE

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies required for building Python packages
# libcurl4-openssl-dev and zlib1g-dev are strictly required for pyBigWig
# build-essential is needed for compiling C extensions
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libhdf5-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements file first to leverage Docker cache layers
COPY requirements.txt .

# Install Python dependencies including gunicorn for production serving
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir gunicorn

# Runtime configuration defaults
ENV PORT=3939 \
    DASH_HOST=0.0.0.0 \
    DASH_DEBUG=false \
    GUNICORN_WORKERS=1 \
    GUNICORN_THREADS=2 \
    GUNICORN_TIMEOUT=180 \
    GUNICORN_MAX_REQUESTS=0 \
    GUNICORN_MAX_REQUESTS_JITTER=0 \
    GUNICORN_ACCESS_LOG=/app/logs/dash-access.log \
    GUNICORN_ERROR_LOG=/app/logs/dash-app.log

EXPOSE 3939

# Use gunicorn as the production server
# Workers set to 2 for basic concurrency; adjust based on CPU cores
# Logs directed to /app/logs mapped volume
CMD ["/bin/sh", "-c", "mkdir -p \"$(dirname ${GUNICORN_PID_FILE:-/app/logs/gunicorn.pid})\" && python scripts/memory_watchdog.py & exec gunicorn app:server -b 0.0.0.0:${PORT:-3939} --workers=${GUNICORN_WORKERS:-1} --threads=${GUNICORN_THREADS:-2} --timeout=${GUNICORN_TIMEOUT:-180} --max-requests=${GUNICORN_MAX_REQUESTS:-0} --max-requests-jitter=${GUNICORN_MAX_REQUESTS_JITTER:-0} --pid ${GUNICORN_PID_FILE:-/app/logs/gunicorn.pid} --access-logfile ${GUNICORN_ACCESS_LOG} --error-logfile ${GUNICORN_ERROR_LOG} --capture-output"]
