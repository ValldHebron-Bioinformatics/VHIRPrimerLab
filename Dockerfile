# vhirprimerlab web app. Also the image the Hugging Face Space builds.
#   docker build -t vhirprimerlab-web .
#   docker run --rm -p 7860:7860 vhirprimerlab-web     # then open http://localhost:7860
FROM python:3.11-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends mafft \
 && rm -rf /var/lib/apt/lists/*

# Hugging Face Spaces run the container as user 1000
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    VHIRPRIMERLAB_JOBS_DIR=/tmp/vhirprimerlab-jobs
WORKDIR /home/user/app

COPY --chown=user pyproject.toml ./
COPY --chown=user vhirprimerlab ./vhirprimerlab
RUN pip install --no-cache-dir --user ".[web]"

EXPOSE 7860
CMD ["uvicorn", "vhirprimerlab.web.app:app", "--host", "0.0.0.0", "--port", "7860"]
