# CPU-only image for running spapros (probe set selection and evaluation).
#
# Build:  docker build -t spapros .
# Run:    docker run --rm -v "$PWD:/work" spapros spapros evaluation ...
#    or:  docker run --rm -v "$PWD:/work" spapros python my_pipeline.py
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    # Headless plotting: spapros saves its figures to files.
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba

WORKDIR /opt/spapros
COPY pyproject.toml README.md LICENSE ./
COPY spapros ./spapros

# xgboost's Linux wheel pulls in nvidia-nccl-cu12 (~200 MB), which is only used for
# multi-GPU training and loaded lazily, so it is removed to keep the image CPU-only.
RUN pip install --upgrade pip \
    && pip install . \
    && pip uninstall -y nvidia-nccl-cu12 \
    && python -c "import spapros, xgboost; print('spapros', spapros.__version__, '| xgboost', xgboost.__version__)" \
    && rm -rf /tmp/*

RUN useradd --create-home --uid 1000 spapros
USER spapros
WORKDIR /work

CMD ["spapros", "--help"]
