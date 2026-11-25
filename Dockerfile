FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential git make \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.12.5 /uv /usr/local/bin/uv

WORKDIR /work
COPY pyproject.toml uv.lock .python-version README.md LICENSE ./
COPY src ./src
RUN uv sync --frozen --extra dev --no-install-project && uv sync --frozen --extra dev

COPY . .
ENV PATH="/work/.venv/bin:${PATH}"

CMD ["make", "quick"]
