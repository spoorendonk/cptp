# rcspp-bac build + local LLM analysis image
# Runs on CPU-only GitHub runners (no GPU required)
#
# Build:  docker build -t rcspp-llm -f ci/Dockerfile .
# Run:    docker run --rm rcspp-llm

# ── Stage 1: build llama.cpp ─────────────────────────────────────────
FROM ubuntu:24.04 AS llama-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        git cmake g++ make ca-certificates && \
    rm -rf /var/lib/apt/lists/*

ARG LLAMA_CPP_VERSION=b5200
RUN git clone --depth 1 --branch ${LLAMA_CPP_VERSION} \
        https://github.com/ggerganov/llama.cpp /llama.cpp && \
    cmake -B /llama.cpp/build -S /llama.cpp \
        -DCMAKE_BUILD_TYPE=Release \
        -DGGML_NATIVE=OFF \
        -DLLAMA_CURL=OFF \
        -DLLAMA_SERVER=ON \
        -DBUILD_SHARED_LIBS=OFF && \
    cmake --build /llama.cpp/build -j$(nproc) --target llama-server llama-cli

# ── Stage 2: runtime image ───────────────────────────────────────────
FROM ubuntu:24.04

# Build toolchain for the project
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc-14 g++-14 cmake make libtbb-dev git \
        curl jq python3 python3-pip ca-certificates && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 && \
    rm -rf /var/lib/apt/lists/*

# Copy llama.cpp binaries
COPY --from=llama-builder /llama.cpp/build/bin/llama-server /usr/local/bin/
COPY --from=llama-builder /llama.cpp/build/bin/llama-cli    /usr/local/bin/

# Download Qwen2.5-Coder-1.5B GGUF (Q4_K_M, ~1GB)
ARG MODEL_URL=https://huggingface.co/Qwen/Qwen2.5-Coder-1.5B-Instruct-GGUF/resolve/main/qwen2.5-coder-1.5b-instruct-q4_k_m.gguf
RUN mkdir -p /models && \
    curl -fSL -o /models/model.gguf "${MODEL_URL}"

# Copy analysis scripts
COPY ci/analyze.sh     /usr/local/bin/analyze.sh
COPY ci/llm-review.sh  /usr/local/bin/llm-review.sh
RUN chmod +x /usr/local/bin/analyze.sh /usr/local/bin/llm-review.sh

WORKDIR /workspace

ENTRYPOINT ["/usr/local/bin/analyze.sh"]
