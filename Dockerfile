# GUIDANCE3 command-line tool in a container: runs with Docker or Podman, nothing else to install.
#
#   docker build -t guidance3 .
#   docker run --rm -v "$PWD":/data guidance3 --seqFile /data/seqs.fasta --msaProgram MAFFT --seqType aa \
#              --outDir /data/guidance3_out --proc_num 4
#
# The image is built natively for the machine that builds it:
#   linux/amd64 (Intel/AMD PCs and servers, Windows): uses the bundled x86-64 Linux binaries in guidance3/programs/linux.
#   linux/arm64 (Apple-silicon Macs, ARM servers): compiles the helper programs from src/ and uses the official
#               IQ-TREE ARM build, so nothing runs under emulation.
ARG IQTREE_VERSION=3.1.4

# ---- per-architecture helper programs -------------------------------------------------------------
FROM python:3.10-slim AS tools-amd64
RUN mkdir /tools

FROM python:3.10-slim AS tools-arm64
ARG IQTREE_VERSION
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential cmake curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY src /src
RUN sh /src/build_tools.sh /tools \
 && curl -fsSL "https://github.com/iqtree/iqtree3/releases/download/v${IQTREE_VERSION}/iqtree-${IQTREE_VERSION}-Linux-arm.tar.gz" \
    | tar -xz -C /tmp \
 && cp /tmp/iqtree-${IQTREE_VERSION}-Linux-arm/bin/iqtree3 /tools/iqtree3

# TARGETARCH is set automatically by Docker BuildKit (amd64 or arm64); only the matching stage is built.
FROM tools-${TARGETARCH} AS tools

# ---- final image -----------------------------------------------------------------------------------
FROM python:3.10-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends mafft prank clustalo libgomp1 \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/guidance3
COPY pyproject.toml README.md ./
COPY guidance3 ./guidance3
RUN pip install --no-cache-dir . \
 && cd / && rm -rf /opt/guidance3

# On arm64, replace the bundled x86-64 binaries with the native ones built above.
WORKDIR /
COPY --from=tools /tools /tools
RUN P=$(python -c "import guidance3, os; print(os.path.join(os.path.dirname(guidance3.__file__), 'programs', 'linux'))") \
 && if [ -f /tools/iqtree3 ]; then \
        cp /tools/msa_set_score /tools/isEqualTree /tools/removeTaxa /tools/features_for_msas "$P/" \
     && cp /tools/iqtree3 "$P/iqtree/bin/iqtree3" \
     && ln -sf /usr/bin/clustalo "$P/clustalo"; \
    fi \
 && rm -rf /tools

WORKDIR /data
ENTRYPOINT ["guidance3"]
CMD ["--help"]
