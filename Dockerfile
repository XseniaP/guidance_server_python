# GUIDANCE3 command-line tool in a container: runs with Docker or Podman, nothing else to install.
#
#   docker build --platform linux/amd64 -t guidance3 .
#   docker run --rm -v "$PWD":/data guidance3 --seqFile /data/seqs.fasta --msaProgram MAFFT --seqType aa \
#              --outDir /data/guidance3_out --proc_num 4
#
# The bundled binaries (IQ-TREE 3, msa_set_score, ...) are x86-64 Linux builds, hence linux/amd64.
FROM python:3.10-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends mafft prank \
 && rm -rf /var/lib/apt/lists/*

# PyTorch is only needed for model training; install the CPU build so pip does not pull the CUDA wheels.
RUN pip install --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cpu \
        torch==2.3.0+cpu torchvision==0.18.0+cpu \
 && pip install --no-cache-dir torchtext==0.18.0

WORKDIR /opt/guidance3
COPY pyproject.toml README.md ./
COPY guidance3 ./guidance3
RUN pip install --no-cache-dir . \
 && rm -rf /opt/guidance3/build

WORKDIR /data
ENTRYPOINT ["guidance3"]
CMD ["--help"]
