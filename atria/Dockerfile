FROM julia:1.8 AS build

ARG VERSION=4.1.1
RUN apt update -y && \
    apt upgrade -y && \
    apt install -y wget pigz pbzip2 gcc g++ && \
    apt-get clean
RUN wget https://github.com/cihga39871/Atria/archive/refs/tags/v${VERSION}.tar.gz && \
    tar -zxvf v${VERSION}.tar.gz && \
    cd Atria-${VERSION} && \
    julia build_atria.jl /app

FROM julia:1.8
ENV PATH=/app/atria-${VERSION}/bin:"$PATH"
COPY --from=build /app /app
RUN apt update -y && \
    apt upgrade -y && \
    apt install -y pigz pbzip2 && \
    apt-get clean
ENTRYPOINT ["atria"]