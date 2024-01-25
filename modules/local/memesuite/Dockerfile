# Use an official base image
FROM biocontainers/biocontainers:v1.2.0_cv1

# Set metadata (optional)
LABEL maintainer="thondeboer@me.com"
LABEL description="Custom Bioinformatics Docker Image"

USER root

WORKDIR /opt
RUN apt-get update && apt-get install -y libexpat1-dev zlib1g-dev \
    && graphicsmagick ghostscript \
    wget -q https://meme-suite.org/meme/meme-software/5.5.5/meme-5.5.5.tar.gz && \
    tar -xvzf meme-5.5.5.tar.gz && \
    cd meme-5.5.5 && \
    ./configure --prefix=/usr/local --enable-build-libxml2 --enable-build-libxslt && \
    make && \
    make install && \
    conda install graphicsmagick
RUN cpan File::Which \
    HTML::PullParser \
    HTML::Template \
    HTML::TreeBuilder \
    JSON \
    XML::Simple \
    XML::Parser::Expat \
    Ssys::Info

# Set a default command or entry point (optional)
CMD ["bash"]