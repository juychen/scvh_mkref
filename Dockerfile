FROM ubuntu:20.04

LABEL maintainer="junyichen@d24h.hk"

# Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

ENV DEBIAN_FRONTEND noninteractive
ENV PATH=/opt/samtools/bin:$PATH
COPY . /code

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	    apt-utils \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		apt-transport-https \
		gsfonts \
		gnupg2 \
        build-essential \
        libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev tini \
	&& rm -rf /var/lib/apt/lists/*

# Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -xf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure --prefix=/opt/samtools && \
    make && \
    make install && \
    echo "export PATH=/opt/samtools/bin:\$PATH" >> ~/.bashrc && \
    . ~/.bashrc && \
    chmod 755 /code

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

ENTRYPOINT [ "/usr/bin/tini","--"]
CMD [ "/bin/bash" ]
