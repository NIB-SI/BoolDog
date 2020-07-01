# docker-32bit-java6
FROM 32bit/debian

# Default to UTF-8 file.encoding
ENV LANG C.UTF-8

ENV JAVA_FOLDER /usr/local/openjdk-6
ENV JAVA_HOME $JAVA_FOLDER/jre1.6.0_45
ENV PATH $JAVA_HOME/bin:$PATH

ENV SQUAD_HOME /squad


WORKDIR $JAVA_FOLDER
COPY ./libraries/jre-6u45-linux-i586.bin .
RUN chmod a+x jre-6u45-linux-i586.bin && \
     ./jre-6u45-linux-i586.bin -yi && \
     rm jre-6u45-linux-i586.bin

RUN echo "deb http://ftp.de.debian.org/debian jessie main" | tee -a /etc/apt/sources.list && \
    apt-get update && \
    apt-get install libxp6 libxtst6 libstdc++5 libc6 -y --force-yes 

COPY ./libraries/compat-libstdc++-296_2.96-147.1_i386.deb .
RUN dpkg -i compat-libstdc++-296_2.96-147.1_i386.deb && \
    rm compat-libstdc++-296_2.96-147.1_i386.deb 

WORKDIR /home
RUN mkdir log && \
    pwd -P && \
    ls && \
    java -version 

VOLUME $SQUAD_HOME
COPY squad2-2/ $SQUAD_HOME/

ENTRYPOINT ["/squad/Squad"]
