# Dockerfile, Image, Container

FROM python:3.8

WORKDIR /arraylib

# Install system dependencies
RUN apt-get update && \
    apt-get install -y bowtie2 

RUN pip install arraylib-solve

COPY . /arraylib

CMD ["/bin/bash" ]