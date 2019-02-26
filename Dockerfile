# DockerFile for the HYPE model.
# The input data root directory should be mounted as /data
# The filedir.txt in the input directory should only contain the single line: /data/
# The output folder should exist within the mounted volume
FROM ubuntu:bionic
MAINTAINER Gijs van den Oord <g.vandenoord@esciencecenter.nl>
RUN apt-get update
RUN apt-get install -y wget git build-essential g++ make gfortran
RUN git clone https://github.com/eWaterCycle/grpc4bmi.git /opt/grpc4bmi
COPY . /opt/hype-bmi
WORKDIR /opt/hype-bmi/src
RUN make CPPFLAGS="-I/opt/grpc4bmi/cpp"
VOLUME /data
WORKDIR /data
CMD ["/opt/hype-bmi/src/hypec","/data/"]
