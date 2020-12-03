# This script contains commands to create docker images for multi-arch
# testing of MACS3 Note that The dependencies to compile and run MACS3
# have been installed in these images; however, MACS3 is not
# installed.

# base images are from Debian official images

# create Docker image for testing MACS3 in ppc64le
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="ppc64le/debian:buster" -t macs_ppc64le_base -f test/Dockerfile.multiarch.base.buster.py37 .

# create Docker image for testing MACS3 in i386
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="i386/debian:buster" -t macs_i386_base -f test/Dockerfile.multiarch.base.buster.py37 .

# create Docker image for testing MACS3 in arm32v7
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="arm32v7/debian:buster" -t macs_arm32v7_base -f test/Dockerfile.multiarch.base.buster.py37 .

# create Docker image for testing MACS3 in arm64v8
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="arm64v8/debian:buster" -t macs_arm64v8_base -f test/Dockerfile.multiarch.base.buster.py37 .

# create Docker image for testing MACS3 in s390x
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="s390x/debian:buster" -t macs_s390x_base -f test/Dockerfile.multiarch.base.buster.py37 .


