FROM quay.io/pypa/manylinux2014_x86_64

# Step 1: build LLVM
RUN yum -y update
RUN yum -y install devtoolset-7 ninja-build cmake3 bzip2-devel
RUN mkdir -p /opt/llvm-codon
RUN git clone --depth 1 -b codon https://github.com/exaloop/llvm-project /github/llvm-src
RUN scl enable devtoolset-7 -- \
    cmake3 -S /github/llvm-src/llvm -G Ninja \
           -B /github/llvm-src/llvm/build \
           -DCMAKE_BUILD_TYPE=Release \
           -DLLVM_INCLUDE_TESTS=OFF \
           -DLLVM_ENABLE_RTTI=ON \
           -DLLVM_ENABLE_ZLIB=OFF \
           -DLLVM_ENABLE_TERMINFO=OFF \
           -DLLVM_TARGETS_TO_BUILD="host" \
           -DLLVM_BUILD_TOOLS=OFF \
           -DLLVM_ENABLE_PROJECTS=clang \
           -DCMAKE_INSTALL_PREFIX=/opt/llvm-codon
RUN scl enable devtoolset-7 -- cmake3 --build /github/llvm-src/llvm/build
RUN cmake3 --install /github/llvm-src/llvm/build

# RUN cd /github/llvm-src && tar cjvf /opt/llvm-codon-$(git rev-parse --short HEAD).tar.bz2 -C /opt llvm-codon/
# CMD cp /opt/llvm-codon-*.tar.bz2 /mnt/
# COPY llvm-codon-55b0b8fa1.tar.bz2 /opt
# RUN tar xf /opt/llvm-codon-55b0b8fa1.tar.bz2 -C /opt

# Step 2: build Codon & Seq
RUN git clone -b develop https://github.com/exaloop/codon /github/codon
RUN cmake3 -S /github/codon -B /github/codon/build \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=/opt/llvm-codon/bin/clang \
    -DCMAKE_CXX_COMPILER=/opt/llvm-codon/bin/clang++ \
    -DLLVM_DIR=/opt/llvm-codon/lib/cmake/llvm \
    -DCMAKE_INSTALL_PREFIX=/opt/codon
RUN cmake3 --build /github/codon/build
RUN cmake3 --install /github/codon/build

RUN git clone -b develop https://github.com/exaloop/seq /github/seq
RUN cmake3 -S /github/seq -B /github/seq/build \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=/opt/llvm-codon/bin/clang \
    -DCMAKE_CXX_COMPILER=/opt/llvm-codon/bin/clang++ \
    -DLLVM_DIR=/opt/llvm-codon/lib/cmake/llvm \
    -DCODON_PATH=/opt/codon \
    -DCMAKE_INSTALL_PREFIX=/opt/codon/lib/codon/plugins/seq
RUN cmake3 --build /github/seq/build
RUN cmake3 --install /github/seq/build

# RUN cd /github/codon && tar cjvf /opt/codon-$(git rev-parse --short HEAD).tar.bz2 -C /opt codon/
# RUN cp /opt/codon-*.tar.bz2 /mnt/

RUN git clone -b master https://github.com/0xTCG/biser /github/biser
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade twine setuptools wheel
RUN cd /github/biser && PATH=/opt/codon/bin:${PATH} python3 setup.py sdist bdist_wheel -p manylinux2014-x86_64
CMD cp -r /github/biser/dist /mnt/dist
# twine upload dist/*
