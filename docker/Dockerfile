FROM laristra/flecsi-third-party:fedora


ARG COVERAGE
ARG SONARQUBE
ARG SONARQUBE_TOKEN
ARG SONARQUBE_GITHUB_TOKEN
ARG CC

#for coverage
ARG CI
ARG TRAVIS
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT
ARG TRAVIS_OS_NAME
# Cleans and install 

RUN rm -rf /home/flecsi/.ccache
COPY flecsph /home/flecsi/flecsph
COPY ccache/ /home/flecsi/.ccache
COPY sonar/ /home/flecsi/.sonar
USER root
RUN chown -R flecsi:flecsi /home/flecsi/flecsph /home/flecsi/.ccache /home/flecsi/.sonar
RUN yum install -y which; exit 0
RUN yum install -y gsl-devel; exit 0 
RUN yum install -y gtest-devel; exit 0
RUN sed -i "s%http://archive.ubuntu.com%http://old-releases.ubuntu.com%g" /etc/apt/sources.list; exit 0
RUN apt-get update -y; exit 0
RUN apt-get -y install gsl-bin libgsl0-dev; exit 0 
RUN apt-get install -y libgtest-dev; exit 0
RUN apt-get -y install software-properties-common; exit 0
RUN apt-get install wget; exit 0
RUN apt-get install libssl-dev; exit 0
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.0-rc1/cmake-3.18.0-rc1.tar.gz; tar zxf cmake-3.18.0-rc1.tar.gz; cd cmake-3.18.0-rc1/; cmake .; make; make install; exit 0
#RUN apt-get install -y cmake; exit 0

#build flecsi
RUN cd /home/flecsi && \ 
    git clone --depth 1 --recursive https://github.com/laristra/flecsi flecsi && \
    cd flecsi && mkdir build && cd build && \ 
    cmake .. -DFLECSI_RUNTIME_MODEL=mpi \
              -DENABLE_LOG=OFF \
              -DENABLE_MPI=ON \
              -DENABLE_OPENMP=ON \
              -DENABLE_MPI_CXX_BINDINGS=ON \
              -DENABLE_CONFORMANCE_STANDARD=c++17 \
              -DLegion_ROOT=/usr/local \
              -DCMAKE_INSTALL_PREFIX=/usr/local \
              -DENABLE_BOOST_PREPROCESSOR=ON  \
              -DENABLE_FLECSIT=OFF \
              -DENABLE_FLECSI_TUTORIAL=OFF   && \ 
    make -j4 && \
    make install 

# Buidl FleCSPH 

ENV LD_LIBRARY_PATH="/usr/local/lib64/:/usr/local/lib/:${LD_LIBRARY_PATH}"
ENV CMAKE_PREFIX_PATH="/usr/local/lib:/usr/local/lib64:${CMAKE_PREFIX_PATH}"

USER flecsi 
RUN cd /home/flecsi/flecsph && \
    mkdir build && cd build && \
    ccache -z && \
    cmake -DCMAKE_BUILD_TYPE=Debug \
          -DENABLE_MPI=ON \
          -DENABLE_UNIT_TESTS=ON \
          -DENABLE_MPI_TESTS=OFF \
          -DENABLE_OPENMP=ON \
          -DENABLE_DOXYGEN=ON \
          -DCMAKE_CXX_FLAGS="-fpermissive" \
          -DCXX_CONFORMANCE_STANDARD=c++17 \
          -DENABLE_BOOST_PREPROCESSOR=ON \
          -DENABLE_LOG=ON \
          -DENABLE_COLOR_UNIT_TESTS=ON \
          -DFleCSI_INCLUDE_DIR=/usr/local/include \
          -DFleCSI_RUNTIME=/usr/local/share/FleCSI/runtime \
          -DENABLE_MPI_THREAD_MULITPLE=ON \
          -DMPIEXEC_PREFLAGS=--oversubscribe \
          ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} .. && \
    HDF5_USE_SHLIB=yes ${SONARQUBE:+build-wrapper-linux-x86-64 --out-dir bw-output} make -j2 && \
    ccache -s && \
    make install DESTDIR=${PWD}/install && rm -rf ${PWD}/install  && \
    make test

# COVERAGE & SONARQUE 
WORKDIR /home/flecsi/flecsph
RUN if [ ${COVERAGE} ]; then \
      $HOME/.local/bin/codecov -F "${CC}"; \
    fi 
WORKDIR /home/flecsi
