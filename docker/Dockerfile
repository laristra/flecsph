FROM laristra/flecsi-third-party:fedora

ARG MPI
ARG COVERAGE
ARG SONARQUBE
ARG SONARQUBE_TOKEN
ARG SONARQUBE_GITHUB_TOKEN
ARG CC
ARG CXX
ARG CXXFLAGS

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

ENV CLOG_ENABLE_STDLOG=0
RUN rm -rf /home/flecsi/.ccache
COPY flecsph /home/flecsi/flecsph
COPY ccache/ /home/flecsi/.ccache
COPY sonar/ /home/flecsi/.sonar
USER root
RUN chown -R flecsi:flecsi /home/flecsi/flecsph /home/flecsi/.ccache /home/flecsi/.sonar
RUN yum install -y which; exit 0
USER flecsi

#build flecsi
WORKDIR /home/flecsi
RUN git clone -b feature/flecsph --depth 1 --recursive https://github.com/laristra/flecsi flecsi
RUN  cd flecsi && mkdir build && cd build
WORKDIR flecsi/build
RUN  cmake .. -DFLECSI_RUNTIME_MODEL=mpi \
              -DENABLE_CLOG=OFF \
              -DENABLE_MPI=ON \
              -DENABLE_OPENMP=ON \
              -DENABLE_MPI_CXX_BINDINGS=ON \
              -DENABLE_CONFORMANCE_STANDARD=c++17 \
              -DLegion_ROOT=/usr/local \
              -DCMAKE_INSTALL_PREFIX=/usr/local \
              -DENABLE_BOOST_PREPROCESSOR=ON  \
              -DENABLE_FLECSIT=OFF \
              -DENABLE_FLECSI_TUTORIAL=OFF
RUN make VERBOSE=1 -j4
USER root
RUN make install
USER flecsi

WORKDIR /home/flecsi/flecsph
RUN mkdir build

WORKDIR build
RUN ccache -z
# FleCSI lib location regarding distro
RUN cmake -DENABLE_MPI=ON \
          -DENABLE_UNIT_TESTS=ON \
          -DENABLE_OPENMP=ON \
          -DENABLE_DOXYGEN=ON \
          -DCMAKE_CXX_FLAGS="-fpermissive" \
          -DCXX_CONFORMANCE_STANDARD=c++17 \
          -DENABLE_BOOST_PREPROCESSOR=ON \
          -DENABLE_CLOG=ON \
          -DENABLE_COLOR_UNIT_TESTS=ON \
          -DFleCSI_INCLUDE_DIR=/usr/local/include \
          -DFleCSI_RUNTIME=/usr/local/share/FleCSI/runtime \
          -DFleCSI_LIBRARY=/usr/local/lib*/libFleCSI.so \
          -DENABLE_MPI_THREAD_MULITPLE=ON \
          ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} ..

RUN HDF5_USE_SHLIB=yes ${SONARQUBE:+build-wrapper-linux-x86-64 --out-dir bw-output} make VERBOSE=1 -j2
RUN ccache -s
RUN make doxygen
RUN make install DESTDIR=${PWD}/install && rm -rf ${PWD}/install
ENV LD_LIBRARY_PATH="/usr/local/lib64/:${LD_LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"
RUN echo $LD_LIBRARY_PATH
RUN make test ARGS="-V"
RUN cd .. && if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov -F "${CC}" --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov -F "${CC}"; \
  fi; \
fi && cd -
RUN cd .. && if [ ${SONARQUBE} ]; then \
  sonar-scanner \
    -Dsonar.projectKey=${TRAVIS_REPO_SLUG##*/} \
    -Dsonar.projectName=${TRAVIS_REPO_SLUG#*/} \
    -Dsonar.projectVersion=${TRAVIS_COMMIT} \
    -Dsonar.branch=/${TRAVIS_BRANCH} \
    -Dsonar.links.homepage=http://${TRAVIS_REPO_SLUG%%/*}.github.io/${TRAVIS_REPO_SLUG#*/} \
    -Dsonar.links.ci=https://travis-ci.org/${TRAVIS_REPO_SLUG} \
    -Dsonar.links.scm=https://github.com/${TRAVIS_REPO_SLUG} \
    -Dsonar.links.issue=https://github.com/${TRAVIS_REPO_SLUG}/issues \
    -Dsonar.sources=. \
    -Dsonar.exclusions=build/CMakeFiles/**,cinch/**,build/doc/doxygen/**,build/include/** \
    -Dsonar.cfamily.build-wrapper-output=build/bw-output \
    -Dsonar.cfamily.gcov.reportsPath=. \
    -Dsonar.host.url=https://sonarcloud.io \
    -Dsonar.organization=${TRAVIS_REPO_SLUG%%/*} \
    $([ ${TRAVIS_PULL_REQUEST} != false ] && echo \
       -Dsonar.analysis.mode=preview \
       -Dsonar.github.pullRequest=${TRAVIS_PULL_REQUEST} \
       -Dsonar.github.repository=${TRAVIS_REPO_SLUG} \
       -Dsonar.github.oauth=${SONARQUBE_GITHUB_TOKEN}) \
    -Dsonar.login=${SONARQUBE_TOKEN}; \
fi && cd -
WORKDIR /home/flecsi
