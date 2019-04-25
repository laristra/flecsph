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
ENV CLOG_ENABLE_STDLOG=0
USER root
COPY flecsph /home/flecsi/flecsph
RUN rm -rf /home/flecsi/.ccache && \
    #cp -r laristra/flecsph /home/flecsi/ && \
    cp -r ccache/ /home/flecsi/.ccache && \
    cp -r sonar/ /home/flecsi/.sonar && \
    chown -R flecsi:flecsi /home/flecsi/flecsph /home/flecsi/.ccache /home/flecsi/.sonar && \
    yum install -y which; exit 0

#build flecsi
RUN cd /home/flecsi && \ 
    git clone -b feature/flecsph --depth 1 --recursive https://github.com/laristra/flecsi flecsi && \
    cd flecsi && mkdir build && cd build && \ 
    cmake .. -DFLECSI_RUNTIME_MODEL=mpi \
              -DENABLE_CLOG=OFF \
              -DENABLE_MPI=ON \
              -DENABLE_OPENMP=ON \
              -DENABLE_MPI_CXX_BINDINGS=ON \
              -DENABLE_CONFORMANCE_STANDARD=c++17 \
              -DLegion_ROOT=/usr/local \
              -DCMAKE_INSTALL_PREFIX=/usr/local \
              -DENABLE_BOOST_PREPROCESSOR=ON  \
              -DENABLE_FLECSIT=OFF \
              -DENABLE_FLECSI_TUTORIAL=OFF   && \ 
    make VERBOSE=1 -j4 && \
    make install 

# Buidl FleCSPH 
RUN cd /home/flecsi/flecsph && \
    mkdir build && cd build && \
    ccache -z && \
    cmake -DENABLE_MPI=ON \
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
          ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} .. && \
    HDF5_USE_SHLIB=yes ${SONARQUBE:+build-wrapper-linux-x86-64 --out-dir bw-output} make VERBOSE=1 -j2 && \
    ccache -s && \
    make doxygen && \
    make install DESTDIR=${PWD}/install && rm -rf ${PWD}/install && \

ENV LD_LIBRARY_PATH="/usr/local/lib64/:${LD_LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

# Tests and coverage 
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
