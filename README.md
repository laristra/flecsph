# SPH on FleCSI 

This project is an implementation of SPH problem using FleCSI framework.

# Getting the Code 

    % git clone --recursive git@gitlab.lanl.gov/jloiseau/FleCSI_SPH.git

# Build 

    % mkdir build
    % cd build 
    % ccmake ../ 
    % make 

# Known problems 

When FleCSI is make and installed I needed to add manually some includes:
- helper.h 
- future.h
- thread_pool.h
- virutal_semaphore.h