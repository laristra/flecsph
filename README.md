![logo](doc/flecsph_logo_bg.png)

# SPH on FleCSI 

This project is an implementation of SPH problem using FleCSI framework.

# Getting the Code 

    % git clone --recursive git@gitlab.lanl.gov/jloiseau/FleCSI_SPH.git

# Requirements

Before trying to compile you need to install on your system: 

- FleCSI thrird party
- FleCSI = compile legion

# Build 

    % mkdir build
    % cd build 
    % ccmake ../ 
    % make 


# Todo 
- Handle the file name from the specialization driver to the mpi task 
- Add non power of two handler in the branch sharing
- See the code structuration 
