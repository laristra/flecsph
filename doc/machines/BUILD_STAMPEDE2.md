## Build FleCSPH on Stampede2 Supercomputer

This document explains how to build FleCSPH on Stampede2 supercomputer. 
Stampede2 is the flagship supercomputer at The University of Texas at Austin's 
Texas Advanced Computing Center (TACC). 
Specification about stampede2 supercomputer can be found 
[here](https://portal.tacc.utexas.edu/user-guides/stampede2).

### Set-up in Stampede2
Stampede2 mounts three file Lustre file systems that are shared across all nodes: 
the home, work, and scratch file systems. Stampede2's startup mechanisms define 
corresponding account-level environment variables `$HOME`, `$SCRATCH`, and `$WORK` 
that store the paths to directories that you own on each of these file systems.
Check [Stampede2 file system](https://portal.tacc.utexas.edu/user-guides/stampede2#files)
for more detail information. This [table](https://portal.tacc.utexas.edu/user-guides/stampede2#table3)
shows the basic characteristics of these file systems.

We encourage to use scratch directory to build FleCSPH. To do this, type
```
cds
```
when you login to your account. Then follow the suggested directory structure 
which is decsribed 
[here](https://github.com/laristra/flecsph/blob/master/README.md#suggested-directory-structure)


### Module List

Below is working module
```
  1) intel/18.0.2      3) impi/18.0.2    5) python3/3.7.0   7) git/2.9.0
  2) libfabric/1.6.1   4) cmake/3.10.2   6) boost/1.68
```
Once you have all these module, remaining procedure is exactly same as described in 
[README.md](https://github.com/laristra/flecsph/blob/master/README.md)

Note that Stampede2 has a `gcc/7.3` module but many libs/modules depend on intel complier.
So, we encourage to use above modules to build FleCSPH. 

To see details about these modules in stampede2, type 
```
module spider <module name>
```

## Contact
If you have questions or find more problems, please contact Hyun Lim (hylim1988@gmail.com)

