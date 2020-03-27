## how to compile
Geant4 simulation of a microdosimetry detector for hadrontherapy

Set the number of particles in the `run.mac` macro (optionally the type of source in `gps.mac`). Then move to the build directory and execute:

    cmake /path/to/source

to compile with the diamond detector. Alternatively the reference silicon detector can be used with:

    cmake /path/to/source -DUSE_SILICON=ON

The program can then be compiled and executed with

    make -j4
    ./myProgram run.mac

To visualize the output you can then run:

    python3 plot_fd.py
    
### ToDo
- [ ] properly write this readme
- [ ] finish diamond microdosimeter
- [ ] add silicon microdosimeter (cfr article)
- [x] switch between the two with and `#ifdef`
- [ ] add a detector-agnostic region with higher precision
- [ ] set the source position via c++
- [ ] improve the GPS macro
- [ ] confront spectra with published ones
- [ ] confront spectra with those obtained by Gabriele
- [ ] remove unused histograms (warnings)
- [ ] check why `delete analysis;` in `main.cc` causes the following error when executing a program from the command line:

    free(): invalid next size (fast)
    Aborted (core dumped)

- [ ] track type of hadron (its Z?) via the tracking action

