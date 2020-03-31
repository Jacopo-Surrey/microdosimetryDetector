### how to compile
Geant4 simulation of a microdosimetry detector for hadrontherapy

Move to the build directory and execute:

    cmake /path/to/source

to compile with the diamond detector. Alternatively the reference silicon detector can be used with (unimplemented as of yet):

    cmake /path/to/source -DUSE_SILICON=ON

If needed, set the number of particles in the `run.mac` macro, and/or the type of source in `gps.mac`, in the build directory

The program can then be compiled and executed with

    make -j4
    ./myProgram run.mac

To visualize the output you can then run:

    python3 plot_fd.py
    
### ToDo
- [ ] check that the cuts in the inner region have been implemented correctly
- [ ] properly write this readme
- [ ] finish diamond microdosimeter
- [ ] set the source position via c++ -- might be pointlessly complex
- [ ] improve the GPS macro
- [ ] confront spectra with published ones
- [ ] confront spectra with those obtained by Gabriele
- [ ] remove unused histograms (warnings)
- [ ] check why `delete analysis;` in `main.cc` causes the following error when executing a program from the command line: `free(): invalid next size (fast)`   `Aborted (core dumped)`
- [ ] track type of hadron (its Z?) via the tracking action
- [ ] add a rotation matrix to the diamond detector placement, so it can be rotated changing a single variable (as it can be rotated in the lab)
- [ ] modify the sensitive volume in the stepping action to always pick the active one
- [x] ~~switch between the two with and `#ifdef`~~ DONE
- [x] ~~remove the if condition inside the for loop when placing the diamond SV (change how the position is calculated)~~ DONE
- [x] ~~add silicon microdosimeter (cfr article)~~ DONE
- [x] ~~add a detector-agnostic region with higher precision~~ DONE
