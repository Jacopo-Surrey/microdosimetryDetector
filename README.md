## how to compile
Geant4 simulation of a microdosimetry detector for hadrontherapy

Set the number of particles in the run.mac macro, then move to the build directory and execute:

    cmake /path/to/source
    make -j4
    ./myProgram run.mac

To visualize the output you can then run:

    python3 plot_fd.py
    
### ToDo
- [ ] properly write this readme
- [ ] finish diamond microdosimeter
- [ ] add silicon microdosimeter (cfr article)
- [ ] switch between the two with and #ifdef
- [ ] add a detector-agnostic region with higher precision
- [ ] set the source position via c++
- [ ] improve the GPS macro
- [ ] confront spectra with published ones
- [ ] confront spectra with those obtained by Gabriele
