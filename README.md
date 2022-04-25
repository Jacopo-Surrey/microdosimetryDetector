Internal repo. The code present here is mostly WIP. Relevant parts are periodically committed to the advanced example Radioprotection, published as part of each official Geant4 release.

To visualize the output you can run:

    python3 plot_fd.py

Afterwards, relevant microdosimetric quantities and sample RBE estimates can be obtained via:

    python3 get_RBE.py
    
### Current issues:
- due to some small changes necessary in main.cc, the program won't compile when using the latest version of Geant4. For the time being, use 10.5.x or 10.6.x instead
- multithreading with two-stage output has not been tested, and the eventID might not correspond properly
- only MicroDiamond, WaterPixel, and Telescope are available, the other detectors are commented out (minimal changes are needed for them to work)
- kinetic energy scoring currently only enabled for WaterPixel; it can be easily added to the other ones, if needed
- TPS-like kinetic energy scoring across the whole water pixel to be added soon
