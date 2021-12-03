Internal repo. The code present here is mostly WIP. Relevant parts are periodically committed to the advanced example Radioprotection, published as part of each official Geant4 release.

To visualize the output you can run:

    python3 plot_fd.py

Afterwards, relevant microdosimetric quantities and sample RBE estimates can be obtained via:

    python3 get_RBE.py
    
### Current issues:
- multithreading with two-stage output has not been tested, and the eventID might not correspond properly
- only MicroDiamond and Telescope are available, the other detectors are commented out (minimal changes are needed for them to work)
- changing cuts via macro doesn't seem to have any effect (this is being looked into)
- is the size of the substrate and the region being set depending on the space occupied by the SV? They might end up outside of it!
