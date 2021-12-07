Internal repo. The code present here is mostly WIP. Relevant parts are periodically committed to the advanced example Radioprotection, published as part of each official Geant4 release.

To visualize the output you can run:

    python3 plot_fd.py

Afterwards, relevant microdosimetric quantities and sample RBE estimates can be obtained via:

    python3 get_RBE.py

### To-do:
- currently SV are set far apart, so that it's unlikely that an event should be registered in more than one at the same time. This means that if the detector is very small, most of the space is empty, and the simulation takes much longer. If hits on different detectors are scored as different events this shouldn't be a problem, and the detectors could be placed much closer together (albeit some data might not be fully independend).
- to achieve the previous point: check how much time I lose assigning multiple SD (one per SV) rather then one for the whole logical volume
    
### Current issues:
- multithreading with two-stage output has not been tested, and the eventID might not correspond properly
- only MicroDiamond and Telescope are available, the other detectors are commented out (minimal changes are needed for them to work)
- changing cuts via macro doesn't seem to have any effect (this is being looked into)
