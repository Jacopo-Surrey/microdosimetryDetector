WIP

To visualize the output you can run:

    python3 plot_fd.py
    
### Current issues:
- only MicroDiamond is available, the other detectors are commented out (minimal changes are needed for them to work)
- changing cuts via macro doesn't seem to have any effect (this is being looked into)
- is the size of the substrate and the region being set depending on the space occupied by the SV? They might end up outside of it!
- the applyChangesToGeometryCmd recreates geometry and physics regardless of whether this is actually needed (to-do: add a bool to keep track of this)
