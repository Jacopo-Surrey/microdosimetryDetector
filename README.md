WIP

To visualize the output you can run:

    python3 plot_fd.py
    
### Current issues:
- running the simulations without water phantom causes a segfault, since it tries to set cut for a non-existent region
- only MicroDiamond is available, the other detectors are commented out (minimal changes are needed for them to work)
- changing cuts via macro doesn't seem to have any effect (this is being looked into)
- is the size of the substrate and the region being set depending on the space occupied by the SV? They might end up outside of it!
