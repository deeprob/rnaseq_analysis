# Docker build and push to shub
```bash
$ sudo singularity build glrna.sif deffile
$ sudo singularity push -U glrna.sif library://USER/REPO/IMG

# Pull image from shub
```bash
$ singularity pull glrna.sif library://jdr6123/test/glrna:TAG
```

# Resources
1. https://sylabs.io/2023/03/installing-singularityce-on-macos-with-apple-silicon-using-utm-rocky/
2. 
