# Docker build and push to shub
```bash
$ # build image
$ sudo singularity build glrna.sif deffile
$ # push image to singularity cloud library
$ sudo singularity push -U glrna.sif library://jdr6123/test/glrna:0.0.1
$ # Pull image from singularity cloud library
$ singularity pull --arch amd64 glrna.sif library://jdr6123/test/glrna:0.0.1
```

# Resources
1. https://sylabs.io/2023/03/installing-singularityce-on-macos-with-apple-silicon-using-utm-rocky/
