# Docker build and push to ghcr.io
```bash
$ docker image build -t glrna:0.0.4 -f ./docker/dockerfile .
$ docker tag glrna:0.0.4 ghcr.io/deeprob/glrna:latest
$ docker push ghcr.io/deeprob/glrna
```

# Resources
1. https://stackoverflow.com/questions/54437030/how-can-i-create-a-docker-image-to-run-both-python-and-r
2. https://support.bioconductor.org/p/131866/
3. https://www.biostars.org/p/9525979/#9525981

