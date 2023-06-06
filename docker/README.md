# Docker build and push to ghcr.io
```bash
$ docker buildx build --load --platform linux/amd64 -t glrna-amd64:0.0.7 -f ./docker/dockerfile .
$ docker buildx build --load --platform linux/arm64 -t glrna-arm64:0.0.4 -f ./docker/dockerfile .
$ docker tag glrna-amd64:0.0.7 ghcr.io/deeprob/glrna-amd64:latest
$ docker tag glrna-arm64:0.0.4 ghcr.io/deeprob/glrna-arm64:latest
$ docker push ghcr.io/deeprob/glrna-amd64:latest
$ docker push ghcr.io/deeprob/glrna-arm64:latest
```

# Resources
1. https://stackoverflow.com/questions/54437030/how-can-i-create-a-docker-image-to-run-both-python-and-r
2. https://support.bioconductor.org/p/131866/
3. https://www.biostars.org/p/9525979/#9525981
4. https://www.docker.com/blog/multi-arch-images/

