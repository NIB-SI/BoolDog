# Dockerised SQUAD

Modified SQUAD from https://www.vital-it.ch/research/software/SQUAD that works with docker. 

To build:

```bash
docker build --tag squad .
```

To use (only tested on Ubuntu, display issues may occur):

```bash
xhost +"local:docker@"
docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -v $PWD/squad2-2/samples:/root/samples squad
```
