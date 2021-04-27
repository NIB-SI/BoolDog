# Dockerised SQUAD

Modified SQUAD from https://www.vital-it.ch/research/software/SQUAD that works with docker. 

## Usage

To build:

```bash
cd squad
docker build --tag squad .
```

## Running on \*nix

(only tested on Ubuntu)

```bash
xhost +"local:docker@"
docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -v $PWD/squad2-2/samples:/root/samples squad
```

## Running on Windows

Install an X server (e.g. vcxsrv).  

```bash
docker run -e DISPLAY=host.docker.internal:0 -v $PWD/squad2-2/samples:/root/samples  squad
```


## Issues

### Graphical Interface
* For an overview of running a GUI with Docker see: [Running Desktop Apps in Docker
](https://betterprogramming.pub/running-desktop-apps-in-docker-43a70a5265c4)
* (Windows) TLDR: To see if your X server is running/configured, a test is to run xeyes:
`docker run --rm -ti -e DISPLAY=host.docker.internal:0 fr3nd/xeyes`

### SQUAD java errors
These are errors that occur when attempting to run SQUAD. So far any Boolean analysis throws errors. Without access to the source code, debugging is nigh impossible. Attempts have been made to resolve missing libraries, but at this point it is unlikely that the current errors will be resolved. 

* (Windows) After loading a network and clicking "Run Analysis":  
`java.lang.UnsatisfiedLinkError: no st in java.library.path` 
* (Ubuntu) After loading a network and clicking "Run Analysis", a `MMalloc` error occurs