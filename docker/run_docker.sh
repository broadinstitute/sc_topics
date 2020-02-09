#!/bin/sh

docker run --rm -ti -v `pwd`:/home/jovyan/work -p 8888:8888 kdgosik/sc_topics
