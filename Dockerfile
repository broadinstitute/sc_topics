FROM jupyter/scipy-notebook


VOLUME /home/jovyan/work

## Add files to container
ADD . /jupyter-scanpy

## Switch to root to install packages
USER root

## Install Python Requirements
RUN pip3 install -r /jupyter-scanpy/requirements.txt

## Notebook user
USER $NB_USER

## Notebook working directory
WORKDIR /home/jovyan/work
