#!/bin/sh

set -ev

VERSION=`cat VERSION.txt`

docker push kdgosik/sc_topics:${VERSION}
docker push kdgosik/sc_topics:latest
