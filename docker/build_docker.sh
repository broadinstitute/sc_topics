#!/bin/sh

set -ev

VERSION=`cat VERSION.txt`

docker build -t kdgosik/sc_topics:$VERSION .
docker build -t kdgosik/sc_topics .

