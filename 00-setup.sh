#!/bin/bash

cd ../../cEMM.machlearn.data

find ./ -maxdepth 2 \( -name *.HEAD -o -name *.BRIK* \) -print0 | xargs -0i cp -d --parents {} ../../cEMM.newpipeline/data/
