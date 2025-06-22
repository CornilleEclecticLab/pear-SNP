#!/usr/bin/env bash

# singularity pull MEGAnE.NOGIT.sif docker://shoheikojima/megane:v1.0.1.beta

# 2024-04-12
singularity build --fakeroot MEGAnE.NOGIT.sif docker://shoheikojima/megane:v1.0.1.beta