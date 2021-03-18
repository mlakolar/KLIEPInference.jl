#!/bin/bash

parallel --delay .2 julia experiment3.jl ${1} ${2} {1} ::: {1..1000}