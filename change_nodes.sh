#!/bin/bash

n=$1

sed -i "s/NODES_job/$n/g" qs_vasp
