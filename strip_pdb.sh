#!/bin/bash

awk '{print $7 " " $8 " " $9}' bdna_rotated.pdb > bdna_rotated.dat
