#!/bin/sh

file=nodes.txt

out=.txt

for i in {30..60..10}

        do

                x=$i${file}

                y=$i${out}

                ./pomcp --problem pocman --outputfile $y --runs 1 --mindoubles 1 --maxdoubles 2 --accuracy 0.1 --horizon 3 --verbose 0 --autoexploration false --timeout 3600 --userave false --ravediscount 1.0 --raveconstant 0.01 --disabletree false --exploration 1 --usetransforms true

       done
