#!/bin/bash
#--------------------------------
# Transform the eps to png.
#--------------------------------
for eps in *.eps
do 
    png=${eps//".eps"/".png"}
    convert -density 300 -background white -flatten "$eps" "$png"
    echo "$eps to $png is done!"
done