#!/bin/bash

# for image in *.png; do convert $image -resize 50% -quality 100 half_size/$image; done
image=$1
convert $image -resize 50% -quality 100 half_size/$image