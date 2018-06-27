#!/bin/bash

# TODO: add description and usage help
#       
ffmpeg -framerate 60 -f image2 -s 2048x2048 -i $1 -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4
