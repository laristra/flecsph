#!/bin/bash

ffmpeg -f 24 -f image2 -s 800x800 -i output_%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4
