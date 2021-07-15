#!/bin/bash
#://trac.ffmpeg.org/wiki/Encode/H.264

[[ $# < 2 ]] && { echo "Usage: ./$(basename $0) <input> <output> [<crf>]"; exit 1; }

# To vary size and quality, change CRF (Constant Rate Factor), where 0 is lossless, 23 is the default, and 51 is worst quality possible.
crf=25
[[ $# > 3 ]] && cfr=$3

ffmpeg -i $1 -c:v libx264 -profile:v high -level:v 4.0 -pix_fmt yuv420p -crf $crf $2 >/tmp/video$1.log
