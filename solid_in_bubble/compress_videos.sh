#!/bin/bash

print_help() {
    cat << EOF
Compress and manipulate video files using H.264 codec.
For details, see https://trac.ffmpeg.org/wiki/Encode/H.264

Usage: ./$(basename "$0") [<options>] [<file> ...]
Options:
  -a, --accelerate=<factor> Accelerate playback in <factor> times.
  -c, --compatibility       Use the highest level of compatibility (~25% extra size).
  -e, --extension=<ext>     Output file extension (mp4 by default).
  -q, --quality=<crf>       Change CRF (Constant Rate Factor),
      where 0 is lossless, 23 is the default, and 51 is worst quality possible.
  -r, --resize=<factor>     Reduce the resolution by <factor> times.
  -s, --suffix=<word>       Output file name suffix (_compressed by default).
  -h, --help                Print this help.
EOF
}

### Coloring
declare -r RED='\033[1;31m'
declare -r WHITE='\033[1;97m'
declare -r NC='\033[0m'
echo $(which -a ffmpeg)
ffmpeg=/usr/local/bin/ffmpeg
### Default options
options=(
    -c:v libx264      # H.264 video codec
    -pix_fmt yuv420p  # safe pixel format
    -an               # remove audio
)
declare -a files
extension='mp4'
suffix='_compressed'
accel=1
resize=1

fatal_err() {
    local msg="$1"
    echo -e "${RED}Fatal Error:${WHITE} $msg${NC}" >&2
    exit 1
}

for arg; do case $arg in
    -a=*|--accelerate=*)    accel=${arg#*=};;
    -c|--compatibility)     compat=1;;
    -e=*|--extension=*)     extension=${arg#*=};;
    -q=*|--quality=*)       options+=("-crf ${arg#*=}");;
    -r=*|--resize=*)        resize=${arg#*=};;
    -s=*|--suffix=*)        suffix=${arg#*=};;
    -h|--help)              print_help; exit 0;;
    -*)                     fatal_err "Unknown option '$arg'.";;
    *)                      files+=("$arg");;
esac; done

[[ ${#files} -gt 0 ]] || fatal_err "No input file specified."

if [[ $compat ]]; then
    options+=(-profile:v baseline -level 3.0)
else
    options+=(-profile:v high -level 4.0)
fi

for file in "${files[@]}"; do
    out=${file%.*}$suffix.$extension
    fps=$(${ffmpeg} -i "$file" 2>&1 | grep -oE "[0-9]+ fps" | sed 's/ fps//')
    local_options=(
        -vf "fps=$fps*$accel,setpts=PTS/$accel,scale=iw/$resize:ih/$resize"
    )
    set -x
    ${ffmpeg} -i "$file" "${options[@]}" "${local_options[@]}" "$out"
    set +x
done
