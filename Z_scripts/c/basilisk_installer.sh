#!/bin/bash
set -x

print_help() {
cat << EOF
Installing of Basilisk
For details, see http://basilisk.fr/src/INSTALL
please, write 'sudo -s', before ./basilisk_installer.sh
Usage: ./$(basename "$0") [<options>]
Options:
  CONFIG, choose compilator config file, default is config.gcc. 
  Note that on 32-bits systems you need to use config.gcc.32bits instead
  BASILISK_PATH, choose path where Basilisk will be installed, default is /opt
  WBASILISK, choose the path to the workspace, default is ~/wbasilisk
  IS_SUPERCOMPUTER, set 'T' if you use supercomputer (there is no graphics support), default is 'F'
P.S. Options can be written in any order
EOF
}

print_help

if [ "$(whoami)" != "root" ]; then
	echo "Sorry, you are not root."
	exit 1
fi

read -p "Press [Enter] key to start installing..."

#default values
BASILISK_PATH=/opt
CONFIG=config.gcc
WBASILISK=~/wbasilisk
LOCALPATH=$PWD
IS_SUPERCOMPUTER=F
#redefining of variables
for ARGUMENT in "$@"
do

    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=)   

    case "$KEY" in
            CONFIG)              CONFIG=${VALUE} ;;
            WBASILISK)           WBASILISK=${VALUE} ;;
	    BASILISK_PATH)       BASILISK_PATH=${VALUE} ;;   
	    IS_SUPERCOMPUTER)    IS_SUPERCOMPUTER=${VALUE} ;;  
            *)   
    esac    


done

#installing of git analog, the fast lexical analyzer generator and make
sudo apt-get install darcs flex make
#fetching of source code from the site and saving it in /opt
cd $BASILISK_PATH

#if installing is stopped, then the 'bailisk' folder has already existed
darcs get --lazy http://basilisk.fr/basilisk

cd basilisk
darcs pull
make clean #?
make       #?
cd src
#exporting Basilisk environments
echo "export BASILISK=$PWD" >> ~/.bashrc
echo "export PATH=\$PATH:$PWD" >> ~/.bashrc
export BASILISK=$PWD
export PATH=$PATH:$BASILISK

sudo ln -sf $CONFIG config
if [ $IS_SUPERCOMPUTER == 'T' ]; then
	echo "export OPENGLIBS='-lfb_osmesa -lGLU -lOSMesa'" 
	echo 'OPENGLIBS=-lfb_osmesa -lGLU -lOSMesa' 
else 
	echo "export OPENGLIBS='-lfb_glx -lGLU -lGLEW -lGL -lX11'"  
	echo 'OPENGLIBS=-lfb_glx -lGLU -lGLEW -lGL -lX11' 
fi

make -k
make

#additional packages

sudo apt-get install gnuplot imagemagick ffmpeg  smpeg-plaympeg graphviz valgrind gifsicle pstoedit
#sudo apt install libglu1-mesa-dev freeglut3-dev
sudo apt-get install libglu1-mesa-dev libglew-dev libgl1-mesa-dev
sudo apt-get install libglu1-mesa-dev libosmesa6-dev
#Visualization
#Online visualisation with Basilisk View requires a separate installation.
#Offline (interactive) visualisation requires the installation of the Basilisk View servers.
#A short recipe for installation of both the interactive and non-interactive versions of Basilisk View is:
cd $BASILISK/gl
make libglutils.a libfb_glx.a

#BASILISK View Server
cd $BASILISK
make bview-servers

#Using Basilisk with python
sudo apt-get install swig libpython-dev
echo "creating a workspace for Basilisk in $WBASILISK"

mkdir $WBASILISK
echo "Testings"
cd $WBASILISK
mkdir bump
cd bump
cp $LOCALPATH/bump_test.c .
$BASILISK/qcc -O2 bump_test.c -o bump_test -lm
./bump_test > out.ppm 2> log
animate out.ppm
set +x

