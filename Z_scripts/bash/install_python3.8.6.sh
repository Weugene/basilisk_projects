pyver=3.8.6
pyver_short=3.8
parver=5.9.0
paraview_bin_path=$HOME/opt/
wget https://www.python.org/ftp/python/${pyver}/Python-${pyver}.tgz
tar xzvf Python-${pyver}.tgz
cd Python-${pyver}
./configure --prefix=$HOME/.local
make -j8
make install

#cd $HOME/.local/lib/python${pyver_short}/lib-dynload
#wget "https://github.com/azwyane/simple_text_editor/raw/master/edi/_bz2.cpython-37m-x86_64-linux-gnu.so"

cd ${paraview_bin_path}/ParaView-${parver}-MPI-Linux-Python${pyver_short}-64bit/lib
mv python${pyver_short} python${pyver_short}_
ln -s ~/.local/lib/python${pyver_short} .
pip${pyver_short} install numpy pandas plotly orca sklearn networkx psutil requests statsmodels tsmoothie --user
