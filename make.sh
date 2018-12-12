#!/bin/bash

set -o nounset

CFITSIO_LIBS=$(pkg-config --libs cfitsio)
CFITSIO_CFLAGS=$(pkg-config --cflags cfitsio)
X11_LIBS=$(pkg-config --libs x11)
X11_CFLAGS=$(pkg-config --cflags x11)
PNG_LIBS=$(pkg-config --libs libpng)
PNG_CFLAGS=$(pkg-config --cflags libpng)
PGPLOT_LIBS="-L $PGPLOT_DIR -lcpgplot -lpgplot -lgfortran $X11_LIBS $PNG_LIBS"
PGPLOT_CFLAGS="-I$PGPLOT_DIR $X11_CFLAGS $PNG_CFLAGS"
FFTW_LIBS=$(pkg-config --libs fftw3)
FFTW_CFLAGS=$(pkg-config --cflags fftw3)

if [ $# -ne 1 ]
then
  echo usage: $0 install_dir
  exit
fi

INSTALL_DIR=$1


gcc -lm -o pfits_addFRB pfits_addFRB.c pfits_setup.c pfits_loader.c T2toolkit.c $CFITSIO_LIBS
install -m 755 pfits_addFRB $INSTALL_DIR

gcc -lm -o pfits_cal pfits_cal.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS $PGPLOT_CFLAGS $PGPLOT_LIBS
install -m 755 pfits_cal $INSTALL_DIR

gcc -lm -o pfits_change pfits_change.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS 
install -m 755 pfits_change $INSTALL_DIR

gcc -lm -o pfits_convertText pfits_convertText.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS
install -m 755 pfits_convertText $INSTALL_DIR

gcc -lm -o pfits_describe pfits_describe.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS 
install -m 755 pfits_describe $INSTALL_DIR

gcc -lm -o pfits_dedisperse pfits_dedisperse.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS $FFTW_CFLAGS $FFTW_LIBS
install -m755 pfits_dedisperse $INSTALL_DIR

gcc -lm -o pfits_fftSearch pfits_fftSearch.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS 
install -m 755 pfits_fftSearch $INSTALL_DIR

gcc -lm -o pfits_frb pfits_frb.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS $PGPLOT_CFLAGS $PGPLOT_LIBS
install -m 755 pfits_frb $INSTALL_DIR

gcc -lm -o pfits_fv pfits_fv.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS $PGPLOT_CFLAGS $PGPLOT_LIBS
install -m 755 pfits_fv $INSTALL_DIR

gcc -lm -o pfits_getZeroDM pfits_getZeroDM.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS
install -m 755 pfits_getZeroDM $INSTALL_DIR

gcc -lm -o pfits_merge pfits_merge.c $CFITSIO_CFLAGS $CFITSIO_LIBS
install -m 755 pfits_merge $INSTALL_DIR

gcc -lm -o pfits_output pfits_output.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS
install -m 755 pfits_output $INSTALL_DIR

gcc -lm -o pfits_plot pfits_plot.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS $PGPLOT_CFLAGS $PGPLOT_LIBS
install -m 755 pfits_plot $INSTALL_DIR

gcc -lm -o pfits_process pfits_process.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS
install -m 755 pfits_process $INSTALL_DIR

gcc -lm -o pfits_sim pfits_sim.c pfits_setup.c pfits_loader.c T2toolkit.c $CFITSIO_CFLAGS $CFITSIO_LIBS
install -m 755 pfits_sim $INSTALL_DIR

gcc -lm -o pfits_simRealData pfits_simRealData.c pfits_loader.c T2toolkit.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS 
install -m 755 pfits_simRealData $INSTALL_DIR

gcc -lm -o pfits_statistcs pfits_statistics.c pfits_setup.c pfits_loader.c $CFITSIO_CFLAGS $CFITSIO_LIBS $PGPLOT_CFLAGS $PGPLOT_LIBS
install -m 755 pfits_statistcs $INSTALL_DIR

gcc -lm -o pfits_zapProfile pfits_zapProfile.c pfits_setup.c $CFITSIO_CFLAGS $CFITSIO_LIBS $PGPLOT_CFLAGS $PGPLOT_LIBS
install -m 755 pfits_zapProfile $INSTALL_DIR
