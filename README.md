# relectro

R package to analyze electrophysiological data. It is still in development and new functionalities are added regularly. The focus is on data from tetrode recording experiments. Most of the package is written in R but some sections are written in C for speed. Analysis can easily be run in parallel with the snow package. If you want to know more, a good starting point is the vignette included with the package (see the vignettes directory in the source code).

You will need to install the fftw3 C library to compute discrete Fourier transforms with relectro. You can install it on a linux machine.

* With Fedora: `dnf install -y fftw3-devel` 
* With Ubuntu: `apt-get install libfftw3-dev`

## Install

* Clone the relectro repository `git clone https://github.com/kevin-allen/relectro.git`

If you want to modify the code of relecto, I recommand using R studio together with the book "R Packages" by Hadley Wickham. Most of the tools used to develop relectro are presented in this book.

* From R studio, click File/Open project... and select relectro.Rproj

If you just want to use the functions and objects of relectro, go in the terminal and 
run the following
* R CMD build relecto
* R CMD INSTALL relectro