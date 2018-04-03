This is an all-in-one installation of HTSlib, SAMtools and BCFtools, which are
compiled with gcc-4.6.4 on CentOS5 largely using the following command lines:

  (cd htslib-1.0; make prefix=/path/to/install_dir install)
  (cd samtools-1.0; make HTSDIR=../htslib-1.0 prefix=/path/to/install_dir install)
  (cd bcftools-1.0; make HTSDIR=../htslib-1.0 prefix=/path/to/install_dir install)

Note that the SAMtools Makefile has been modified such that libcurses is
statically linked. Libz, libc, libm and libpthread are dynamically linked.
In addition, the lib/pkgconfig directory has been removed as it keeps the
absolute installation path which is system dependent.
