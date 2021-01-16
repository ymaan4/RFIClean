# RFIClean
### Mitigation of periodic and spiky RFI from filterbank data.

RFIClean excises periodic RFI (broadband as well as narrow-band) in
the Fourier domain, and then mitigates narrow-band spectral line RFI
as well as broadband bursty time-domain RFI using robust statistics.
Currently, RFIClean anticipates the input data either in SIGPROC
filterbank format or the GMRT's pulsar data format.

RFIClean was primarily designed to efficiently search and mitigate
periodic RFI from GMRT time-domain data. Over the time, RFIClean has evolved
to mitigate any spiky (in time or frequency) RFI as well, and from any SIGPROC
filterbank format data file. It is written primarily in C. For handling the
filterbank format I/O, RFIClean uses several modules from Duncan Lorimer's SIGPROC
(https://sourceforge.net/p/sigproc/wiki/Home/; thank you Dunc) which are included
here in the src/ext/ folder (some of these codes are also modified suitably).

### Installation:
* For compiling RFIClean, a `Makefile` is included in the package.
* If you want to install the executable at a location other than `RFIClean/bin/`, then change `MYBIN` in the `Makefile` accordingly.
* RFIClean has the following dependencies: `FFTW3` and `PGPLOT`. If these are not included in the regular library path then amend `LIBS` in the `Makefile` accordingly.
* To compile, run `make`. For installing the executable in your favourite location, run `make install`.
* Once installed, use `rficlean -h` or `rficlean --help` to see the usage information.


### Bash-script based parallel processing
For faster processing, use the bash script in bin/ to use RFIClean on different parts of a single data file simultaneously and then combine the output products at the end --- parallel processing in a rather crude but very efficient way.

### Diagnostic plot
RFIClean produces a diagnostic plot showing which Fourier frequencies are
mitigated from the data and how frequently. The plot data are also output
in a file.

### RFIClean paper
Details of the methods used in RFIClean as well as some of the early
scientific contributions of RFIClean can be found [in this paper.](https://arxiv.org/abs/2012.11630)

#
If you find RFIClean useful, it would be great if you could cite the following paper that describes the underlying method in detail: https://arxiv.org/abs/2012.11630 . Thanks!

Yogesh Maan  <ymaan4[@]gmail.com>
