# RFIClean
#### Mitigation of periodic and spiky RFI from filterbank data.

RFIClean was primarily designed to efficiently search and mitigate
periodic RFI from GMRT time-domain data. Over the time, RFIClean has evolved
to mitigate any spiky (in time or frequency) RFI as well, and from any SIGPROC
filterbank format data file. It is written primarily in C. For handling the
filterbank format I/O, RFIClean uses several modules from Duncan Lorimer's SIGPROC
(https://sourceforge.net/p/sigproc/wiki/Home/; thank you Dunc) which are included
here in the src/ext/ folder (some of these codes are also modified suitably).

#### Dependencies:
FFTW3
PGPLOT

* Once installed, use `rficlean -h` or `rficlean --help` to see the usage information.


### Bash-script based parallel processing
For faster processing, use the bash script in bin/ to use RFIClean on different parts of a single data file simultaneously and then combine the output products at the end --- parallel processing in a rather crude but efficient way.

### Diagnostic plot
RFIClean produces a diagnostic plot showing which Fourier frequencies are
mitigated from the data and how frequently. The plot data are also output
in a file.


Yogesh Maan  <ymaan4[@]gmail.com>
