# rfiClean
Mitigation of periodic and spiky RFI from filterbank data.

rfiClean was primarily designed to efficiently search and mitigate
periodic RFI from GMRT time-domain data. Over the time, rfiClean has evolved
to mitigate any spiky (in time or frequency) RFI as well, and from any SIGPROC
filterbank format data file. It is written primarily in C. For handling the
filterbank format I/O, rfiClean uses several modules from Duncan Lorimer's SIGPROC
(https://sourceforge.net/p/sigproc/wiki/Home/ ; thanks you Dunc) which are included
here in the src/ext/ folder.

# Dependencies:
FFTW3

PGPLOT



