IDIR =include
VPATH=src/:src/ext
ODIR=src/obj
LDIR =lib
BINDIR=bin
MYBIN=/home/maan/pulsar_softwares/bin/

CC=gcc
CFLAGS=-I$(IDIR) -Wno-unused-result -O3 -march=native
LIBS=-lm -lfftw3 -lcpgplot

_DEPS = header.h  rficlean.h  version.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = bcast_header.o  pack_unpack.o  rficlean.o  scaledata.o  strings_equal.o  cleanit.o  read_block.o  rficlean_data.o  send_stuff.o swap_bytes.o  nsamples.o  read_header.o sizeof_file.o plot_data.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c  $(DEPS)
	@mkdir -p $(BINDIR)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BINDIR)/rficlean: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: install
install:$(BINDIR)/rficlean
	install $(BINDIR)/rficlean $(MYBIN)/

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core


