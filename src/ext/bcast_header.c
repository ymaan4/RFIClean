#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "rficlean.h"
void bcast_header()
{
  /* go no further if header is not needed */
  if (headerless) return;
  /* broadcast the header parameters */
  send_string("HEADER_START");
  if (!strings_equal(source_name,"")) {
    send_string("source_name");
    send_string(source_name);
  }
  send_int("telescope_id",telescope_id); 
  send_int("machine_id",machine_id);
  send_coords(src_raj,src_dej,az_start,za_start);
  if (nchans==1) {
    send_int("data_type",2);
    send_double("refdm",refdm); 
  } else {
    send_int("data_type",1);
  }
  send_double("fch1",fch1);
  send_double("foff",foff);
  send_int("nchans",nchans);
  send_int("nbits",obits);
  send_double ("tstart",tstart); 
  send_double("tsamp",tsamp);
  send_int("nifs",nifs);
  send_int("barycentric",barycentric);
  send_int("nbeams",nbeams);
  send_int("ibeam",ibeam);
  send_string("HEADER_END");
}
