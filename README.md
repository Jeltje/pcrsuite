# pcrsuite
PCRSuite cgi and html

Prerequisites: primer3_core installed (http://primer3.sourceforge.net/)

Note: PCRSuite was last tested with primer3_core release 1.1.4, later versions are not guaranteed to work.

The following lines in the PCRSuite.cgi script should be changed to something local.

my$file="/usr/local/apache/share/eur/htdocs/fgg/kgen/primer/stats.txt";
  use constant UPLOAD_DIR => "/usr/local/apache/share/eur/htdocs/fgg/kgen/primer/output/";
  my$PRIMER_BIN = '/srv/mblab/bin/pcr_pipeline/bin/primer3_core';
