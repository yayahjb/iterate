#
#  cryst1-2-cif.awk
#
#  simple awk script to convert a file with PDB CRYST1 records
#  to a CIF for input to iterate.  Extracted from pdb2cif
#  by Bourne, Bernstein and Bernstein.  See
#     http://ndbserver.rutgers.edu/software/pdb2cif
#
#  Warning -- this script assumes the entry id is in columns 73-76
#  If this is not true, then you make get a blank entry_id
#
BEGIN {
  printf ("\n data_cells\n\n loop_\n")
  printf ("_cell.entry_id\n")
  printf ("_cell.length_a\n")
  printf ("_cell.length_b\n")
  printf ("_cell.length_c\n")
  printf ("_cell.angle_alpha\n")
  printf ("_cell.angle_beta\n")
  printf ("_cell.angle_gamma\n")
  printf ("_cell.volume\n")
  printf ("_cell.Z_PDB\n")
  printf ("_cell.space_group_name_H-M  # not defined in dictionary\n\n\n")
}
#==========================================================================
#  keyword CRYST1
#
#
{
  if ($1 == "CRYST1") {
  #
  #  Contains a b c alpha beta gamma SG Z
  #

  # calculate cell volume

  {
    ca = cos(substr( ($0),34, 7) * 0.0174532)
    cb = cos(substr( ($0),41, 7) * 0.0174532)
    cc = cos(substr( ($0),48, 7) * 0.0174532)
    cz = (1.0 - (ca*ca - cb*cb - cc*cc) + (2.0*ca*cb*cc))
    vol = (substr( ($0), 7, 9) *\
           substr( ($0),16, 9) *\
           substr( ($0),25, 9) * (sqrt(cz)))
  }
  # localize space group and Z

  {
    sg = substr( ($0), 56, 11)
    Z  = substr( ($0), 67, 4 )
  }
  if (vol-1 > .01) {
  printf (" %s",substr( ($0), 73, 4))
  printf (" %9.3f", substr( ($0), 7, 9))
  printf (" %9.3f", substr( ($0),16, 9))
  printf (" %9.3f", substr( ($0),25, 9))
  printf (" %7.2f", substr( ($0),34, 7))
  printf (" %7.2f", substr( ($0),41, 7))
  printf (" %7.2f\n", substr( ($0),48, 7))
  printf (" %10.1f", vol)
  printf (" %3d ", Z)
  printf (" '%11s'\n", sg)
  }
  }
}
END {}
