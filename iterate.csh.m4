#!/bin/csh
# iterate.csh
#
# Herbert J. Bernstein, Bernstein + Sons
# Lawrence C. Andrews, Thuridion, Inc.
#
# 29 September 1996
#
# This is a service script for the iterate.html web page
# It must be placed in an appropriate cgi-bin directory on
# the server pointed to by iterate.html
#
#
# To operate correctly, the programs tr and sed must be in the
# default path and the /bin/echo version of echo must follow
# system V conventions sufficiently to produce an empty line
# call, below
#
/bin/echo "Content-type: text/html"
/bin/echo 
echo "<HEAD>"
echo "<TITLE>G6 Lattice Identification"
echo "</TITLE>"
echo "</HEAD>"
echo "<BODY>"
tr '\&' '\n'  |sed "s/^./set &/"  > /tmp/outstr$$
#cat /tmp/outstr$$
source /tmp/outstr$$
rm /tmp/outstr$$
echo "<H3># G6 Lattice Identification</H3>"
echo $Centering > /tmp/instr$$
echo "<P>#  Centering: " $Centering 
echo $A $B $C $Alpha $Beta $Gamma >>/tmp/instr$$
echo "<P># Cell: " $A $B $C $Alpha $Beta $Gamma
echo $sigA $sigB $sigC $sigAlpha $sigBeta $sigGamma >>/tmp/instr$$
echo "<P># Sigmas: " $sigA $sigB $sigC $sigAlpha $sigBeta $sigGamma
echo "<P><H3># Results of ITERATE Run</H3>"
setenv ITERATE_QUERY NO
setenv OUTPUT_STYLE $OutputStyle
echo "<PRE>"
BINPATH() < /tmp/instr$$
rm /tmp/instr$$
#cat /tmp/instr$$ 
echo "</PRE>"
echo "</BODY>"
