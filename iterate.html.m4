<!doctype html public "-//IETF//DTD HTML 2.0//EN">
<HTML>
<HEAD>
<TITLE>
WWW G6 Bravais Lattice Determination
</TITLE> 
</HEAD> 
<BODY>  
<H2> G6 Bravais Lattice
Determination Interface </H2>
<P> by
<P> Lawrence C. Andrews, Thuridion, Inc.,
<A HREF=mailto:andrews@thuridion.com>andrews@thuridion.com</A> and
<BR>Herbert J. Bernstein, Bernstein+Sons,
<A HREF=mailto:yaya@bernstein-plus-sons.com>yaya@bernstein-plus-sons.com</A>
<FORM method=POST ACTION=FULLHTDOCS()>
<BR>
<P>
<STRONG>
Please read the NOTICE below before use of this web page
</STRONG>
<P>
<INPUT type="submit">
<INPUT type="reset">
<BR>
Output Style: 
<SELECT name="OutputStyle" size="2">
<option selected value="TEXT">text
<option value="CIF">CIF
</SELECT>
<H3>Select the crystal lattice centering:</H3> 
<SELECT name="Centering" size="4"> 
<option selected value="P">P (primitive)
<option value="A"> A (a-centered)
<option value="B"> B (b-centered)
<option value="C"> C (c-centered)
<option value="F"> F (all-faces-centered)
<option value="I"> I (body-centered)
<option value="R"> R (rhombohedral-obverse)
<option value="H"> H (hexagonal primitive)
</SELECT>
<H3>Specify the cell edge lengths and angles:</H3>
<BR>
<table>
<td>_cell.length_a <td><INPUT TYPE="text" NAME="A" VALUE="10." SIZE="9"> 
<td>_cell.angle_alpha <td> <INPUT TYPE="text" NAME="Alpha" VALUE="90." SIZE="9"><BR><tr>  
<td>_cell.length_b <td><INPUT TYPE="text" NAME="B" VALUE="10." SIZE="9"> 
<td>_cell.angle_beta  <td> <INPUT TYPE="text" NAME="Beta" VALUE="90." SIZE="9"><BR><tr>
<td>_cell.length_c <td><INPUT TYPE="text" NAME="C" VALUE="10." SIZE="9"> 
<td>_cell.angle_gamma <td> <INPUT TYPE="text" NAME="Gamma" VALUE="90." SIZE="9"><BR><tr>
</table>


<H3>Specify the cell edge length esd's and angle esd's:</H3>
<BR>
<table>
<td>_cell.length_a_esd <td> <INPUT TYPE="text" NAME="sigA" VALUE=".15" SIZE="9">   
<td>_cell.angle_alpha_esd  <td> <INPUT TYPE="text" NAME="sigAlpha" VALUE=".2" SIZE="9"><BR><tr>
<td>_cell.length_b_esd<td> <INPUT TYPE="text" NAME="sigB" VALUE=".15" SIZE="9"> 
<td>_cell.angle_beta_esd   <td> <INPUT TYPE="text" NAME="sigBeta" VALUE=".2" SIZE="9"><BR><tr> 
<td>_cell.length_c_esd<td> <INPUT TYPE="text" NAME="sigC" VALUE=".15" SIZE="9">   
<td>_cell.angle_gamma_esd  <td> <INPUT TYPE="text" NAME="sigGamma" VALUE=".2" SIZE="9"><BR><tr>
</table>
<hr>
<INPUT type="hidden" NAME="Flush" VALUE="DUMMY">
<INPUT type="submit">
<INPUT type="reset">
</Form> <hr>

<H2>NOTICE</H2>
<P>
<P>
<STRONG>
Some of the software and documents included within this
software package are the intellectual property of various
parties, and placement in this package does not in anyway
imply that any such rights have in any way been waived
or diminished.
<P>
With respect to any software or documents for which a
copyright exists, ALL RIGHTS ARE RESERVED TO THE OWNERS
OF SUCH COPYRIGHT.
<P>
Even though the authors of the various documents
and software found here have made a good faith
effort to ensure that the documents are correct and
that the software performs according to its
documentation, and we would greatly appreciate
hearing of any problems you
may encounter, the programs and documents
any files created by the programs are provided
**AS IS** without any warrantee as to
correctness, merchantability or fitness
for any particular or general use.
<P>
THE RESPONSIBILITY FOR ANY ADVERSE
CONSEQUENCES FROM THE USE OF PROGRAMS
OR DOCUMENTS OR ANY FILE OR FILES CREATED BY
USE OF THE PROGRAMS OR DOCUMENTS
LIES SOLELY WITH THE USERS OF THE PROGRAMS
OR DOCUMENTS
OR FILE OR FILES AND NOT WITH AUTHORS OF
THE PROGRAMS OR DOCUMENTS.                      
</STRONG>
<P>
<hr>
<H2> Access to the source of ITERATE </H2>
<P>
This program and related scripts are available as
<A HREF="iterate.shar"> a self-extracting shell-script archive</A> or as
<A HREF="iterate.cshar"> a self-extracting C-shell-script archive.</A>
<H2> What Does This Web Page Do? </H2>
<P>
In simple terms, what this page does is to find the cells which are
&quot;close&quot; to the cell given, in order to help find the
Bravais lattice of highest symmetry consistent with the cell.
<P>A central problem in the solution of every crystal structure
is to determine the correct Bravais lattice of the crystal.
The Bravais lattices as they are usually listed
are:<P>
<table>
<td>aP   <td>triclinic (anorthic) primitive <BR><tr>
<td>mP   <td>monoclinic primitive <BR><tr>
<td>mS   <td>monoclinic side-centered (usually C-centered) <BR><tr>
<td>oP   <td>orthorhombic primitive <BR><tr>
<td>oS   <td>orthorhombic side-centered <BR><tr>
<td>oF   <td>orthorhombic face-centered <BR><tr>
<td>oI   <td>orthorhombic body-centered <BR><tr>
<td>hP   <td>hexagonal primitive <BR><tr>
<td>hR   <td>hexagonal rhombohedrally-centered <BR><tr>
<td>tP   <td>tetragonal primitive <BR><tr>
<td>tI   <td>tetragonal body-centered <BR><tr>
<td>cP   <td>cubic primitive <BR><tr>
<td>cF   <td>cubic face-centered <BR><tr>
<td>cI   <td>cubic body-centered <BR><tr>
</table>

<P>Failure to find the highest correct symmetry has several consequences, the
worst of which is that the structure may not be solved. The least of the
consequences is that Richard Marsh may 
publish a paper that points out the error,
corrects it, and finds a better solution 
to the structure. Many methods have been
described for finding the correct Bravais lattice. A summary of the published
methods was published in the paper that 
described the G6 formalism (which is used
in the program on this web page).

<P>&quot;Lattices and Reduced Cells as 
Points in 6-Space and Selection of Bravais
Lattice Type by Projections.&quot; 
Lawrence C. Andrews and Herbert J. Bernstein, Acta
Crystallographica, A44, 1009-1018 (1988).

<P>The program on this Web page implements a 
search in G6 for the various Bravais
lattices that the user's cell may fit. For each lattice type, the best metric
match is reported. If the higher symmetry type is actually correct, then that is
likely to be the best cell from which to start further refinement. However, the
possibility exists that one of the rejected cells (which did not match as well)
was actually the correct one to use. The reason for this ambiguity is
experimental error and its propagation in the transformations of the lattices in
the program. Fortunately, the rejected cells are usually quite similar to the
accepted one.

<P>A note on standard deviations: First, even in the best of circumstances,
standard deviations of unit cell dimensions 
from 4-circle diffractometer data are
always underestimated (by at least a factor 
of 2). In addition, the points chosen
for the determination are often not well distributed (for example all in the
first octant of orthorhombic lattices). These less than optimal choices cause
substantial systematic error. The experimental errors are amplified in the
mathematical conversions between various 
lattices that any lattice search program
must perform.  It is not a rare occurrence for angles to be incorrect by 0.5
degrees in initial unit cell determinations.
<P> <STRONG>Note:</STRONG> Even in most well determined unit
cells, the actual errors in the edge lengths is 0.2 to 0.5 parts per thousand.
(Note that reproducibility of the measurements is substantially better, leading
to the illusion that diffractometers 
produce excellent unit cell parameters). Use
of standard deviations that are too small is a common reason for failure of
Bravais lattice searches. For small molecules, 0.1 Angstroms is a reasonable
error for the edge lengths, for proteins, 0.4 to 0.5 (or even more for
preliminary measurements). Accurate unit cell parameters must by determined by a
number of more complex methods and must `include' extrapolation to remove
systematic effects. For an excellent summary, 
see &quot;Xray Structure Determination&quot;,
G.H.Stout and L.H.Jensen, Wiley, 1989.
</BODY>
</HTML>
