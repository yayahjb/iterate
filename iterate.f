C     ITERATE -- Program for G6 Bravais Lattice Determination
C
C     by
C
C     Lawrence C. Andrews, Thuridion, Inc.,
C     andrews@thuridion.com
C 
C     and
C
C     Herbert J. Bernstein, Bernstein+Sons,
C     yaya@bernstein-plus-sons.com
C
C     This program finds the cells which are "close" to the cell given, 
C     in order to help find the Bravais lattice of highest symmetry 
C     consistent with the cell.
C     
C     A central problem in the solution of every crystal structure
C     is to determine the correct Bravais lattice of the crystal.
C     The Bravais lattices as they are usually listed are:
C
C     aP   triclinic (anorthic) primitive 
C     mP   monoclinic primitive 
C     mS   monoclinic side-centered (usually C-centered)
C     oP   orthorhombic primitive 
C     oS   orthorhombic side-centered 
C     oF   orthorhombic face-centered 
C     oI   orthorhombic body-centered 
C     hP   hexagonal primitive 
C     hR   hexagonal rhombohedrally-centered 
C     tP   tetragonal primitive 
C     tI   tetragonal body-centered 
C     cP   cubic primitive 
C     cF   cubic face-centered 
C     cI   cubic body-centered 
C     
C     
C     Failure to find the highest correct symmetry has several consequences, 
C     the worst of which is that the structure may not be solved. The 
C     least of the consequences is that Richard Marsh may publish a paper 
C     that points out the error, corrects it, and finds a better solution 
C     to the structure. Many methods have been described for finding the 
C     correct Bravais lattice. A summary of the published methods was 
C     published in the paper that described the G6 formalism (which is used
C     in this program).
C     
C     "Lattices and Reduced Cells as Points in 6-Space and Selection of 
C     Bravais Lattice Type by Projections", Lawrence C. Andrews and 
C     Herbert J. Bernstein, Acta Crystallographica, A44, 1009-1018 (1988).
C
C     This program accepts cell parameters and esd's and produces a list
C     of cells "close" to the cell given
C
C**********************************************************************C
      SUBROUTINE BADCAL (A,B)
      CHARACTER *6 A,B
C----------------------------------------------------------------------C
      include 'ITERATE.cmn'
      WRITE (*,*)  ' '//hm//' BAD SUBROUTINE CALL TO ',B
      WRITE (*,*)  ' '//hm//' CALLING NAME =',A
      STOP
      END


CC***********************************************************************
C      SUBROUTINE CPYVN (N,V1,V2)
C      REAL V1(N),V2(N)
C
C-----------------------------------------------------------------------
C      DO 1000 I=1,N
C 1000 V2(I) = V1(I)
C      END

C***********************************************************************
      subroutine bldprj (maxprj,nproj,itdesg,chrlat,pjn,prj,test)

C This function builds the projectors of Paciorek and Bonin, J. Appl.
C Cryst., 25, (5) pp 632-637. Internal checks of the correctness of
C the projectors are made before it exits.
C
C The output values are just transferred from the stored parameters
C in data statements. Note that if ngtype is negative, then the
C projector is not output. This is because those are ones that
C are easy to find later by searching. Either they are cases where
C two projectors are exactly the same subspace (and therefore the
C same projector) or else they are just simple exchange of axes
C (which are dealt with in mkrefl in program iterate). If the routine
C is extracted for other uses, then it may be correct to set the
C values all positive.

C   maxprj sets the maximum number of projectors to build
C   nproj  is the actual number that bldprj constructs
C   itdesg is the numeric designation of the Niggli type in the
C          International Tables for Crystallography
C   chrlat is the returned 2 character designators for the lattice
C          type (also call Pearson symbols)
C   pjn    is the normalizer of the projector -- the integer values
C          stored in the returned matrix need to be divided by the
C          normalizer to make the actual projector
C   prj    the integer part of the projector -- divide by the normalizer
C          (pjn) to get the actual projector


      include 'ITERATE.cmn'
      real pjn(maxprj)
      real prj(36,MAXPRJ)
      integer itdesg(maxprj)
      character *2 chrlat(maxprj)
      integer ngtype(42)
      real projct(36,42)
      real zprj(6,6,42)
      equivalence (projct,zprj)
      real pjnorm(42)
      character *2 lattyp(42)
      character *6 test

      data ngtype(1) /3/
      data lattyp(1) /'cP'/
      data pjnorm(1) /3./
      data( projct(i,1),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(2) /5/
      data lattyp(2) /'cI'/
      data pjnorm(2) /39./
      data (projct(i,2),i=1,36) /
     1   9, 9, 9,-6,-6,-6,
     2   9, 9, 9,-6,-6,-6,
     3   9, 9, 9,-6,-6,-6,
     4  -6,-6,-6, 4, 4, 4,
     5  -6,-6,-6, 4, 4, 4,
     6  -6,-6,-6, 4, 4, 4 /


      data ngtype(3) /1/
      data lattyp(3) /'cF'/
      data pjnorm(3) /6./
      data (projct(i,3),i=1,36) /
     1  1,1,1,1,1,1,
     2  1,1,1,1,1,1,
     3  1,1,1,1,1,1,
     4  1,1,1,1,1,1,
     5  1,1,1,1,1,1,
     6  1,1,1,1,1,1 /



      data ngtype(4) /11/
      data lattyp(4) /'tP'/
      data pjnorm(4) /2./
      data (projct(i,4),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(5) /-21/
      data lattyp(5) /'tP'/
      data pjnorm(5) /2./
      data (projct(i,5),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(6) /15/
      data lattyp(6) /'tI'/
      data pjnorm(6) /4./
      data (projct(i,6),i=1,36) /
     1  1, 1, 0,-1,-1, 0,
     2  1, 1, 0,-1,-1, 0,
     3  0, 0, 4, 0, 0, 0,
     4 -1,-1, 0, 1, 1, 0,
     5 -1,-1, 0, 1, 1, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(7) /6/
      data lattyp(7) /'tI'/
      data pjnorm(7) /26./
      data (projct(i,7),i=1,36) /
     1  6, 6, 6,-4,-4,-4,
     2  6, 6, 6,-4,-4,-4,
     3  6, 6, 6,-4,-4,-4,
     4 -4,-4,-4, 7, 7,-6,
     5 -4,-4,-4, 7, 7,-6,
     6 -4,-4,-4,-6,-6,20 /


      data ngtype(8) /-7/
      data lattyp(8) /'tI'/
      data pjnorm(8) /26./
      data (projct(i,8),i=1,36) /
     1  6, 6, 6,-4,-4,-4,
     2  6, 6, 6,-4,-4,-4,
     3  6, 6, 6,-4,-4,-4,
     4 -4,-4,-4,20,-6,-6,
     5 -4,-4,-4,-6, 7, 7,
     6 -4,-4,-4,-6, 7, 7 /


      data ngtype(9) /18/
      data lattyp(9) /'tI'/
      data pjnorm(9) /26./
      data (projct(i,9),i=1,36) /
     1  8, 0, 0, 4, 8, 8,
     2  0,13,13, 0, 0, 0,
     3  0,13,13, 0, 0, 0,
     4  4, 0, 0, 2, 4, 4,
     5  8, 0, 0, 4, 8, 8,
     6  8, 0, 0, 4, 8, 8 /


      data ngtype(10) /12/
      data lattyp(10) /'hP'/
      data pjnorm(10) /3./
      data (projct(i,10),i=1,36) /
     1  1, 1, 0, 0, 0,-1,
     2  1, 1, 0, 0, 0,-1,
     3  0, 0, 3, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1,-1, 0, 0, 0, 1 /


      data ngtype(11) /-22/
      data lattyp(11) /'hP'/
      data pjnorm(11) /3./
      data (projct(i,11),i=1,36) /
     1  3, 0, 0, 0, 0, 0,
     2  0, 1, 1,-1, 0, 0,
     3  0, 1, 1,-1, 0, 0,
     4  0,-1,-1, 1, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(12) /9/
      data lattyp(12) /'hR'/
      data pjnorm(12) /5./
      data (projct(i,12),i=1,36) /
     1  1,1,0,1,1,1,
     2  1,1,0,1,1,1,
     3  0,0,5,0,0,0,
     4  1,1,0,1,1,1,
     5  1,1,0,1,1,1,
     6  1,1,0,1,1,1  /


      data ngtype(13) /2/
      data lattyp(13) /'hR'/
      data pjnorm(13) /3./
      data (projct(i,13),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,1,1,1,
     5  0,0,0,1,1,1,
     6  0,0,0,1,1,1 /


      data ngtype(14) /-4/
      data lattyp(14) /'hR'/
      data pjnorm(14) /3./
      data (projct(i,14),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,1,1,1,
     5  0,0,0,1,1,1,
     6  0,0,0,1,1,1 /


      data ngtype(15) /24/
      data lattyp(15) /'hR'/
      data pjnorm(15) /53./
      data (projct(i,15),i=1,36) /
     1  27,  3,  3,  6,-18,-18,
     2   3, 18, 18,-17, -2, -2,
     3   3, 18, 18,-17, -2, -2,
     4   6,-17,-17, 19, -4, -4,
     5 -18, -2, -2, -4, 12, 12,
     6 -18, -2, -2, -4, 12, 12 /


      data ngtype(16) /32/
      data lattyp(16) /'oP'/
      data pjnorm(16) /1./
      data (projct(i,16),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(17) /36/
      data lattyp(17) /'oS'/
      data pjnorm(17) /2./
      data (projct(i,17),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(18) /-38/
      data lattyp(18) /'oS'/
      data pjnorm(18) /2./
      data (projct(i,18),i=1,36) /
     1  1, 0, 0, 0, 0,-1,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1, 0, 0, 0, 0, 1 /


      data ngtype(19) /13/
      data lattyp(19) /'oS'/
      data pjnorm(19) /2./
      data (projct(i,19),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,2 /


      data ngtype(20) /-23/
      data lattyp(20) /'oS'/
      data pjnorm(20) /2./
      data (projct(i,20),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(21) /-40/
      data lattyp(21) /'oS'/
      data pjnorm(21) /2./
      data (projct(i,21),i=1,36) /
     1  2, 0, 0, 0, 0, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(22) /16/
      data lattyp(22) /'oF'/
      data pjnorm(22) /10./
      data (projct(i,22),i=1,36) /
     1  3, 3, 0,-2,-2,-2,
     2  3, 3, 0,-2,-2,-2,
     3  0, 0,10, 0, 0, 0,
     4 -2,-2, 0, 3, 3,-2,
     5 -2,-2, 0, 3, 3,-2,
     6 -2,-2, 0,-2,-2, 8  /


      data ngtype(23) /26/
      data lattyp(23) /'oF'/
      data pjnorm(23) /13./
      data (projct(i,23),i=1,36) /
     1  4, 0, 0, 2, 4, 4,
     2  0,13, 0, 0, 0, 0,
     3  0, 0,13, 0, 0, 0,
     4  2, 0, 0, 1, 2, 2,
     5  4, 0, 0, 2, 4, 4,
     6  4, 0, 0, 2, 4, 4  /


      data ngtype(24) /8/
      data lattyp(24) /'oI'/
      data pjnorm(24) /13./
      data (projct(i,24),i=1,36) /
     1  3, 3, 3,-2,-2,-2,
     2  3, 3, 3,-2,-2,-2,
     3  3, 3, 3,-2,-2,-2,
     4 -2,-2,-2,10,-3,-3,
     5 -2,-2,-2,-3,10,-3,
     6 -2,-2,-2,-3,-3,10  /


      data ngtype(25) /19/
      data lattyp(25) /'oI'/
      data pjnorm(25) /6./
      data (projct(i,25),i=1,36) /
     1  2,0,0,0,2,2,
     2  0,3,3,0,0,0,
     3  0,3,3,0,0,0,
     4  0,0,0,6,0,0,
     5  2,0,0,0,2,2,
     6  2,0,0,0,2,2  /


      data ngtype(26) /42/
      data lattyp(26) /'oI'/
      data pjnorm(26) /2./
      data (projct(i,26),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(27) /33/
      data lattyp(27) /'mP'/
      data pjnorm(27) /1./
      data (projct(i,27),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,1,0,
     6  0,0,0,0,0,0  /


      data ngtype(28) /-35/
      data lattyp(28) /'mP'/
      data pjnorm(28) /1./
      data (projct(i,28),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,1,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0  /


      data ngtype(29) /-34/
      data lattyp(29) /'mP'/
      data pjnorm(29) /1./
      data (projct(i,29),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,1  /


      data ngtype(30) /39/
      data lattyp(30) /'mS'/
      data pjnorm(30) /2./
      data (projct(i,30),i=1,36) /
     1  1, 0, 0, 0, 0,-1,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 2, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1, 0, 0, 0, 0, 1  /


      data ngtype(31) /-41/
      data lattyp(31) /'mS'/
      data pjnorm(31) /2./
      data (projct(i,31),i=1,36) /
     1  2, 0, 0, 0, 0, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5  0, 0, 0, 0, 2, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(32) /-37/
      data lattyp(32) /'mS'/
      data pjnorm(32) /2./
      data (projct(i,32),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 2, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(33) /10/
      data lattyp(33) /'mS'/
      data pjnorm(33) /2./
      data (projct(i,33),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,1,1,0,
     5  0,0,0,1,1,0,
     6  0,0,0,0,0,2  /


      data ngtype(34) /-14/
      data lattyp(34) /'mS'/
      data pjnorm(34) /2./
      data (projct(i,34),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,1,1,0,
     5  0,0,0,1,1,0,
     6  0,0,0,0,0,2  /


      data ngtype(35) /-20/
      data lattyp(35) /'mS'/
      data pjnorm(35) /2./
      data (projct(i,35),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,1,1,
     6  0,0,0,0,1,1  /


      data ngtype(36) /-25/
      data lattyp(36) /'mS'/
      data pjnorm(36) /2./
      data (projct(i,36),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,1,1,
     6  0,0,0,0,1,1  /


      data ngtype(37) /28/
      data lattyp(37) /'mS'/
      data pjnorm(37) /10./
      data (projct(i,37),i=1,36) /
     1  5, 0, 0, 0, 5, 0,
     2  0,10, 0, 0, 0, 0,
     3  0, 0,10, 0, 0, 0,
     4  0, 0, 0, 2, 0, 4,
     5  5, 0, 0, 0, 5, 0,
     6  0, 0, 0, 4, 0, 8  /


      data ngtype(38) /-30/
      data lattyp(38) /'mS'/
      data pjnorm(38) /10./
      data (projct(i,38),i=1,36) /
     1  10, 0, 0, 0, 0, 0,
     2   0, 5, 0, 5, 0, 0,
     3   0, 0,10, 0, 0, 0,
     4   0, 5, 0, 5, 0, 0,
     5   0, 0, 0, 0, 2, 4,
     6   0, 0, 0, 0, 4, 8  /


      data ngtype(39) /-29/
      data lattyp(39) /'mS'/
      data pjnorm(39) /10./
      data (projct(i,39),i=1,36) /
     1  5, 0, 0, 0, 0, 5,
     2  0,10, 0, 0, 0, 0,
     3  0, 0,10, 0, 0, 0,
     4  0, 0, 0, 2, 4, 0,
     5  0, 0, 0, 4, 8, 0,
     6  5, 0, 0, 0, 0, 5  /


      data ngtype(40) /43/
      data lattyp(40) /'mI'/
      data pjnorm(40) /20./
      data (projct(i,40),i=1,36) /
     1  11, 1, 0, 1,-9,-4,
     2   1,11, 0,-9, 1,-4,
     3   0, 0,20, 0, 0, 0,
     4   1,-9, 0,11, 1,-4,
     5  -9, 1, 0, 1,11,-4,
     6  -4,-4, 0,-4,-4,16  /


      data ngtype(41) /17/
      data lattyp(41) /'mI'/
      data pjnorm(41) /10./
      data (projct(i,41),i=1,36) /
     1  3, 3, 0,-2,-2,-2,
     2  3, 3, 0,-2,-2,-2,
     3  0, 0,10, 0, 0, 0,
     4 -2,-2, 0, 8,-2,-2,
     5 -2,-2, 0,-2, 8,-2,
     6 -2,-2, 0,-2,-2, 8  /


      data ngtype(42) /27/
      data lattyp(42) /'mI'/
      data pjnorm(42) /3./
      data (projct(i,42),i=1,36) /
     1  1,0,0,0,1,1,
     2  0,3,0,0,0,0,
     3  0,0,3,0,0,0,
     4  0,0,0,3,0,0,
     5  1,0,0,0,1,1,
     6  1,0,0,0,1,1  /
C-----------------------------------------------------------------------

      nprob = 0
      if (test .ne. 'BLDPRJ') then
         write (*,*)
     *   ' '//hm//' test string was not BLDPRJ in that routine'
         stop
      endif

C check the projectors for internal correctness

      do 4000 iproj=1,42
C there are only 42 niggli lattice types (ignoring triclinic)
        if (ngtype(iproj).gt.43) then
            write (*,*)
     *      ' '//hm//' bad ngtype ',iproj,ngtype(iproj)
            nprob = nprob + 1
        endif
C check that the lattice type is present
        if (lattyp(iproj).eq.' ') then
           write (*,*) ' '//hm//' blank lattyp ',iproj
           nprob = nprob + 1
        endif
C check that the normalizers are in the range in Paciorek and Bonin
        if (pjnorm(iproj) .le. 0 .or. pjnorm(iproj).gt. 60) then
           write (*,*) ' '//hm//' bad pjnorm ',iproj,' ',pjnorm(iproj)
           nprob = nprob + 1
        endif
C check that the projectors are symmetrical matrices
        do 1000 i=1,5
        do 1000 j=i+1,6
          if (zprj(j,i,iproj) .ne. zprj(i,j,iproj)) then
             write (*,*) ' '//hm//' bad projector ',iproj,' ',i,' ',j
             write (*,*) zprj(j,i,iproj),' ',zprj(i,j,iproj)
             nprob = nprob + 1
          endif
 1000   continue
        do 2000 i=1,36
C check that the actual projector has no value greater than 1.0
          if (abs(projct(i,iproj)/(pjnorm(iproj))).gt. 1.00001) then
             write (*,*) ' '//hm//' bad projector, value > 1.0 '
             write (*,*) ' '//hm,iproj,' pjnorm ',pjnorm(iproj),
     2    ' i ',i,' ',projct(i,iproj)
             write (*,*)
             nprob = nprob + 1
          endif
 2000    continue
C check that the projector is positive definite
        do 3500 i=1,3
        sum = 0.0
        do 3000 j=1,6
          sum = sum + abs(zprj(i,j,iproj))
 3000   continue
        if (sum .le. 0) then
           write (*,*) ' '//hm//' zero xyz row, ',iproj
           nprob = nprob + 1
        endif
 3500   continue
 4000 continue

C actually output the projectors

      NPROJ = 0
      do 5100 i=1,MIN(MAXPRJ,42)
         IF (NGTYPE(I) .GT. 0) THEN
            NPROJ = NPROJ + 1
            itdesg(NPROJ) = ngtype(i)
            chrlat(NPROJ) = lattyp(i)
            pjn(NPROJ)    = pjnorm(i)
         do 5000 ii=1,36
            prj(ii,nproj) = projct(ii,i)
 5000    continue
         ENDIF
 5100 continue

      if (nprob .gt. 0) then
         write (*,*) ' '//hm//' ',
     *     nprob,' problems were found with projectors'
         stop
      endif

      itemp = projct(1,1)
      projct(1,1) = itemp

      end

C**********************************************************************C
      SUBROUTINE BLDTRE (MXTREE,NVEC,X,IDIN,TREE,TEST)
C  BLDTRE is called once for each point to be loaded into its internal
C  data structure (TREE).  It builds the tree structure of Kalantari
C  and McDonald (IEEE Transactions on Software Engineering, v. SE-9,
C  pp. 631-634,1983) for the extremely fast retrieval of coordinates.


C   MXTREE is the largest index that is allowed in the array TREE.
C   TREE   must be the order of 9-10 times the number of points to be
C          included.
C   X      is an input point's coordinates.
C   IDIN   is an arbitrary integer input, which will often be an array
C          index to be retrieved later. BLDTRE does not examine IDIN.
C   TEST   must be the string 'BLDTRE' -- it is used to make sure
C          that the number of formal parameters is correct.

C  TREE is used by NEARST and by INSPHR to find the nearest neighbor
C  to a probe point.  To initialize (or reinitialize) a TREE, set
C  TREE(1) equal to 0.0


      include 'ITERATE.cmn'
      CHARACTER*6 TEST
      LOGICAL DEBUG
      REAL TREE(MXTREE)
      INTEGER RMAX
      DATA LINK,RMAX,ID,ICHILD  /1,2,3,4/
      DATA DEBUG /.false./
C----------------------------------------------------------------------C

      IF (TEST .NE. 'BLDTRE' .AND. TEST .NE. 'bldtre')
     2     CALL BADCAL (TEST,'BLDTRE')

C   THE NODE SIZE IS 4 PLUS THE SIZE OF THE VECTOR
      NODSIZ = 4+NVEC
      IPOINT = 2
      IF (TREE(1) .GT. 0) THEN
          IFREE = TREE(1)
      ELSE
         IFREE = 2
      ENDIF
      tree(ifree) = 0
 1000 CONTINUE
      IF (DEBUG) WRITE (*,*)
     2 ' '//hm//' AFTER 1000 IN BLDTRE, IPOINT,TREE(IPOINT) ',
     3   IPOINT,TREE(IPOINT)
      IF (DEBUG) WRITE (*,*) ' '//hm//' IFREE,TREE(1),TREE(2) ',
     2     IFREE,INT(TREE(1)),INT(TREE(2))

      IF (TREE(IPOINT) .EQ. 0) THEN
         IF (DEBUG) WRITE (*,*)
     *   ' '//hm//' A NEW NODE IS BEING ALLOCATED'
         IPOINT = IFREE
         TREE(IPOINT) = -1
         TREE(IPOINT+LINK) = -1
         TREE(IPOINT+ID) = IDIN
         CALL CPYVN (NVEC,X,TREE(IPOINT+ICHILD))
         TREE(1) = IFREE + NODSIZ
         RETURN
      ELSEIF (TREE(IPOINT) .EQ. -1) THEN
         IF (DEBUG) WRITE (*,*)
     *   ' '//hm//' RIGHT CHILD OF NODE IS BEING FILLED'
         TREE(IPOINT) = IFREE
         TREE(IFREE+LINK) = -1
         TREE (IFREE+RMAX) = -1.0
         TREE(IFREE+ID) = IDIN
         CALL CPYVN (NVEC,X,TREE(IFREE+ICHILD))
         IFREE = IFREE + NODSIZ
         TREE(1) = IFREE
         RETURN
      ELSE
         IRIGHT = TREE(IPOINT)
         DL = TREELN (NVEC,X,TREE(IPOINT+ICHILD))
         DR = TREELN (NVEC,X,TREE(IRIGHT+ICHILD))
         IF (DEBUG) WRITE (*,*) ' '//hm//' DL,DR ',DL,DR
         IF (DR .GT. DL) THEN
             IF (DEBUG) WRITE (*,*)
     *       ' '//hm//' ',DR,DL,IPOINT,LINK,RMAX
             IF (TREE(IPOINT+LINK) .LE. 0) THEN
                TREE(IPOINT+RMAX) = DL
                TREE(IPOINT+LINK) = IFREE
                IPOINT = IFREE
             ELSE
                TREE(IPOINT+RMAX) =
     2                 MAX(DL,TREE(IPOINT+RMAX))
                IPOINT = TREE(IPOINT+LINK)
             ENDIF

         ELSE
             IF (TREE(IRIGHT+LINK) .LE. 0) THEN
                 TREE(IRIGHT+RMAX) = DR
                 TREE(IRIGHT+LINK) = IFREE
                 IPOINT = IFREE
              ELSE
                 TREE(IRIGHT+RMAX) =
     2                  MAX(DR,TREE(IRIGHT+RMAX))
                 IPOINT = TREE(IRIGHT+LINK)
              ENDIF

         ENDIF
         GO TO 1000
      ENDIF
      END


C***********************************************************************
      SUBROUTINE CHKVEC(V)
C  Check that a g6 vector represents a valid cell. Currently, it only
C  checks that the cell edges are real

      include 'ITERATE.cmn'
      REAL V(6)
C-----------------------------------------------------------------------
      NBAD = 0
      DO 1000 I=1,3
         IF(V(I).LE. 0.0) THEN
            WRITE (*,*) ' '//hm//' BAD VECTOR, I=',I,' ',V(I)
            NBAD = NBAD + 1
         ENDIF
 1000 CONTINUE
      IF (NBAD .GT. 0) STOP
      END

C**********************************************************************C
      SUBROUTINE CPYVN (NVEC,X,Y)
C----COPY A VECTOR X INTO A VECTOR Y
      DIMENSION X(NVEC), Y(NVEC)
C----------------------------------------------------------------------C
      DO 1000 I=1,NVEC
      Y(I) = X(I)
 1000 CONTINUE
      END

C***********************************************************************
      SUBROUTINE CTOG6 (CV,CVE,G,GE,SIZE,ERRSIZ,RATIO,TEST)
C  Convert from a unit cell (edge lengths and angles) to a g6 vector,
C  also computing the errors in the g6 vector. The length of the vector
C  and the error in the length of the vector is what iterate is really
C  going to use.

      include 'ITERATE.cmn'
      CHARACTER *6 TEST
      REAL COSI(4:6)
      REAL C(6),CV(6),CE(6),CVE(6),G(6),GE(6)
C-----------------------------------------------------------------------

      IF (TEST .NE. 'CTOG6 ') THEN
         WRITE (*,*) ' '//hm//' TEST IS WRONG IN CTOG6'
         STOP
      ENDIF
      RAD = ATAN2(0.0,-1.0) / 180.0
      DO 1000 I=1,6
          C(I) = CV(I)
          CE(I) = CVE(I)
 1000 CONTINUE
      DO 1100 I=4,6
         C(I) = C(I) * RAD
         CE(I) = CE(I) * RAD
 1100 CONTINUE

      DO 1500 I=1,3
         J = I + 3
         G(I) = C(I)*C(I)
         COSI(J) = COS(C(J))
         IF (C(I).NE.0.0) THEN
            G(J)=2.0*C(1)*C(2)*C(3)*COSI(J)/C(I)
         ELSE
            G(J) = 0.0
         ENDIF
         GE(I) = 2.*ABS(C(I)*CE(I))
C         WRITE (*,*) ' '//hm//' GE(I) ',I,' ',GE(I)
 1500 CONTINUE

      DO 2000 I=1,3
         J = I + 3
         IF (J.EQ.4) THEN
            I1 = 2
            I2 = 3
         ELSEIF (J.EQ.5) THEN
            I1 = 1
            I2 = 3
         ELSE
            I1 = 1
            I2 = 2
         ENDIF
         GE(J) = 2.0*SQRT(G(I1)*(COSI(J)*CE(I2))**2 +
     2                    G(I2)*(COSI(J)*CE(I1))**2 +
     3                    G(I1)*G(I2)*(SIN(C(J))*CE(J))**2)

C
C   NOTE THE UNITS IMBALANCE ABOVE
C

C      WRITE (*,*) ' '//hm//' GE(J) ',J,' ',GE(J)
 2000 CONTINUE

      SIZE = 0.0
      ERRSIZ = 0.0
      DO 3000 I=1,6
         SIZE = SIZE + G(I)*G(I)
         ERRSIZ = ERRSIZ + GE(I)*GE(I)
 3000 CONTINUE

      SIZE = SQRT(SIZE)
      ERRSIZ = SQRT(ERRSIZ)
      RATIO = ERRSIZ / SIZE
      END

C***********************************************************************
      FUNCTION DOTVN (N,V1,V2)
C compute a dot product
      REAL V1(N),V2(N)
C-----------------------------------------------------------------------
      DOTVN = 0.0
      DO 1000 I=1,N
         DOTVN = DOTVN + V1(I)*V2(I)
 1000 CONTINUE
      END

C***********************************************************************
      SUBROUTINE DRMV6 (V1,M,V2)
      REAL V1(6),V2(6)
      REAL M(36)
      DOUBLE PRECISION SUM
C-----------------------------------------------------------------------
      DO 3000 I=1,6
      SUM = 0.0D0
      DO 2000 J=1,6
         SUM = SUM + DBLE(M(6*(I-1)+J))*DBLE(V1(J))
 2000 CONTINUE
      V2(I) = SUM
 3000 CONTINUE
      END


C***********************************************************************
      logical function G6TOC (G,C,TEST)
C compute the normal unit cell parameters from a given g6 vector

      include 'ITERATE.cmn'
      CHARACTER *6 TEST
      REAL G(6),C(6)
C-----------------------------------------------------------------------
      IF (TEST .NE. 'G6TOC ') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN G6TOC'
         STOP
      ENDIF
      g6toc = .true.
      DO 900 I=1,3
         IF (G(I) .LE. 0.0) THEN
            WRITE (*,*) ' '//hm//' G(I)<=0, I=',I,'  ',G(I)
            g6toc = .false.
         ENDIF
  900 CONTINUE

      DO 1000 I=1,3
 1000 C(I) = SQRT(G(I))
      AC = 0.5*G(4)/C(2)/C(3)
      IF (ABS(AC) .LE. 1.0) THEN
         C(4) = 57.296*ACOS(AC)
      ELSE
         C(4) = 0.0
         WRITE (*,*) ' '//hm//' ARG>1.0 C(4) ',AC,1.0-ABS(AC)
            g6toc = .false.
      ENDIF

      AC = 0.5*G(5)/C(1)/C(3)
      IF (ABS(AC) .LE. 1.0) THEN
         C(5) = 57.296*ACOS(AC)
      ELSE
         C(5) = 0.0
         WRITE (*,*) ' '//hm//' ARG>1.0 C(5) ',AC,1.0-ABS(AC)
            g6toc = .false.
      ENDIF

      AC = 0.5*G(6)/C(1)/C(2)
      IF (ABS(AC) .LE. 1.0) THEN
         C(6) = 57.296*ACOS(AC)
      ELSE
         C(6) = 0.0
         WRITE (*,*) ' '//hm//' ARG>1.0 C(6) ',AC,1.0-ABS(AC)
            g6toc = .false.
      ENDIF
      END

C***********************************************************************
      subroutine g6tor3 (g6,m3)
C compute the normal unit cell parameters from a given g6 vector

      include 'ITERATE.cmn'
      logical pcmnt_
      integer i,j
      real g6(6,6),m3(3,3)
C-----------------------------------------------------------------------

      do 2000 i=1,3
         do 1000 j=1,3
            if (g6(i,j) .lt. -1.0e-6) then
               if (ostyle.ne.'CIF ') then
                 write (*,*)
     *           ' '//hm//' negative element in upper left of g6'
               else
                 cifres=pcmnt_(' negative square in g6 matrix')
               endif
            elseif (g6(i,j) .lt. 1.0e-6) then
               m3(i,j) = 0.0
            else
               m3(i,j) = sqrt(g6(i,j))
            endif
 1000 continue
 2000 continue


      do 3000 i=1,3
      call gtr3sn(m3(i,1),m3(i,2),m3(i,3), g6(i,4),g6(i,5),g6(i,6))
 3000 continue

      if (abs(g6(5,5)) .gt. 1.0e-6) then
         if (g6(5,5)*(m3(1,1)*m3(3,3)+m3(1,3)*m3(3,1)) .lt. 0.0) then
            do 4000 i= 1,3
               m3(3,i) = -m3(3,i)
 4000       continue
         endif
      elseif (abs(g6(5,4)) .gt. 1.0e-6) then
         if (g6(5,4)*(m3(1,2)*m3(3,3)+m3(1,3)*m3(3,2)) .lt. 0.0) then
            do 4100 i= 1,3
               m3(3,i) = -m3(3,i)
 4100       continue
         endif
      elseif (abs(g6(5,6)) .gt. 1.0e-6) then
         if (g6(5,6)*(m3(1,1)*m3(3,2)+m3(1,2)*m3(3,1)) .lt. 0.0) then
            do 4200 i= 1,3
               m3(3,i) = -m3(3,i)
 4200       continue
         endif
      endif

      if (abs(g6(6,6)) .gt. 1.0e-6) then
         if (g6(6,6)*(m3(1,1)*m3(2,2)+m3(1,2)*m3(2,1)) .lt. 0.0) then
            do 5000 i= 1,3
               m3(2,i) = -m3(2,i)
 5000       continue
         endif
      elseif (abs(g6(6,5)) .gt. 1.0e-6) then
         if (g6(6,5)*(m3(1,1)*m3(2,3)+m3(1,3)*m3(2,1)) .lt. 0.0) then
            do 5100 i= 1,3
               m3(2,i) = -m3(2,i)
 5100       continue
         endif
      elseif (abs(g6(6,4)) .gt. 0.0) then
         if (g6(6,4)*(m3(1,2)*m3(2,3)+m3(1,3)*m3(2,2)) .lt. 0.0) then
            do 5200 i= 1,3
               m3(2,i) = -m3(2,i)
 5200       continue
         endif
      endif
      end

      subroutine gtr3sn(e1,e2,e3, g4,g5,g6)
      if (e1 .ne. 0.0) then
         e2 = unitsn(g6)*e2
         e3 = unitsn(g5)*e3
      elseif (e3 .ne. 0.0) then
         e3 = unitsn(g4)*e3
      endif
      end

C***********************************************************************
      LOGICAL FUNCTION INPCEL (LATSYM,CV,CE,eof)
C get the input lattice type, cell, and errors in the cell parameters

      include 'ITERATE.cmn'
      logical char_
      logical numb_
      EXTERNAL OKCELL
      LOGICAL OKCELL
      CHARACTER *1 LATSYM
      REAL CV(6),CE(6)
      logical eof
C-----------------------------------------------------------------------

      eof = .false.
      inpcel = .true.
      cifeid = '.'
      cifsgs = 'P'
 1000 CONTINUE
      IF (querst.ne.'NO')
     *  WRITE (*,*) ' '//hm//' Input Xtal Lattice Centering '
      if (istyle.ne.'CIF ') then
        READ (*,'(A1)',end=9000) LATSYM
      else
        cifres = char_('_cell.entry_id',cifeid)
        if (.not.cifres) cifeid = '.'
        cifres = char_('_cell.space_group_name_H-M',cifsgs)
        LATSYM = 'P'
        if(cifres) LATSYM=cifsgs(1:1)
        if(cifsgs.eq.' ') cifsgs = 'P'
      endif
      IF (LATSYM .GE. 'a' .AND. LATSYM .LE. 'z')
     2 LATSYM = CHAR(ICHAR(LATSYM)-ICHAR('a')+ICHAR('A'))
      IF (querst.ne.'NO')
     *  WRITE (*,*) ' '//hm//' Input Cell Parameters'
      if (istyle.ne.'CIF ') then
      READ (*,*,end=9000) CV
      else
      do ii = 1,6
      CE(II) = 0.
      enddo
      cifres = numb_('_cell.length_a',cv(1),ce(1))
      if (.not.cifres) goto 9000
      cifres = numb_('_cell.length_b',cv(2),ce(2))
      cifres = numb_('_cell.length_c',cv(3),ce(3))
      cifres = numb_('_cell.angle_alpha',cv(4),ce(4))
      cifres = numb_('_cell.angle_beta',cv(5),ce(5))
      cifres = numb_('_cell.angle_gamma',cv(6),ce(6))
      endif
      INPCEL = OKCELL(LATSYM,CV,' TALK')
      IF (.NOT. INPCEL) GO TO 1000

      IF (querst.ne.'NO')
     * WRITE (*,*)
     * ' '//hm//' Input Standard Deviations of Cell Parameters'
      if (istyle.ne.'CIF ') then
      READ (*,*,end=9000) CE
      endif
      DO II = 1,6
      CE(II) = MAX(CE(II),ABS(CV(II))*5.E-7,1.E-4)
      ENDDO
      return
 9000 eof=.true.
      inpcel = .false.
      END

C**********************************************************************C
      FUNCTION INSPHR (MXTREE,NVEC,X,RADMAX,TREE,
     2    MXLIST,NLIST,LIST,IDLIST,TEST)
C         After the TREE has been constructed using BLDTRE,
C         INSPHR is used to retrieve all of the points within
C         RADMAX of the point X.  MXTREE is the maximum size of
C         TREE.  The indices of the found points are returned in
C         the array LIST; NLIST are returned, up to a maximum of
C         MXLIST.  For instance, TREE(LIST(3)) is the vector of
C         the third point found in the list.  IDLIST contains the
C         corresponding list of the input ID's.  If no points are
C         found within RADMAX of X, then INSPHR and NLIST are
C         returned as 0; otherwise they are returned as the
C         number of points found.  If more than MXLIST points were
C         found, then NLIST is returned equal to MXLIST, and
C         INSPHR is returned equal to -MXLIST.  TEST must be the
C         string 'INSPHR'.
C         See also BLDTRE and NEARST.
      include 'ITERATE.cmn'
      CHARACTER*6 TEST
      LOGICAL DEBUG
      INTEGER ISTAK(1000)
      REAL TREE(MXTREE)
      INTEGER LIST(MXLIST),IDLIST(MXLIST)
      INTEGER RMAX	
      DATA DEBUG /.FALSE./
      DATA LINK,RMAX,ID,ICHILD /1,2,3,4/
      DATA RIGHT,LEFT,END /111,112,113/
C----------------------------------------------------------------------C

      IF (TEST .NE. 'INSPHR' .AND. TEST .NE. 'insphr')
     2    CALL BADCAL (TEST,'INSPHR')

      ISTKP = 0
      NLIST = 0
      IPOINT = 2
      CURMIN = RADMAX
      DIR = LEFT
      DIRPRV = RIGHT
      IPREV = IPOINT -1
      if (tree(1).le.0.0) go to 8000
 1000 CONTINUE
      IF (IPREV .EQ. IPOINT .AND. DIRPRV .EQ. DIR) THEN
         WRITE (*,*) ' '//hm//' INTERNAL ERROR IN INSPHR '
         WRITE (*,*)
     *   ' '//hm//' TREE POINTER DIDN''T CHANGE',IPOINT,' ',DIR
         STOP
      ELSEIF (IPOINT .EQ. 0) THEN
         WRITE (*,*)
     *   ' '//hm//' INTERNAL ERROR IN INSPHR, IPOINT = 0'
         STOP
      ENDIF
      IPREV = IPOINT
      DIRPRV = DIR
      IF (DEBUG) WRITE (*,*)  
     *  ' '//hm//' IN INSPHR 1000, IPOINT = ',IPOINT
      IF (TREE(IPOINT) .EQ. 0) THEN
         IF (DEBUG) WRITE (*,*)
     *  ' '//hm//' AT AN END WITH IPOINT = ',IPOINT
         DIR = END
      ELSEIF (DIR .EQ. RIGHT) THEN
         IRIGHT = TREE(IPOINT)
         IF (DEBUG) WRITE (*,*)
     *   ' '//hm//' WENT RIGHT WITH IPOINT ', IPOINT
         DR = TREELN (NVEC,TREE(IRIGHT+ICHILD),X)
         IF (DR .LT. CURMIN) THEN
            NLIST = NLIST + 1
            IF (NLIST .GT. MXLIST) GO TO 8000
            LIST(NLIST) = IRIGHT + ICHILD
            IDLIST(NLIST) = TREE(IRIGHT+ID)
         ENDIF
         IF (TREE(IRIGHT+LINK) .LE. 0) THEN
            IF (DEBUG) WRITE (*,*)
     *      ' '//hm//' ON RIGHT BRANCH, UNSTACK A POINT'
            DIR = END
         ELSEIF (TREE(IRIGHT+RMAX)+CURMIN .GT. DR) THEN
            IPOINT = TREE(IRIGHT+LINK)
            DIR = LEFT
         ELSE
             DIR = END
         ENDIF
      ELSE
         IF (DEBUG) WRITE (*,*) ' '//hm//' WENT LEFT, IPOINT ',IPOINT
         DIR = LEFT
         IF (TREE(IPOINT) .GT. 0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//' STACK ONE '
            CALL TRSTCK (IPOINT,ISTAK,ISTKP)
         ENDIF
         DL = TREELN (NVEC,TREE(IPOINT+ICHILD),X)
         IF (DL .LT. CURMIN) THEN
            NLIST = NLIST + 1
            IF (NLIST .GT. MXLIST) GO TO 8000
            LIST(NLIST) = IPOINT+ICHILD
            IDLIST(NLIST) = TREE(IPOINT+ID)
         ENDIF
         IF (TREE(IPOINT+LINK) .LE. 0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//' NO LEFT LINK, GO BACK'
            DIR = END
         ELSEIF (TREE(IPOINT+RMAX) .LT. 0.0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//' NO DESCENDING LEFT TREE'
            DIR = END
         ELSEIF (TREE(IPOINT+RMAX)+CURMIN .GT. DL) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//
     2        ' GOING TO GO DOWN ONE LEVEL ',
     3         ' IPOINT AND UPDATE ',IPOINT,' ',INT(TREE(IPOINT+LINK))
            IPOINT = TREE(IPOINT+LINK)
         ELSE
           IF (DEBUG) WRITE (*,*) ' '//hm//' NO CLOSER POINTS ON LEFT '
           IF (DEBUG) WRITE (*,*) ' '//hm//' CURMIN,TREE(IPOINT+RMAX),
     2     DL ',CURMIN,TREE(IPOINT+RMAX),DL
           DIR = END
         ENDIF
      ENDIF
      IF (DIR .EQ. END) THEN
         IF (IUNSTK(IPOINT,ISTAK,ISTKP) .LE. 0) GO TO 8000
         DIR = RIGHT
      ENDIF
      GO TO 1000
 8000 CONTINUE
      IF (NLIST .LE. MXLIST) THEN
         INSPHR = NLIST
      ELSE
         INSPHR = -MXLIST
         NLIST = MXLIST
      ENDIF
      END

C**********************************************************************C
      FUNCTION IUNSTK(NEXT,ISTAK,ISTKP)
C helper function for INSPHR and NEARST
      include 'ITERATE.cmn'
      INTEGER ISTAK(1000)
      LOGICAL DEBUG
      DATA DEBUG /.FALSE./
C----------------------------------------------------------------------C
      IF (DEBUG) WRITE (*,*) ' '//hm//' IUNSTK,NEXT ',NEXT
      IF (ISTKP .GT. 0) THEN
         NEXT = ISTAK(ISTKP)
         ISTKP = ISTKP-1
      ELSE
         NEXT = 0
      ENDIF
      IUNSTK = NEXT
      END

C***********************************************************************
      SUBROUTINE MKNORM (VI,Mnorm,VOUT,TEST)
C converts an input g6 vector to "normalized" form (Gruber's
C terminology) or "standard presentation" (Andrews and Bernstein's
C terminology)  and the corresponding transformation matrix
      include 'ITERATE.cmn'
      CHARACTER *6 TEST
      REAL VIN(6),VI(6),VOUT(6)
      real Mnorm(36),MAT(36),M1(36),MAT3(9)
      LOGICAL AGAIN
C-----------------------------------------------------------------------

      IF (TEST .NE. 'MKNORM') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN MKNORM'
         STOP
      ENDIF

      CALL CPYVN(6,VI,VIN)

      CALL RUNTMN(6,Mnorm)
      NCYCLE = 0
 1000 CONTINUE
      NCYCLE = NCYCLE + 1
      AGAIN =.FALSE.
      CALL ZEROS (36,MAT)
      IF ( (ABS(VIN(1)).GT.ABS(VIN(2))) .OR.
     2 (VIN(1).EQ.VIN(2) .AND. ABS(VIN(4)).GT.ABS(VIN(5))) ) THEN
         MAT(2) = 1
         MAT(7) = 1
         MAT(15) = 1
         MAT(23) = 1
         MAT(28) = 1
         MAT(36) = 1
         AGAIN = .TRUE.
      ELSEIF ( (ABS(VIN(2)).GT.ABS(VIN(3))) .OR.
     2 (VIN(2).EQ.VIN(3) .AND. ABS(VIN(5)).GT.ABS(VIN(6))) ) THEN
         MAT(1) = 1
         MAT(14) = 1
         MAT(9) = 1
         MAT(22) = 1
         MAT(35) = 1
         MAT(30) = 1
         AGAIN = .TRUE.
      ENDIF

      IF (AGAIN) THEN
         CALL mm6(MAT,Mnorm,M1)
         CALL CPYVN(36,M1,Mnorm)
         CALL RMV6(VIN,MAT,VOUT)
         CALL CPYVN(6,VOUT,VIN)
      ENDIF
      IF (AGAIN .AND. NCYCLE.LT.4) GO TO 1000

      NUMNEG = 0
      DO 2000 I=4,6
         IF (VIN(I).LT.0.0) NUMNEG = NUMNEG + 1
 2000 CONTINUE
      CALL RUNTMN(3,MAT3)
      DO 4000 I=4,6
      IF (NUMNEG.EQ.1) THEN
         IF(VIN(I).GE.0.) MAT(3*(I-4)+I-3) = -1.
C        MAT(6*(I-1)+I) = -SIGN(1.0,VIN(I))
      ELSEIF (NUMNEG.EQ.2) THEN
         IF(VIN(I).LT.0.) MAT(3*(I-4)+I-3) = -1.
C        MAT(6*(I-1)+I) = SIGN(1.0,VIN(I))
      ENDIF
 4000 CONTINUE
      call r3tog6(MAT3,MAT)
      CALL mm6(MAT,Mnorm,M1)
      CALL CPYVN(36,M1,Mnorm)
      CALL RMV6(VIN,MAT,VOUT)
      CALL CPYVN(6,VOUT,VIN)
      END

C***********************************************************************
      SUBROUTINE MKPRIM (LATSYM,GIN,M,GOUT,TEST)
C converts and input g6 vector to one corresponding to a primitive
C lattice and the corresponding transformation matrix

      include 'ITERATE.cmn'
      CHARACTER *6 TEST
      CHARACTER LATSYM
      REAL GIN(6),GOUT(6)
      REAL M(36)
C-----------------------------------------------------------------------

      IF (TEST .NE. 'MKPRIM') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN MKPRIM'
         STOP
      ENDIF
      CALL ZEROS (36,M)

      IF (LATSYM .EQ. 'P') THEN
         CALL RUNTMN(6,M)

      ELSEIF (LATSYM .EQ. 'I') THEN
         M(1) = 1
         M(8) = 1
         DO 1000 I=13,18
 1000    M(I) = 0.25
         M(20) = 1
         M(22) = 0.5
         M(24) = 0.5
         M(25) = 1
         M(29) = 0.5
         M(30) = 0.5
         M(36) = 1

      ELSEIF (LATSYM .EQ. 'F') THEN
         M(1) = 0.25
         M(2) = 0.25
         M(6) = 0.25
         M(7) = 0.25
         M(9) = 0.25
         M(11) = 0.25
         M(14) = 0.25
         M(15) = 0.25
         M(16) = 0.25
         M(21) = 0.5
         DO 2000 I=22,24
 2000    M(I)= 0.25
         M(26) = 0.5
         DO 2100 I=28,30
 2100    M(I) = 0.25
         M(31) = 0.5
         DO 2200 I=34,36
 2200    M(I) = 0.25
      ELSEIF (LATSYM .EQ. 'A') THEN
         M(1) = 1
         M(8) = 1
         M(14) = 0.25
         M(15) = 0.25
         M(16) = 0.25
         M(20) = 1
         M(22) = 0.5
         M(29) = 0.5
         M(30) = 0.5
         M(36) = 1
      ELSEIF (LATSYM .EQ. 'B') THEN
         M(1) = 1
         M(8) = 1
         M(13) = 0.25
         M(15) = 0.25
         M(17) = 0.25
         M(22) = 0.5
         M(24) = 0.5
         M(25) = 1
         M(29) = 0.5
         M(36) = 1
      ELSEIF (LATSYM .EQ. 'C') THEN
         M(1) = 1
         M(7) = 0.25
         M(8) = 0.25
         M(12) = 0.25
         M(15) = 1
         M(22) = 0.5
         M(23) = 0.5
         M(29) = 1
         M(31) = 1
         M(36) = 0.5
      ELSEIF (LATSYM .EQ. 'R') THEN
         DO 2300 I=1,36
 2300    M(I) = 1./9.
         M(1) = 4./9.
         M(5) = 2./9.
         M(6) = 2./9.
         M(11) = -1./9.
         M(12) = -1./9.
         M(14) = 4./9.
         M(16) = -2./9.
         M(17) = -1./9.
         M(18) = 2./9.
         M(19) = 2./9.
         M(20) = -4./9.
         M(21) = 2./9.
         M(22) = -1./9.
         M(23) = -2./9.
         M(25) = -4./9.
         M(26) = -4./9.
         M(27) = 2./9.
         M(28) = -1./9.
         M(30) = -5./9.
         M(31) = -4./9.
         M(32) = 2./9.
         M(33) = 2./9.
         M(34) = 2./9.
      ELSE
         WRITE (*,*) ' '//hm//' DID NOT FIND LATTICE SYMBOL ',LATSYM
         STOP
      ENDIF
      CALL RMV6(GIN,M,GOUT)
      END


C      CALL MKREFL (RATIO,MXTREE,TREE,NVMAX,V,MATREF,NV,GRED,'MKREFL')

C*********************************************************************** SUBROUTINE mm6(M1,M2,M3)
      SUBROUTINE MKREFL
     2 (DEBUG,RATIO,MXTREE,TREE,NVMAX,V,MATREF,NV,GRED,TEST)

C MKREFL performs the iterations to search out the various
C representations of a single lattice by a set of different unit cells.
C There are many ways to do this iteration, and this one may well
C not be optimal. It has been found to work well in practice as long
C as cutoffs are not too strict. It is clear that, in general, it does
C not find all of the possible unit cells within a particular radius
C in g6.

      include 'ITERATE.cmn'
      PARAMETER (MXSWTC=2)
      CHARACTER *6 TEST
      REAL TREE(MXTREE), V(6,NVMAX), GRED(6)
      REAL VT(6),VTT(6)
      real mt2k(36),mti2k(36),mt1k(36)
      real REFL(36,24)
      real SWTCH(36,MXSWTC)
      real MATREF(36,NVMAX)
      LOGICAL DEBUG
C-----------------------------------------------------------------------

      DATA (SWTCH(I,1),I=1,36)/
     1 1,0,0, 0,0,0,
     2 0,1,0, 0,0,0,
     3 1,0,1, 0,1,0,
     4 0,0,0, 1,0,1,
     5 2,0,0, 0,1,0,
     6 0,0,0, 0,0,1 /

      DATA (SWTCH(I,2),I=1,36)/
     1 1,0,0, 0,0,0,
     2 0,1,0, 0,0,0,
     3 1,1,1, 1,1,1,
     4 0,2,0, 1,0,1,
     5 2,0,0, 0,1,1,
     6 0,0,0, 0,0,1 /


      DATA (REFL(I,1),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,1,0,0,
     6 0,0,0,0,1,0,
     7 0,0,0,0,0,1 /

      DATA (REFL(I,2),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,-1,0,0,
     6 0,0,0,0,-1,0,
     7 0,0,0,0,0,1 /

      DATA (REFL(I,3),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,1,0,0,
     6 0,0,0,0,-1,0,
     7 0,0,0,0,0,-1 /

      DATA (REFL(I,4),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,-1,0,0,
     6 0,0,0,0,1,0,
     7 0,0,0,0,0,-1 /



      DATA (REFL(I,5),I=1,36) /
     2 0,1,0,0,0,0,
     3 1,0,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,0,1,0,
     6 0,0,0,1,0,0,
     7 0,0,0,0,0,1 /

      DATA (REFL(I,6),I=1,36) /
     2 0,1,0,0,0,0,
     3 1,0,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,0,-1,0,
     6 0,0,0,-1,0,0,
     7 0,0,0,0,0,1 /

      DATA (REFL(I,7),I=1,36) /
     2 0,1,0,0,0,0,
     3 1,0,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,0,-1,0,
     6 0,0,0,1,0,0,
     7 0,0,0,0,0,-1 /

      DATA (REFL(I,8),I=1,36) /
     2 0,1,0,0,0,0,
     3 1,0,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,0,1,0,
     6 0,0,0,-1,0,0,
     7 0,0,0,0,0,-1 /



      DATA (REFL(I,9),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,0,1,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,1,0,0,
     6 0,0,0,0,0,1,
     7 0,0,0,0,1,0 /

      DATA (REFL(I,10),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,0,1,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,-1,0,0,
     6 0,0,0,0,0,-1,
     7 0,0,0,0,1,0 /

      DATA (REFL(I,11),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,0,1,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,-1,0,0,
     6 0,0,0,0,0,1,
     7 0,0,0,0,-1,0 /

      DATA (REFL(I,12),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,0,1,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,1,0,0,
     6 0,0,0,0,0,-1,
     7 0,0,0,0,-1,0 /



      DATA (REFL(I,13),I=1,36) /
     2 0,0,1,0,0,0,
     3 0,1,0,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,0,1,
     6 0,0,0,0,1,0,
     7 0,0,0,1,0,0 /

      DATA (REFL(I,14),I=1,36) /
     2 0,0,1,0,0,0,
     3 0,1,0,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,0,-1,
     6 0,0,0,0,-1,0,
     7 0,0,0,1,0,0 /

      DATA (REFL(I,15),I=1,36) /
     2 0,0,1,0,0,0,
     3 0,1,0,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,0,-1,
     6 0,0,0,0,1,0,
     7 0,0,0,-1,0,0 /

      DATA (REFL(I,16),I=1,36) /
     2 0,0,1,0,0,0,
     3 0,1,0,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,0,1,
     6 0,0,0,0,-1,0,
     7 0,0,0,-1,0,0 /



      DATA (REFL(I,17),I=1,36) /
     2 0,1,0,0,0,0,
     3 0,0,1,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,1,0,
     6 0,0,0,0,0,1,
     7 0,0,0,1,0,0 /

      DATA (REFL(I,18),I=1,36) /
     2 0,1,0,0,0,0,
     3 0,0,1,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,-1,0,
     6 0,0,0,0,0,-1,
     7 0,0,0,1,0,0 /

      DATA (REFL(I,19),I=1,36) /
     2 0,1,0,0,0,0,
     3 0,0,1,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,-1,0,
     6 0,0,0,0,0,1,
     7 0,0,0,-1,0,0 /

      DATA (REFL(I,20),I=1,36) /
     2 0,1,0,0,0,0,
     3 0,0,1,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,1,0,
     6 0,0,0,0,0,-1,
     7 0,0,0,-1,0,0 /



      DATA (REFL(I,21),I=1,36) /
     2 0,0,1,0,0,0,
     3 1,0,0,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,0,0,1,
     6 0,0,0,1,0,0,
     7 0,0,0,0,1,0 /

      DATA (REFL(I,22),I=1,36) /
     2 0,0,1,0,0,0,
     3 1,0,0,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,0,0,-1,
     6 0,0,0,-1,0,0,
     7 0,0,0,0,1,0 /

      DATA (REFL(I,23),I=1,36) /
     2 0,0,1,0,0,0,
     3 1,0,0,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,0,0,-1,
     6 0,0,0,1,0,0,
     7 0,0,0,0,-1,0 /

      DATA (REFL(I,24),I=1,36) /
     2 0,0,1,0,0,0,
     3 1,0,0,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,0,0,1,
     6 0,0,0,-1,0,0,
     7 0,0,0,0,-1,0 /

      IF (TEST .NE. 'MKREFL') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN MKREFL'
         STOP
      ENDIF

      TREE(1) = 0
      TREE(2) = 0

      NV = 0
      DMIN = 0.1*SQRT(DOTVN(6,GRED,GRED))*RATIO
      GRMIN = GRED(1)**2+GRED(2)**2+GRED(3)**2
      DO 4000 ICYCLE=1,MXSWTC+1
      DO 2000 IREFL=1,24
         CALL RMV6(GRED,REFL(1,IREFL),VT)
C#
C#
C      accumulate matrices starting from the reduced vector
         call cpyvn (36,refl(1,irefl),mt2k)
C#
C#
C There is a theoretical upper bound on how much the sum of the
C edge lengths can change (3.0) in a single transformation. So
C by allowing their square to only change by a factor of 9.0
C (incremented to allow for error), we can limit the search.
         IF (VT(1)**2+VT(2)**2+VT(3)**2.GT. 10.0*GRMIN) THEN
            IF (DEBUG)WRITE (*,*) ' '//hm//' REJECT 1 ',GRMIN,VT
         ELSEIF (ABS(VT(4)**2/VT(2)/VT(3)) .GT. 1.2 .OR.
     2           ABS(VT(5)**2/VT(1)/VT(3)) .GT. 1.2 .OR.
     3           ABS(VT(6)**2/VT(1)/VT(2)) .GT. 1.2 ) THEN
            IF (DEBUG) WRITE (*,*)
     *      ' '//hm//' REJECT 2 ',VT(4)**2/VT(2)/VT(3),
     2                                   VT(5)**2/VT(1)/VT(3),
     3                                   VT(6)**2/VT(1)/VT(2)




         ELSE

C set things up so that on the first cycle, the original vector
C is stored, and nothing else happens. after that, the rest of the
C switch matrices are used (from the 4000 loop).
            IF (ICYCLE.EQ.1) THEN
               MXINNR = 1
C#
C#
               call cpyvn (36,mt2k,mti2k)
C#
C#
            ELSE
               MXINNR = 24
               CALL RMV6(VT,SWTCH(1,ICYCLE-1),VTT)
C#
C#
               call mm6(swtch(1,icycle-1),mt2k,mti2k)
C#
C#
               CALL CPYVN(6,VTT,VT)
            ENDIF
            DO 1000 INNER=1,MXINNR
               IF (ICYCLE.GT.1) THEN
                  CALL RMV6(VT,REFL(1,INNER),VTT)
                  CALL CPYVN(6,VTT,VT)
C#
C#
                  call mm6(refl(1,inner),mti2k,mt1k)
                  call cpyvn (36,mt1k,mti2k)
C#
C#
               else
                  call cpyvn (36,mti2k,mt1k)
               ENDIF
               NV1 = NV + 1
               IF (NEARST(DEBUG,
     2           MXTREE,6,VT,DMIN,TREE,NV1,ID,'NEARST') .EQ. 0)
     3          THEN
                  NV = NV + 1
                  CALL BLDTRE (MXTREE,6,VT,NV,TREE,'BLDTRE')
                  CALL CPYVN(6,VT,V(1,NV))
C#
C#
                  call cpyvn(36,mt1k,matref(1,nv))
C#
C#
                  IF (NV .EQ. NVMAX) GO TO 8000
               ENDIF
 1000       CONTINUE
         ENDIF
 2000 CONTINUE
 4000 CONTINUE
 8000 CONTINUE
      if (debug) then
         WRITE (*,*) ' '//hm//' NV IN MKREFL ',NV
      endif

      END

C***********************************************************************
      SUBROUTINE mm6(M1,M2,M3)
C multiply two matrices (6x6), both in a linear array

      REAL M1(36),M2(36),M3(36)
C-----------------------------------------------------------------------
      CALL ZEROS (36,M3)
      DO 3000 I36=1,36
      IROW = (I36+5)/6
      ICOL = MOD(I36-1,6)+1
      K = ICOL-6
      DO 2000 J=6*IROW-5,6*IROW
        K = K + 6
           M3(I36) = M3(I36) + M1(J)*M2(K)
 1000   CONTINUE
 2000 CONTINUE
 3000 CONTINUE
      END


C***********************************************************************
      LOGICAL FUNCTION NEARRD (A,R,SIG)

C test whether a cell (really a g6 vector) is nearly Buerger reduced

      include 'ITERATE.cmn'
      REAL A(6),R(6), B(6)
      LOGICAL DEBUG
      DATA DEBUG /.FALSE./
C R IS THE REDUCED CELL VECTOR
C check that the vector is near to reduced
C-----------------------------------------------------------------------
      DO 1000 I=1,3
         IF (A(I) .LT. 1.0) THEN
            NEARRD = .FALSE.
            RETURN
         ENDIF
 1000 CONTINUE
      DO 2000 I=1,6
         B(I) = A(I)
 2000 CONTINUE

      IF (B(1) .GT. B(2)) THEN
         BT = B(1)
         B(1) = B(2)
         B(2) = BT
         BT = B(4)
         B(4) = B(5)
         B(5) = BT
      ENDIF
      IF (B(2) .GT. B(3)) THEN
         BT = B(2)
         B(2) = B(3)
         B(3) = BT
         BT = B(5)
         B(5) = B(6)
         B(6) = BT
      ENDIF
      IF (B(1) .GT. B(2)) THEN
         BT = B(1)
         B(1) = B(2)
         B(2) = BT
         BT = B(4)
         B(4) = B(5)
         B(5) = BT
      ENDIF



      NEARRD = .TRUE.
      DO 4000 I=1,3
         IF (B(I)-R(I) .GT. 5.0*SIG) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//
     2        ' NOT NEAR BUERGER ',I,B(I),R(I)
            NEARRD = .FALSE.
            RETURN
         ENDIF
 4000 CONTINUE

      IF (ABS(B(4))/SQRT(B(2)*B(3)) .GT. 1.0+5.0*SIG) THEN
         IF (DEBUG) WRITE (*,*)' '//hm//' ALPHA BAD '
         NEARRD = .FALSE.
         RETURN
      ELSEIF (ABS(B(5))/SQRT(B(1)*B(3)) .GT. 1.0+5.0*SIG) THEN
         IF (DEBUG) WRITE (*,*)' '//hm//' BETA BAD '
         NEARRD = .FALSE.
         RETURN
      ELSEIF (ABS(B(6))/SQRT(B(1)*B(2)) .GT. 1.0+5.0*SIG) THEN
         IF (DEBUG) WRITE (*,*)' '//hm//' GAMMA BAD '
         NEARRD = .FALSE.
         RETURN
      ENDIF



      END

C**********************************************************************C
      FUNCTION NEARST (debug,MXTREE,NVEC,X,RADMAX,TREE,IP,IDOUT,TEST)
C         The parameters are the same as those of INSPHR, except
C         that only the nearest point is found.  If no points are
C         found within RADMAX of X, then IP and NEARST are
C         returned as 0; otherwise they are returned as the index
C         in TREE of the coordinates of the vector of the point
C         closest to X.  TEST must be the string 'NEARST'.  ID is
C         returned as the ID of the nearest point (if the value
C         input as IDIN in BLDTRE was the index)
C         See also BLDTRE and INSPHR
      include 'ITERATE.cmn'
      CHARACTER*6 TEST
      LOGICAL DEBUG
      parameter (maxstk=1000)
      INTEGER ISTAK(maxstk)
      REAL TREE(MXTREE)
C      DATA DEBUG /.FALSE./
      DATA LINK,RMAX,ID,ICHILD /1,2,3,4/
      DATA RIGHT,LEFT,END /111,112,113/
C----------------------------------------------------------------------C

      IF (TEST .NE. 'NEARST' .AND. TEST .NE. 'nearst')
     2    CALL BADCAL (TEST,'NEARST')

      ISTKP = 0
      IP = 0
      IPOINT = 2
      CURMIN = RADMAX
      DIR = LEFT
      DIRPRV = RIGHT
      IPREV = IPOINT -1
      if (tree(1).le.0.0) go to 8000
 1000 CONTINUE
      IPREV = IPOINT
      DIRPRV = DIR
      IF (DEBUG) WRITE (*,*)  ' '//hm//
     2   ' IN NEARST 1000, IPOINT = ',IPOINT
      IF (TREE(IPOINT) .EQ. 0) THEN
         IF (DEBUG) WRITE (*,*) ' '//hm//
     2   ' AT AN END WITH IPOINT = ',IPOINT
         DIR = END
      ELSEIF (DIR .EQ. RIGHT) THEN
         IRIGHT = TREE(IPOINT)
         IF (DEBUG) WRITE (*,*) ' '//hm//
     2     ' WENT RIGHT WITH IPOINT ', IPOINT
         DR = TREELN (NVEC,TREE(IRIGHT+ICHILD),X)
         IF (DR .LT. CURMIN) THEN
            CURMIN = DR
            IP = IRIGHT + ICHILD
            IDOUT = TREE(IRIGHT+ID)
         ENDIF
         IF (TREE(IRIGHT+LINK) .LE. 0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//
     2        ' ON RIGHT BRANCH, UNSTACK A POINT'
            DIR = END
         ELSEIF (TREE(IRIGHT+RMAX)+CURMIN .GT. DR) THEN
            IPOINT = TREE(IRIGHT+LINK)
            DIR = LEFT
         ELSE
             DIR = END
         ENDIF
      ELSE
         IF (DEBUG) WRITE (*,*) ' '//hm//' WENT LEFT, IPOINT ',IPOINT
         DIR = LEFT
         IF (TREE(IPOINT) .GT. 0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//' STACK ONE '
            CALL TRSTCK(IPOINT,ISTAK,ISTKP)
         ENDIF
         DL = TREELN (NVEC,TREE(IPOINT+ICHILD),X)
         IF (DL .LT. CURMIN) THEN
            CURMIN = DL
            IP = IPOINT+ICHILD
            IDOUT = TREE(IPOINT+ID)
         ENDIF
         IF (TREE(IPOINT+LINK) .LE. 0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//' NO LEFT LINK, GO BACK'
            DIR = END
         ELSEIF (TREE(IPOINT+RMAX) .LT. 0.0) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//' NO DESCENDING LEFT TREE'
            DIR = END
         ELSEIF (TREE(IPOINT+RMAX)+CURMIN .GT. DL) THEN
            IF (DEBUG) WRITE (*,*) ' '//hm//
     2         ' GOING TO GO DOWN ONE LEVEL ',
     3         ' IPOINT AND UPDATE ',IPOINT,' ',INT(TREE(IPOINT+LINK))
            IPOINT = TREE(IPOINT+LINK)
         ELSE
           IF (DEBUG) WRITE (*,*) ' '//hm//' NO CLOSER POINTS ON LEFT '
           IF (DEBUG) WRITE (*,*) ' '//hm//' CURMIN,TREE(IPOINT+RMAX),
     2       DL ',CURMIN,TREE(IPOINT+RMAX),DL
           DIR = END
         ENDIF
      ENDIF
      IF (DIR .EQ. END) THEN
         if (debug) write (*,*) ' '//hm//
     2     ' call unstk ',ipoint,istkp,maxstk
         IF (IUNSTK(IPOINT,ISTAK,ISTKP) .LE. 0) GO TO 8000
         DIR = RIGHT
      ENDIF
      GO TO 1000
 8000 CONTINUE
      NEARST = IP
      END

C***********************************************************************
      LOGICAL FUNCTION OKCELL (LATSYM,CV,talk)
C okcell determines if a particular cell makes sense. The variable
C talk is used to determine if output is printed.

      include 'ITERATE.cmn'
      logical pcmnt_
      PARAMETER (NLATT=16)
      CHARACTER *1 SYMLST(NLATT)
      CHARACTER *1 LATSYM
      REAL CV(6)
      character *(*) talk
      DATA SYMLST /'P','A','B','C','I','F','R','H',
     2             'p','a','b','c','i','f','r','h'/
C-----------------------------------------------------------------------
      OKCELL = .TRUE.


      IF (CV(4).GT.175.0 .OR. CV(5).GT.175.0 .OR. CV(6).GT.175.0) THEN
         if (talk .eq. ' TALK')
     2   WRITE (*,*)
     *   ' '//hm//' THE LATTICE ANGLES MUST BE LESS THAN 175 DEGREES'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(4).LT.5.0 .OR. CV(5).LT.5.0 .OR. CV(6).LT.5.0) THEN
         if (talk .eq. ' TALK')
     2   WRITE (*,*)
     *  ' '//hm//' THE LATTICE ANGLES MUST EXCEED 5.0 DEGREES'
         OKCELL = .FALSE.
      ENDIF

      DO 1000 I=1,NLATT
         IF (LATSYM.EQ.SYMLST(I)) THEN
            GO TO 1100
         ENDIF
 1000 CONTINUE
         if (talk .eq. ' TALK')
     2 WRITE (*,*)
     * ' '//hm//' XTAL CENTERING SYMBOL ',LATSYM,' IS NOT IMPLEMENTED'
      OKCELL = .FALSE.
 1100 CONTINUE

      IF (CV(4).GE.CV(5)+CV(6)) THEN
         if (talk .eq. ' TALK')
     2      WRITE (*,*)
     *      ' '//hm//' ERROR, ALPHA EXCEEDS BETA PLUS GAMMA'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(5).GE.CV(4)+CV(6)) THEN
         if (talk .eq. ' TALK')
     2   WRITE (*,*) ' '//hm//' ERROR, BETA EXCEEDS ALPHA PLUS GAMMA'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(6).GE.CV(4)+CV(5)) THEN
         if (talk .eq. ' TALK')
     2   WRITE (*,*) ' '//hm//' ERROR, GAMMA EXCEEDS ALPHA PLUS BETA'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(4)+CV(5)+CV(6) .GT. 355.0) THEN
         if (talk .eq. ' TALK') then
         WRITE (*,*) ' '//hm//
     2     ' THE SUM OF THE CELL ANGLES MUST BE LESS THAN'
         WRITE (*,*) ' '//hm//' 355.0 DEGREES'
         endif
         OKCELL = .FALSE.
      ENDIF

      DO 2000 I=1,3
         IF (CV(I) .LT. 1.0) THEN
            WRITE (*,*) ' '//hm//' THE AXIAL LENGTHS MUST EXCEED 1.0'
            OKCELL = .FALSE.
         ENDIF
 2000 CONTINUE

      IF (LATSYM .EQ. 'R' .OR. LATSYM .EQ. 'r' .OR.
     *  LATSYM .EQ. 'H' .OR. LATSYM .EQ. 'h') THEN
        IF (ABS(CV(4)-90.).GT.1. .OR.
     *    ABS(CV(5)-90.).GT.1. .OR.
     *    ABS(CV(6)-120.).GT.1. .OR.
     *    ABS(CV(1)-CV(2)).GT..5) THEN
          IF (LATSYM.EQ.'R' .OR. LATSYM.EQ.'r') THEN
            if (talk .eq. ' TALK') then
            if (ostyle.ne.'CIF ') then
            WRITE (*,*) ' '//hm//' RHOMBOHEDRAL CENTERING REQUIRES'
            WRITE (*,*) ' '//hm//
     2        ' A HEXAGONAL CELL (A=B, ALPHA=BETA=90,'
            WRITE (*,*) ' '//hm//' GAMMA =120)'
            else
            cifres=pcmnt_(' Non-hexagonal R lattice '//cifeid)
            endif
            endif
            IF (ABS(CV(1)-CV(2)) .LT. .5 .AND.
     *        ABS(CV(2)-CV(3)) .LT. .5 .AND.
     *        ABS(CV(1)-CV(3)) .LT. .5 .AND.
     *        ABS(CV(4)-CV(5)) .LT. 1. .AND.
     *        ABS(CV(5)-CV(6)) .LT. 1. .AND.
     *        ABS(CV(4)-CV(6)) .LT. 1.) THEN
              LATSYM = 'P'
              if (ostyle.ne.'CIF ') then
              WRITE (*,*) ' '//hm//
     2          ' PROCESSING AS PRMITIVE RHOMBOHEDRAL'
              else
              cifres=pcmnt_(' Processing as primitve rhombohedral')
              endif
            ELSE
              OKCELL = .FALSE.
            ENDIF
          ELSE
            if (talk .eq. ' TALK') then
            if (ostyle.ne.'CIF ') then
            WRITE (*,*) ' '//hm//' A HEXAGONAL CELL REQUIRES'
            WRITE (*,*) ' '//hm//' A=B, ALPHA=BETA=90,GAMMA =120'
            else
            cifres=pcmnt_(' non-hexagonal cell '//cifeid)
            endif
            endif
            OKCELL = .FALSE.
          ENDIF
        ENDIF
        IF (LATSYM .EQ. 'H' .OR. LATSYM .EQ. 'h') THEN
          WRITE(*,*) ' '//hm//' PROCESSING H AS P'
          LATSYM = 'P'
        ENDIF
      ENDIF
      RETURN
      END


      subroutine r3tog6(e3,g6)

C give a 3-space transformation matrix, determine the corresponding one
C in g6

      real e3(9),g6(36)
C-----------------------------------------------------------------------
C------------------------
C   upper left 3x3 block
C------------------------
      g6(1)=e3(1)**2
      g6(2)=e3(2)**2
      g6(3)=e3(3)**2

      g6(7)=e3(4)**2
      g6(8)=e3(5)**2
      g6(9)=e3(6)**2

      g6(13)=e3(7)**2
      g6(14)=e3(8)**2
      g6(15)=e3(9)**2
C------------------------
C   upper right 3x3 block
C------------------------
      g6(4)=e3(2)*e3(3)
      g6(5)=e3(1)*e3(3)
      g6(6)=e3(1)*e3(2)

      g6(10)=e3(5)*e3(6)
      g6(11)=e3(4)*e3(6)
      g6(12)=e3(4)*e3(5)

      g6(16)=e3(8)*e3(9)
      g6(17)=e3(7)*e3(9)
      g6(18)=e3(7)*e3(8)
C------------------------
C   lower left 3x3 block
C------------------------
      g6(19)=2*e3(4)*e3(7)
      g6(20)=2*e3(5)*e3(8)
      g6(21)=2*e3(6)*e3(9)

      g6(25)=2*e3(1)*e3(7)
      g6(26)=2*e3(2)*e3(8)
      g6(27)=2*e3(3)*e3(9)

      g6(31)=2*e3(1)*e3(4)
      g6(32)=2*e3(2)*e3(5)
      g6(33)=2*e3(3)*e3(6)
C------------------------
C   lower right 3x3 block
C------------------------
      g6(22)=e3(5)*e3(9) + e3(8)*e3(6)
      g6(23)=e3(4)*e3(9) + e3(7)*e3(6)
      g6(24)=e3(4)*e3(8) + e3(4)*e3(5)

      g6(28)=e3(2)*e3(9) + e3(8)*e3(3)
      g6(29)=e3(1)*e3(9) + e3(7)*e3(3)
      g6(30)=e3(1)*e3(8) + e3(7)*e3(2)

      g6(34)=e3(2)*e3(6) + e3(5)*e3(3)
      g6(35)=e3(1)*e3(6) + e3(4)*e3(3)
      g6(36)=e3(1)*e3(5) + e3(4)*e3(2)

      end


C***********************************************************************
      PROGRAM RED

      include 'ITERATE.cmn'
      logical init_
      logical data_
      logical ocif_
      logical pdata_
      logical pfile_
      logical ploop_
      logical pnumb_
      logical pcmnt_
      logical pchar_
      INTEGER MAXPRJ
      PARAMETER (MAXPRJ=42)
      PARAMETER (NVMAX=1000)
      PARAMETER (MXTREE=11*NVMAX)
      INTEGER ITDESG(MAXPRJ)
      CHARACTER *2 CHRLAT(MAXPRJ)
      REAL PJNORM(MAXPRJ)
C      common /cmprjn/ pjnorm
      real PRJ(36,MAXPRJ)
C      common /cmprj/ prj
      REAL P(36),AP(36)
      REAL CV(6),CE(6),G(6),GE(6),TG(6),AG(6),COUT(6)
      REAL MPRIM(36)
      real m3ptrd(9)
      real MRED(36)
      REAL GOUT(6),GRED(6),CRED(6)
      LOGICAL INPCEL
      LOGICAL NEARRD
      EXTERNAL INPCEL
      CHARACTER LATSYM

      REAL TREE(MXTREE)
C      common /cmtree/ tree
      REAL V(6,NVMAX)
C      common /cmv/ v
      real MATREF(36,NVMAX),m1(36),m2(36)
      real retcel(6),m3(9),m3t(9),m3ti(9)
      REAL VBEST(6),AVBEST(6)

      REAL DOTVN
      EXTERNAL DOTVN
      INTEGER NPROJ,I,J
      REAL RATIO,SIZE,ERRSIZ
      character*8 cellst(6)
      logical eof,dbgmkr
      logical g6toc, test
      logical debug
      data debug /.FALSE./


C-----------------------------------------------------------------------
      CALL GETENV('ITERATE_QUERY',querst)
      CALL GETENV('OUTPUT_STYLE',ostyle)
      CALL GETENV('INPUT_STYLE',istyle)
      if(istyle.eq.'CIF '.or.ostyle.eq.'CIF ') then
        cifres = init_(5,6,21,0)
      endif
      if(istyle.eq.'CIF ') then
        querst = 'NO'
        cifres = ocif_(' ')
        cifres = data_(' ')
      endif
      hm = ' '
      if (ostyle.eq.'CIF ') then
        hm = '#'
        tabl_ = .false.
        cifres = pfile_(' ')
        cifres = pdata_('G6_SEARCH')
        cifres = ploop_('_cell.entry_id')
        cifres = ploop_('_cell.id')
        cifres = ploop_('_cell.space_group_name_H-M')
        cifres = ploop_('_cell.Bravais_lattice_symbol')
        cifres = ploop_('_cell.length_a')
        cifres = ploop_('_cell.length_b')
        cifres = ploop_('_cell.length_c')
        cifres = ploop_('_cell.angle_alpha')
        cifres = ploop_('_cell.angle_beta')
        cifres = ploop_('_cell.angle_gamma')
        cifres = ploop_('_cell.unreduced_length_a')
        cifres = ploop_('_cell.unreduced_length_b')
        cifres = ploop_('_cell.unreduced_length_c')
        cifres = ploop_('_cell.unreduced_angle_alpha')
        cifres = ploop_('_cell.unreduced_angle_beta')
        cifres = ploop_('_cell.unreduced_angle_gamma')
      endif
      do 1 i=2,mxtree
        tree(i) = 2**30
    1 continue
      iunt0 = 0
      iunt1 = 1
      iunt2 = 2
      iunt3 = 3
      iunt10 = 10
      dbgmkr = .false.
      IF (querst.ne.'NO') THEN
        WRITE (*,*) ' '//hm//' BEFORE BLDPRJ'
      ENDIF
      CALL BLDPRJ (MAXPRJ,NPROJ,ITDESG,CHRLAT,PJNORM,PRJ,'BLDPRJ')
C      write (*,*) ' '//hm//' nproj = ',nproj


C loop until the input contains no more data


 1000 continue

      IF (INPCEL(LATSYM,CV,CE,eof)) THEN
         CALL CTOG6(CV,CE,G,GE,SIZE,ERRSIZ,RATIO,'CTOG6 ')
         if (ostyle.ne.'CIF ') then
           WRITE (*,'(A,A)') '  '//hm//' Input Lattice Symbol  ',
     2       latsym
           WRITE (*,*)
           WRITE (*,'(''  '//hm//' INPUT CELL AND ERRORS   '',
     2      ''    INPUT VECTOR AND ERRORS'')')
           WRITE (*,*)
         else
           cifres = pchar_(' ',cifeid)
           cifres = pchar_(' ','.')
           cifres = pchar_(' ',cifsgs)
           cifres = pchar_(' ',latsym)
         endif
         CALL WRCLV6 (CV,CE,G,GE,'WRCLV6')
         if (ostyle.ne.'CIF ') then
           WRITE (*,*)
         else
           do ii = 1,6
           cifres = pchar_(' ','.')
           enddo
         endif
         CALL WRSIZE (SIZE,ERRSIZ,RATIO,'WRSIZE')
         CALL MKPRIM (LATSYM,G,MPRIM,GOUT,'MKPRIM')
         call cpyvn (36,mprim,m1)
         CALL CHKVEC(GOUT)
         CALL RUNTMN(6,MRED)
         CALL REDUCE (GOUT,MRED,GRED,'REDUCE')
         call mm6 (mred,m1,m2)
         CALL CHKVEC(GRED)
         test = G6TOC (GRED,CRED,'G6TOC ')
         SIZE = SQRT(DOTVN(6,GRED,GRED))
         ERRSIZ = RATIO * SIZE
         if (ostyle.ne.'CIF ') then
           WRITE (*,*)
           WRITE (*,'(''  '//hm//' REDUCED CELL  '')')
           WRITE (*,'(''  '//hm//' Red. Cell   '',6F10.3)')CRED
           WRITE (*,'(''  '//hm//' Red. Vector '',6F9.2)') GRED
         endif
         CALL MKREFL (dbgmkr,RATIO,MXTREE,TREE,NVMAX,
     2     V,MATREF,NV,GRED,'MKREFL')
         if (querst .ne. 'NO') then
            if(ostyle.ne.'CIF ')write (*,*) ' '//hm//' nv=',nv
         endif
         do 1200 iv=1,nv
            call mm6(matref(1,iv),m2,m1)
            call cpyvn (36,m1,matref(1,iv))
 1200    continue
         if (debug) then
         do 1300 iv=1,nv
            write (*,'(1x,a3,i5,6f8.3)')
     *      ' '//hm//' ',iv,(v(ip,iv),ip=1,6)
 1300    continue
         endif


C after reducing the input cell and iterating to find various 
C nearly reduced cells, test the found cells using the projectors
C of Paciorek and Bonin

         nmatch = 0
         DO 4000 I=1,NPROJ
            nrej = 0
            DBEST = 1.0E20
            DBESTO = DBEST
            nbest = 0
C subtract the projector from the unit matrix (to give the "prep")
C which when multiplied times a vector gives the vector component
C not in the subspace defined by the projector
            DO 3100 J=1,36
               P(J) = REAL(PRJ(J,I))/PJNORM(I)
 3100          AP(J) = -P(J)
            DO 3200 J=1,36,7
               AP(J) = 1.0 + AP(J)
 3200       CONTINUE
            DO 3300  IV=1,NV
            CALL RMV6 (V(1,IV),AP,AG)
            CALL RMV6 (V(1,IV),P,TG)
            DTEST = DOTVN(6,AG,AG)/dotvn(6,tg,tg)
            IF (DTEST .GE. DBEST) THEN

            ELSEIF (TG(1) .LT. 1.0) THEN
            ELSEIF (TG(2) .LT. 1.0) THEN
            ELSEIF (TG(3) .LT. 1.0) THEN

            ELSEIF (TG(1)*TG(2)*TG(3) + 0.125*TG(6)*TG(4)*TG(5)
     2        + 0.125*TG(5)*TG(6)*TG(4) - 0.25*TG(5)*TG(2)*TG(5)
     3        - 0.25*TG(1)*TG(4)*TG(4) - 0.25*TG(6)*TG(6)*TG(3)
     4        .LE. 0.0) THEN
C           if the metric tensor is negative, the cell is bad
              NREJ = NREJ + 1
              IF (DEBUG) WRITE (*,*) ' '//hm//' FAILED METRIC TENSOR'

            ELSEIF (.NOT. NEARRD(V(1,IV),GRED,ERRSIZ)) THEN
              IF (DEBUG) WRITE (*,*) ' '//hm//' FAILED NEARRD'
            ELSE
               DBEST = DTEST
               CALL CPYVN(6,TG,VBEST)
               CALL CPYVN(6,AG,AVBEST)
               NBEST = IV
               IF (DEBUG) WRITE (*,*) ' '//hm//
     2           ' NBEST,DBEST ',NBEST,DBEST
            ENDIF
 3300       CONTINUE
 3400       CONTINUE

C output those cases that are acceptable

            cutoff = amax1(10.0,amin1(errsiz,999.0))
            if (nbest .ne. 0) dbest = sqrt(DOTVN(6,AVBEST,AVBEST))
            IF (nbest .ne. 0 .and. DBEST .LE. cutoff) THEN
               nmatch = nmatch + 1
               if(ostyle.ne.'CIF ') then
                 WRITE (*,*)
                 WRITE (*,*)
                 write (*,'(1x,I3,A,A,A,F7.2,A,A,A,I2,A)')
     2             i,'   ',chrlat(i),' '//hm//' ',dbest,
     3             ' = Distance Projected','    ',
     4             'Internat. Tables#(',itdesg(i),')'
               else
                 cifres = pchar_(' ',cifeid)
                 cifres = pnumb_(' ',float(i),0.)
                 cifres = pchar_(' ',cifsgs)
                 cifres = pchar_(' ',chrlat(i))
               endif
               if(ostyle.ne.'CIF') then
                 WRITE (*,'(''  '//hm//'    Vector       '',6F8.1)')
     2             VBEST
               endif
               if (G6TOC (VBEST,COUT,'G6TOC ')) then
                 if(ostyle.ne.'CIF ') then
                 WRITE (*,'(2x,'''//hm//
     2             '    cell'',/12x,2(3F8.2,'' ''))') COUT
                 else
                 do ii = 1,6
                   if(cout(ii).ge.0.)itemp = cout(ii)*100.+.5
                   if(cout(ii).lt.0.)itemp = cout(ii)*100.-.5
                   cifres = pnumb_(' ',float(itemp)/100.,0.)
                 enddo
                 endif
                 call g6tor3 (matref(1,nbest),m3t)
                 call inver (m3t,m3ti)
                 call unredc (itdesg(i),'   ',vbest,retcel,m3,cellst)
                 if(ostyle.eq.'CIF ') then
                   do ii = 1,6
                     cifres = pchar_(' ',cellst(ii))
                   enddo
                 endif
                 call trnspz (m3,m3t)
                 call wrcent(chrlat(i),retcel)
                 if (debug) then
                    write (*,'(3(1x,'' '//hm//' '',10x,3i3,/))')
     *              (int(m3(im3)),im3=1,9)
                    write (*,*)
                    write (*,*) ' '//hm//' matref ,nbest = ',
     2                nbest,' of ',nv
                    write (*,'(6(1x,'' '//hm//
     2                ' '',3x,3f8.3,2x,3f8.3/))')
     3                (matref(im,nbest),im=1,36)
                 endif
                 call g6tor3 (matref(1,nbest),m3ptrd)
                 call matmul (m3t,m3ptrd,m3ti)
                 if(ostyle.ne.'CIF ') then
                    write (*,*)' '//hm//'    transformation from ',
     2               '3-space original cell'
                    if (index('PR',chrlat(i)(2:2)) .ne. 0) then
                      write (*,*) ' '//hm//'    to final primitive cell'
                    else
                      write (*,*) ' '//hm//'    to final centered cell'
                    endif
                    write (*,'(3(1x,'' '//hm//' '',8x,3f7.3/))') m3ti
                 else
                   cifres= pcmnt_(char(0))
                 endif
               endif
            ENDIF
 4000    CONTINUE
      ENDIF
      if (nmatch .eq. 0) then
        if(ostyle.ne.'CIF ') then
            WRITE (*,*) ' '//hm//
     2        '*****************************************'
     3        //'*************'
            write (*,*) ' '//hm//
     2        ' NO MATCHES WERE FOUND WITH THE SPECIFIED '
     3        //'UNCERTAINTIES'
            WRITE (*,*) ' '//hm//
     2        '*****************************************'
     3        //'*************'
        else
            cifres = pcmnt_(char(0))
            cifres = pcmnt_(
     * ' NO MATCHES WERE FOUND WITH THE SPECIFIED UNCERTAINTIES')
        endif
      endif
      IF (querst.ne.'NO') THEN
        if (.not. eof) THEN
          if(ostyle.ne.'CIF ') then
            WRITE (*,*)
            WRITE (*,*)
          else
            cifres = pcmnt_(char(0))
            cifres = pcmnt_(' ')
          endif
          go to 1000
        endif
      ENDIF
      if (istyle.eq.'CIF ') then
        if(loop_) goto 1000
      endif
      if (ostyle.eq.'CIF ' .or. istyle.eq.'CIF ') then
        call close_
      endif
      END

C***********************************************************************
      SUBROUTINE REDUCE (VI,M,VOUT,TEST)

C for a given input vector, determine the reduced vector and the
C transformation from one to the other

      include 'ITERATE.cmn'
      CHARACTER *6 TEST

      REAL VI(6),VIN(6), VOUT(6)
      LOGICAL AGAIN
      real M(36),M1(36),M2(36),mnorm(36)
      REAL ZEROS(6)
      DATA ZEROS /6*0.0/
C-----------------------------------------------------------------------

      IF (TEST .NE. 'REDUCE') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN REDUCE'
         STOP
      ENDIF
      CALL CPYVN(6,VI,VIN)
      NCYCLE = 0
 1000 CONTINUE
      LAST = 0
      CALL MKNORM (VIN,Mnorm,VOUT,'MKNORM')
      call mm6 (mnorm,m,m2)
      CALL CPYVN(36,m2,m)
      CALL CPYVN(6,VOUT,VIN)
      CALL RUNTMN (6,M1)
      IF (ABS(VIN(4)) .GT. ABS(VIN(2))) THEN
          M1(14) = 1.0
          M1(16) = -SIGN (1.0,VIN(4))
          M1(20) = -2.0*SIGN(1.0,VIN(4))
          M1(30) = M1(16)
          AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
          LAST = 5
      ELSEIF (ABS(VIN(5)) .GT. ABS(VIN(1))) THEN
         M1(13) = 1
         M1(17) = -SIGN (1.0,VIN(5))
         M1(24) = M1(17)
         M1(25) = 2*M1(17)
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 6
      ELSEIF (ABS(VIN(6)) .GT. ABS(VIN(1))) THEN
         M1(7) = 1
         M1(12) = -SIGN(1.0,VIN(6))
         M1(23) = M1(12)
         M1(31) = 2*M1(12)
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 7
      ELSEIF (VIN(4)+VIN(5)+VIN(6)+ABS(VIN(1))+ABS(VIN(2)) .LT. 0.0)
     2   THEN
         DO 2000 I=13,18
 2000    M1(I) = 1
         M1(20) = 2
         M1(24) = 1
         M1(25) = 2
         M1(30) = 1
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 8
      ELSEIF ( (VIN(4).EQ.VIN(2) .AND. 2.0*VIN(5).LT.VIN(6)) .OR.
     2         (VIN(4).EQ.-ABS(VIN(2)) .AND. VIN(6).LT. 0.0) )
     3   THEN
         M1(14) = 1
         M1(16) = -SIGN(1.0,VIN(4))
         M1(20) = 2*M1(16)
         M1(30) = M1(16)
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 15
      ELSEIF ( (VIN(5).EQ.VIN(1) .AND. 2.0*VIN(4).LT.VIN(6)) .OR.
     2         (VIN(5).EQ.-ABS(VIN(1)) .AND. VIN(6).LT.0.0) )
     3    THEN
         M1(13) = 1
         M1(17) = -SIGN(1.0,VIN(5))
         M1(24) = M1(17)
         M1(25) = 2*M1(17)
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 16
      ELSEIF ( (VIN(6).EQ.VIN(1) .AND. 2.0*VIN(4).LT.VIN(5)) .OR.
     2         (VIN(6).EQ.-VIN(1) .AND. VIN(5).LT.0.0) ) THEN
         M1(7) = 1
         M1(12) = -SIGN(1.0,VIN(6))
         M1(23) = M1(12)
         M1(31) = 2*M1(12)
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 17
      ELSEIF ( (VIN(4)+VIN(5)+VIN(6)+ABS(VIN(1))+ABS(VIN(2)).EQ.0.0)
     2 .AND. ( 2.0*(ABS(VIN(1))+VIN(5))+VIN(6).GT.0.0) ) THEN
         DO 3000 I=13,18
 3000    M1(I) = 1
         M1(20) = 2
         M1(24) = 1
         M1(25) = 2
         M1(30) = 1
         AGAIN = .TRUE.
         CALL mm6(M1,M,M2)
         CALL CPYVN (36,M2,M)
         CALL RMV6(VIN,M1,VOUT)
         LAST = 18
      ELSE
         AGAIN = .FALSE.
         CALL CPYVN (6,VIN,VOUT)
         call cpyvn (36,m2,m)
      ENDIF

      CALL MKNORM (VOUT,Mnorm,VIN,'MKNORM')
         CALL mm6 (mnorm,m,m2)
         CALL CPYVN(36,M2,M)
         CALL CPYVN(6,VIN,VOUT)

      IF (VIN(1).LT. 0.0 .OR. VIN(2).LT.0.0 .OR. VIN(3).LT.0.0) THEN
         WRITE (*,*) ' '//hm//' NEG. SQ. AXIS ',NCYCLE
         CALL WRVEC6(VIN,ZEROS,'WRVEC6')
         CALL WRVEC6(VOUT,ZEROS,'WRVEC6')
         if (istyle.ne.'CIF ') READ (*,*)
      ENDIF

      NCYCLE = NCYCLE + 1
      IF (NCYCLE .LT. 25 .AND. AGAIN) GO TO 1000
      END

C***********************************************************************
      SUBROUTINE RMV6 (V1,M,V2)
      REAL V1(6),V2(6)
      REAL M(36)
C-----------------------------------------------------------------------
      DO 3000 I=1,6
      SUM = 0.0
      DO 2000 J=1,6
         SUM = SUM + M(6*(I-1)+J)*V1(J)
 2000 CONTINUE
      V2(I) = SUM
 3000 CONTINUE
      END

      function root (a)
      if (a.ge.0.0) then
         root = sqrt(a)
      else
         root = 0
      endif
      end

C***********************************************************************
      SUBROUTINE RUNTMN (N,M)
      INTEGER N
      REAL M(N,N)
C-----------------------------------------------------------------------
      DO 1000 I=1,N
      DO 1000 J=1,N
 1000 M(I,J) = 0.0
      DO 2000 I=1,N
 2000 M(I,I) = 1.0
      END


      function sqr (a)
      sqr = a*a
      end

C**********************************************************************C
      FUNCTION TREELN (NVEC,A,B)
C-----GET THE SEPARATION BETWEEN THE ENDS OF TWO VECTORS
C used by bldtre, insphr and nearst
      DIMENSION A(NVEC),B(NVEC)
C----------------------------------------------------------------------C
      SUM = 0.0
      DO 1000 I=1,NVEC
         SUM = SUM + (A(I)-B(I))**2
 1000 CONTINUE
      TREELN = SQRT(SUM)
      END

C**********************************************************************C
      SUBROUTINE TRSTCK (NEXT,ISTAK,ISTKP)

C helper routine for NEARST and INSPHR

      include 'ITERATE.cmn'
      INTEGER ISTAK(1000)
      LOGICAL DEBUG
      DATA DEBUG /.FALSE./
C----------------------------------------------------------------------C
      IF (DEBUG) WRITE (*,*) ' '//hm//
     2    ' IN TRSTCK, STACK POINTER,IPOINT ',
     3    ISTKP,NEXT
      ISTKP = ISTKP + 1
      ISTAK(ISTKP) = NEXT
      END

C**********************************************************************C
      subroutine unitmx(n,a)
C produce a unit matrix of order n
      real a(n,n)
C-----------------------------------------------------------------------
      do 1000 i=1,n
      do 1000 j=1,n
 1000 a(i,j) = 0.0
      do 2000 i=1,n
 2000 a(i,i) = 1.0
      end

C**********************************************************************C
      function unitsn (a)
C-----------------------------------------------------------------------
      if (a.ge.0.0) then
         unitsn = 1
      else
         unitsn = -1
      endif
      end

C**********************************************************************C
      subroutine unredc (itcase,rfcase,v,cell,m,cellst)

C given a g6 vector and a particular bravais lattice (designated either
C by itcase or rfcase, that is by the International tables or by Niggli
C and Roof's designations), compute the standard (often non-reduced)
C unit cell

      real m(9),v(6),cell(6)
      character *3 rfcase
      character *8 cellst(6),ctemp
      integer itcase
      real redv(6),altcel(6),v1(3),v2(3),v3(3),a(3),b(3),c(3),vtemp(3)
C-----------------------------------------------------------------------


      PI = 4.0*ATAN(1.0)
      MAXINT = 32768.0
      do 1000 I=1,6
        Cell(I) = Maxint
        altcel(i) = maxint
        cellst(i) = '.'
 1000 continue
      do 1100 I = 1,9
         M(I) = 0.0
 1100 continue
      do 1200 I=1,6
         redv(I) = root(v(I))
 1200 continue
      redv(4) = v(4)/(redv(2)*redv(3))
      redv(5) = v(5)/(redv(1)*redv(3))
      redv(6) = v(6)/(redv(1)*redv(2))
      IF (itcase .eq. 3 .or. rfcase .EQ. '44A') THEN
          Cell(1) = redv(1)
          CALL UNITMX(3,M)

      ELSEIF (itcase .eq. 5 .or. rfcase .EQ. '44B') THEN
          Cell(1) = root(4.0/3.0*v(1))
          CALL UNITMX(3,M)
          M(7) = 1.0
          M(2) = 1.0
          M(6) = 1.0

      ELSEIF (itcase .eq. 1 .or. rfcase .EQ. '44C') THEN
          Cell(1) = root(2*v(1))
          do 4400 i = 1,9
            M(I) = 1.0
 4400     CONTINUE
          M(4) = -1.0
          M(8) = -1.0
          M(3) = -1.0

      ELSEIF (itcase .eq. 11 .or. rfcase .EQ. '45A') THEN
         Cell(1) = redv(1)
         Cell(3) = redv(3)
         CALL UNITMX(3,M)

      ELSEIF (itcase .eq. 21 .or. rfcase .eq. '45B') THEN
         Cell(1) = redv(2)
         Cell(3) = redv(1)
                          M(4) = 1.0
                                           M(8) = 1.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 15 .or. rfcase .eq. '45C') THEN
         Cell(1) = redv(1)
         Cell(3) = root(4*v(3)-2*v(1))
         M(1) = 1.0
         M(5) = 1.0
         M(9) = 2.0
         M(3) = 1.0
         M(6) = 1.0

      ELSEIF (itcase .eq. 6 .or. rfcase .eq. '45D') THEN
         Cell(1) = root(2*v(1)+v(6))
         Cell(3) = root(2*v(1)+v(4))
         do 4500 I= 1,9
            M(I) = 1.0
 4500    CONTINUE
         do 4510 I = 1,9,4
            M(I) = 0.0
 4510    CONTINUE

      ELSEIF (itcase .eq. 7 .or. rfcase .eq. '45d') THEN
         Cell(1) = root(2*v(1)+v(6))
         Cell(3) = root(2*v(1)+v(4))
         do 4550 I=1,9,4
            M(I) = 1.0
 4550    CONTINUE
         M(4) = 0.0
         M(8) = 0.0
         M(3) = 0.0

      ELSEIF (itcase .eq. 18 .or. rfcase .eq. '45E') THEN
         Cell(1) = redv(1)
         Cell(3) = root( 2*v(3)-0.5*v(1))
         M(4) = -1.0
         M(7) = 1.0
         M(2) = 1.0
         M(5) = -1.0
         M(8) = -1.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 12 .or. rfcase .eq. '48A') THEN
         Cell(1) = redv(1)
         Cell(3) = redv(3)
         CALL UNITMX(3,M)

      ELSEIF (itcase .eq. 22 .or. rfcase .eq. '48B') THEN
         Cell(1) = redv(3)
         Cell(3) = redv(1)
         M(4) = 1.0
         M(8) = 1.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 9 .or. rfcase .eq. '49B') THEN
         Cell(1) = redv(3)
         Cell(4) = 1-2*Sqr(redv(1)/(2*redv(3)))
         altcel(1) = redv(1)
         altcel(3) = root(9*v(3)-3*v(1))
         M(1) = 1.0
         M(2) = -1.0
         M(5) = 1.0
         M(3) = -1.0
         M(6) = -1.0
         M(9) = 3.0

      ELSEIF (itcase .eq. 2 .or. itcase .eq. 4
     2     .or. rfcase .eq. '49x') THEN
         Cell(1) = redv(1)
         Cell(4) = v(4)/(2*v(1))
         altcel(1) = root(0.5*v(3)-0.25*v(4))
         altcel(3) = root(9*v(1)-3*Sqr(altcel(1)))
         M(1) = 1.0
         M(4) = -1.0
         M(2) = -1.0
         M(8) = 1.0
         do 4900 I= 3,9,3
            M(I) = -1.0
 4900    CONTINUE

      ELSEIF (itcase .eq. 24 .or. rfcase .eq. '49E') THEN
         Cell(1) = redv(3)
         Cell(4) = v(4)/(2*v(2))
         altcel(1) = root(0.5*v(3)-0.25*v(4))
         altcel(3) = redv(1)
         M(1) = 1.0
         M(4) = 2.0
         M(7) = 1.0
         M(5) = -1.0
         M(8) = 1.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 32 .or. rfcase .eq. '50C') THEN
         do 5000 I=1,3
            Cell(I) = redv(I)
 5000    CONTINUE
         CALL UNITMX(3,M)

      ELSEIF (itcase .eq. 36 .or. rfcase .eq. '50A') THEN
         Cell(1) = redv(1)
         Cell(2) = root(4*v(3)-v(1))
         Cell(3) = redv(2)
         M(1) = 1.0
         M(2) = -1.0
         M(8) = -2.0
         M(6) = 1.0

      ELSEIF (itcase .eq. 38 .or. rfcase .eq. '50B') THEN
         Cell(1) = redv(1)
         Cell(3) = redv(3)
         Cell(2) = root(4*v(2)-v(1))
         M(1) = 1.0
         M(2) = -1.0
         M(5) = -2.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 13 .or. rfcase .eq. '50D') THEN
         Cell(1) = root(2*v(1)+v(6))
         Cell(2) = root(2*v(1)-v(6))
         Cell(3) = redv(3)
         M(1) = 1.0
         M(4) = 1.0
         M(2) = -1.0
         M(5) = 1.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 23 .or. rfcase .eq. '50E') THEN
         Cell(1) = root(2*v(3)+v(4))
         Cell(2) = root(2*v(3)-v(4))
         Cell(3) = redv(1)
         M(4) = 1.0
         M(7) = 1.0
         M(5) = -1.0
         M(8) = 1.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 40 .or. rfcase .eq. '50F') THEN
         Cell(1) = redv(2)
         Cell(2) = root(4*v(3)-v(2))
         Cell(3) = redv(1)
         M(4) = 1.0
         M(5) = -1.0
         M(8) = -2.0
         M(3) = -1.0

      ELSEIF (itcase .eq. 16 .or. rfcase .eq. '51A') THEN
         Cell(1) = root(2*v(1)+v(6))
         Cell(2) = root(2*v(1)-v(6))
         Cell(3) = root(4*v(3)-Sqr(Cell(1)))
         M(1) = 1.0
         M(4) = -1.0
         M(5) = 1.0
         M(5) = 1.0
         M(8) = 2.0
         M(3) = -1.0
         M(6) = -1.0

      ELSEIF (itcase .eq. 26 .or. rfcase .eq. '51B') THEN
         Cell(1) = redv(1)
         Cell(2) = root(4*v(2)-v(1))
         Cell(3) = root(4*v(3)-v(1))
         M(1) = -1.0
         M(4) = 2.0
         M(2) = -1.0
         M(8) = 2.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 8 .or. rfcase .eq. '52A') THEN
         Cell(1) = root(2*v(1)+v(6))
         Cell(2) = root(2*v(1)+v(5))
         Cell(3) = root(2*v(1)+v(4))
         M(1) = 1.0
         M(7) = 1.0
         M(2) = 1.0
         M(5) = 1.0
         M(6) = 1.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 19 .or. rfcase .eq. '52B') THEN
         Cell(1) = redv(1)
         Cell(2) = root(2*v(2)-v(4))
         Cell(3) = root(2*v(3)+v(4)-v(1))
         M(1) = -1.0
         M(2) = -1.0
         M(5) = 1.0
         M(8) = 1.0
         M(6) = -1.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 42 .or. rfcase .eq. '52C') THEN
         Cell(1) = redv(1)
         Cell(2) = redv(2)
         Cell(3) = root(4*v(3)-v(1)-v(2))
         M(1) = 1.0
         M(5) = 1.0
         M(3) = -1.0
         M(6) = -1.0
         M(9) = -2.0

      ELSEIF (itcase .eq. 33 .or. rfcase .eq. '53A') THEN
         Cell(1) = redv(1)
         Cell(2) = redv(2)
         Cell(3) = redv(3)
         Cell(5) = v(5)/(2*Cell(1)*Cell(3))
         CALL UNITMX(3,M)

      ELSEIF (itcase .eq. 35 .or. rfcase .eq. '53B') THEN
         Cell(1) = redv(2)
         Cell(2) = redv(1)
         Cell(3) = redv(3)
         Cell(5) = v(4)/(2*Cell(1)*Cell(3))
         M(4) = 1.0
         M(2) = 1.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 34 .or. rfcase .eq. '53C') THEN
         Cell(1) = redv(1)
         Cell(2) = redv(3)
         Cell(3) = redv(2)
         Cell(5) = v(6)/(2*Cell(1)*Cell(3))
         M(1) = 1.0
         M(8) = 1.0
         M(6) = 1.0

      ELSEIF (itcase .eq. 39 .or. rfcase .eq. '54A') THEN
         Cell(1) = root(4*v(2)-v(1))
         Cell(2) = redv(1)
         Cell(3) = redv(3)
         Cell(5) = v(4)/(Cell(1)*Cell(3))
         M(1) = 1.0
         M(4) = 2.0
         M(2) = 1.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 41 .or. rfcase .eq. '54B') THEN
         Cell(1) = root(4*v(3)-v(2))
         Cell(2) = redv(2)
         Cell(3) = redv(1)
         Cell(5) = v(5)/(Cell(1)*Cell(3))
         M(4) = 1.0
         M(7) = 2.0
         M(5) = 1.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 37 .or. rfcase .eq. '54C') THEN
         Cell(1) = root(4*v(3)-v(1))
         Cell(2) = redv(1)
         Cell(3) = redv(2)
         Cell(5) = v(4)/(Cell(1)*Cell(3))
         M(1) = 1.0
         M(7) = 2.0
         M(2) = 1.0
         M(6) = 1.0

      ELSEIF (itcase .eq. 10 .or. itcase .eq. 14
     2     .or. rfcase .eq. '55A') THEN
       call v2cart(v,v1,v2,v3)
       call vecsum(v1,v2,A)
       call vecdif(v2,v1,B)
       call cpyvec(v3,c)
       Cell(1) = Sqrt(Dot(A,A))
       Cell(2) = Sqrt(Dot(B,B))
       Cell(3) = Sqrt(Dot(C,C))
       Cell(5) = -ABS(Dot(A,C)/Cell(1)/Cell(3))
       IF (Cell(1) .EQ. Cell(2)) THEN

          CALL UNITMX(3,M)
          M(4) = 1.0
          M(2) = -1.0
       ELSE

          M(4) = 1.0
          M(7) = 1.0
          M(5) = -1.0
          M(8) = 1.0
          M(3) = 1.0
       ENDIF


      ELSEIF (itcase .eq. 20 .or. itcase .eq. 25
     2    .or. rfcase .eq. '55B') THEN
       call v2cart(v,v1,v2,v3)
       call vecsum(v3,v2,A)
       call vecdif(v3,v2,B)
       call cpyvec (v1,c)
       Cell(1) = Sqrt(Dot(A,A))
       Cell(2) = Sqrt(Dot(B,B))
       Cell(3) = Sqrt(Dot(C,C))
       Cell(5) = -ABS(Dot(A,C)/Cell(1)/Cell(3))
       IF ((v(4).LT.0.0) .AND. (v(5).GE.0.0) .AND. (V(6).LT.0.0)) THEN

          M(4) = 1.0
          M(7) = 1.0
          M(5) = 1.0
          M(8) = -1.0
          M(3) = -1.0
       ELSE

          M(4) = 1.0
          M(7) = 1.0
          M(5) = -1.0
          M(8) = 1.0
          M(3) = 1.0
       ENDIF

      ELSEIF (itcase .eq. 28 .or. rfcase .eq. '56A') THEN
         call v2cart(v,A,C,v3)
         do 5600 I=1,3
            v3(I)=2*v3(I)
 5600    continue
         call vecdif(v3,A,B)
         Cell(1)=Sqrt(Dot(A,A))
         Cell(2)=Sqrt(Dot(B,B))
         Cell(3)=Sqrt(Dot(C,C))
         Cell(5)=-ABS(Dot(A,C)/Cell(1)/Cell(3))
         M(1) = -1.0
         M(2) = -1.0
         M(8) = 2.0
         M(6) = 1.0

      ELSEIF (itcase .eq. 30 .or. rfcase .eq. '56B') THEN
         call v2cart (v,C,A,v3)
         do 5620 I=1,3
            v3(I)=2*v3(I)
 5620    continue
         call vecdif(v3,A,B)
         Cell(1)=Sqrt(Dot(A,A))
         Cell(2)=Sqrt(Dot(B,B))
         Cell(3)=Sqrt(Dot(C,C))
         Cell(5)=-ABS(Dot(A,C)/Cell(1)/Cell(3))
         M(4) = -1.0
         M(5) = -1.0
         M(8) = 2.0
         M(3) = 1.0

      ELSEIF (itcase .eq. 29 .or. rfcase .eq. '56C') THEN
         call v2cart (v,A,v2,C)
         do 5640 I=1,3
            v2(I)=2*v2(I)
 5640    continue
         call vecdif(v2,A,B)
         Cell(1)=Sqrt(Dot(A,A))
         Cell(2)=Sqrt(Dot(B,B))
         Cell(3)=Sqrt(Dot(C,C))
         Cell(5)=-ABS(Dot(A,C)/Cell(1)/Cell(3))
         M(1) = -1.0
         M(2) = -1.0
         M(5) = 2.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 43 .or. rfcase .eq. '57A') THEN
         call v2cart(v,A,C,v3)
         do 5700 I=1,3
            v3(I)=2*v3(I)
 5700    continue
         call vecsum(v3,C,v3)
         call vecsum(v3,A,B)
         Cell(1)=Sqrt(Dot(A,A))
         Cell(2)=Sqrt(Dot(B,B))
         Cell(3)=Sqrt(Dot(C,C))
         Cell(5)=-ABS(Dot(A,C)/Cell(1)/Cell(3))
         M(1) = 1.0
         M(2) = 1.0
         M(5) = 1.0
         M(8) = 2.0
         M(6) = 1.0

      ELSEIF (itcase .eq. 17 .or. rfcase .eq. '57B') THEN
         call v2cart(v,v1,v2,v3)
         call vecsum(v1,v2,B)
         call vecsum(v3,v2,C)
         call vecsum(v3,v1,A)
         Cell(1)=Sqrt(Dot(A,A))
         Cell(2)=Sqrt(Dot(B,B))
         Cell(3)=Sqrt(Dot(C,C))
         Cell(5)=-ABS(Dot(A,C)/Cell(1)/Cell(3))
         M(4) = 1.0
         M(7) = 1.0
         M(2) = 1.0
         M(5) = 1.0
         M(3) = 1.0
         M(9) = 1.0

      ELSEIF (itcase .eq. 27 .or. rfcase .eq. '57C') THEN
         call v2cart(v,B,v2,v3)
         call vecsum(v2,v3,vtemp)
         call vecdif(vtemp,B,C)
         call vecdif(v2,v3,A)
         Cell(1)=Sqrt(Dot(A,A))
         Cell(2)=Sqrt(Dot(B,B))
         Cell(3)=Sqrt(Dot(C,C))
         Cell(5)=-ABS(Dot(A,C)/Cell(1)/Cell(3))
         M(4) = 1.0
         M(7) = -1.0
         M(2) = 1.0
         M(3) = -1.0
         M(6) = 1.0
         M(9) = 1.0
      ENDIF
      do 8000 I = 4,6
         IF (Cell(I)  .NE.  Maxint) THEN
            Cell(I) = 180.0/PI*acos(Cell(I))
            write(cellst(i),'(f8.2)') cell(i)
         ENDIF
 8000 CONTINUE
      do 8100 i = 1,3
         if(cell(i).ne.Maxint)
     *     write(cellst(i),'(f8.2)') cell(i)
 8100 CONTINUE
      do 8500 I = 4,6
         IF (altcel(I)  .NE.  Maxint) THEN
            altcel(I) = 180.0/PI*acos(altcel(I))
         ENDIF
 8500 CONTINUE
      do 8600 I = 1,6
      do 8700 II = 1,8
      if (cellst(I)(II:II).ne.' ') goto 8800
 8700 continue
      cellst(I)='.'
      return
 8800 ctemp = cellst(I)(II:8)
      cellst(I)=ctemp
 8600 continue
      return

      END

C**********************************************************************C
      subroutine  v2Cart(v ,v1,v2,v3 )

C Compute the 3-space Cartesianizing transformation matrix corresponding
C to a particular g6 vector. The base vectors of the transformation
C are returned

      real v(6),v1(3),v2(3),v3(3),cell(6)
      real mat8
C-----------------------------------------------------------------------
      Cell(1) = Root(v(1))
      Cell(2) = Root(v(2))
      Cell(3) = Root(v(3))
      Cell(4) = 0.5*v(4)/(Cell(2)*Cell(3))
      Cell(5) = 0.5*v(5)/(Cell(1)*Cell(3))
      Cell(6) = 0.5*v(6)/(Cell(1)*Cell(2))
      SinAl = Sqrt(1-Sqr(Cell(4)))
      SinBe = Sqrt(1-Sqr(Cell(5)))
      SinGa = Sqrt(1-Sqr(Cell(6)))

      v1(1) = Cell(1)
      v1(2) = 0
      v1(3) = 0

      v2(1) = Cell(2) * Cell(6)
      v2(2) = Cell(2) * SinGa
      v2(3) = 0

      v3(1) = Cell(3) * Cell(5)
      Mat8 = (Cell(4)-Cell(5)*Cell(6)) / SinGa
      v3(2) = Mat8 * Cell(3)
      v3(3) = Cell(3)*Sqrt( Sqr(SinBe)-Sqr(Mat8) )
      end


C***********************************************************************
      SUBROUTINE WRCELL (C,CE,TEST)
      include 'ITERATE.cmn'
      REAL C(6),CE(6)
      CHARACTER *6 TEST
C-----------------------------------------------------------------------
      IF (TEST .NE. 'WRCELL') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN WRCELL'
         STOP
      ENDIF

      SUM = 0.0
      DO 1000 I=1,6
 1000 SUM = SUM + CE(I)
      if (ostyle.ne.'CIF ') then
      IF (SUM .EQ. 0) THEN
         WRITE (*,'(1x,a3,1X,6F10.3)') ' '//hm//' ',C
      ELSE
         DO 3000 I=1,6
            IF (CE(I) .GT. 0.0) THEN
               WRITE (*,'(1x,a3,1X,F10.3,3X,F10.3)')
     *         ' '//hm//' ',C(I),CE(I)
            ELSE
               WRITE (*,'(1x,a3,1X,F10.3,3X,F10.3)') ' '//hm//' ',C(I)
            ENDIF
 3000    CONTINUE
      ENDIF
      WRITE (*,*)
      endif
      END


C***********************************************************************
      subroutine wrcent(lat,retcel)
      include 'ITERATE.cmn'
      character  *2 lat
      real retcel(6)
      character *79 line
C-----------------------------------------------------------------------


      if (lat(2:2) .eq. 'P') then
         line = ' primitive'
      elseif (lat(2:2) .eq. 'S') then
         line = ' side-centered '
      elseif (lat(2:2) .eq. 'I') then
         line = ' body-centered '
      elseif (lat(2:2) .eq. 'F') then
         line = ' face-centered '
      elseif (lat(2:2) .eq. 'R') then
         line = ' as rhomboh. '
      endif

      write (line(16:),'(f10.3)') retcel(1)



      if (lat(1:1) .eq. 'm') then
         write (line(26:),'(f10.3)') retcel(2)
         write (line(36:),'(f10.3)') retcel(3)
         write (line(49:),'(a)') 'beta'
         write (line(56:),'(f10.3)') retcel(5)
      elseif (lat(1:1) .eq. 'o') then
         write (line(26:),'(f10.3)') retcel(2)
         write (line(36:),'(f10.3)') retcel(3)
      elseif (lat(1:1) .eq. 't') then
         write (line(26:),'(f10.3)') retcel(3)
      elseif (lat .eq. 'hR') then
         write (line(29:),'(a)') 'alpha'
         write (line(36:),'(f10.3)') retcel(4)
      elseif (lat(1:1) .eq. 'h') then
         write (line(26:),'(f10.3)') retcel(3)
      elseif (lat(1:1) .eq. 'c') then
      endif
      if (ostyle.ne.'CIF ') write (*,'(a)') '  '//hm//' '//line
      end

C***********************************************************************
      SUBROUTINE WRMATR (N,A)
      include 'ITERATE.cmn'
      REAL A(N,N)
C-----------------------------------------------------------------------
      DO 1000 I=1,N
         WRITE (*,*) ' '//hm,(A(J,I),J=1,N)
 1000 CONTINUE
      END

C***********************************************************************
      SUBROUTINE WRSIZE (SIZE,ERRSIZ,RATIO,TEST)
      include 'ITERATE.cmn'
      CHARACTER *6 TEST
C-----------------------------------------------------------------------
      IF (TEST .NE. 'WRSIZE') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN WRSIZE'
         STOP
      ENDIF
      if (ostyle.ne.'CIF ') WRITE
     2  (*,'(''  '//hm//
     3  ' INPUT VECTOR SIZE AND ERROR AND RATIO '',3F10.2)')
     4  SIZE,ERRSIZ,RATIO
      END

C***********************************************************************
      SUBROUTINE WRVEC6(V,VE,TEST)
      include 'ITERATE.cmn'
      REAL V(6),VE(6)
      CHARACTER *6 TEST
C-----------------------------------------------------------------------

      IF (TEST .NE. 'WRVEC6') THEN
         WRITE (*,*) ' '//hm//' TEST WAS WRONG IN WRVEC6'
         STOP
      ENDIF

      SUM = 0.0
      DO 1000 I=1,6
 1000 SUM = SUM + ABS(VE(I))
      if (ostyle.ne.'CIF ') then
      IF (SUM .EQ. 0) THEN
         WRITE (*,'(1x,a3,1X,6F10.2)') ' '//hm//' ',V
      ELSE
         DO 3000 I=1,6
            WRITE (*,'(1x,a3,1X,F10.3,3X,F10.3)') ' '//hm//' '
     2        ,V(I),VE(I)
 3000    CONTINUE
      ENDIF
      WRITE (*,*)
      endif
      END

C**********************************************************************C
      SUBROUTINE ZEROS (N,V)
      REAL V(N)
C----------------------------------------------------------------------C
      DO 1000 I=1,N
 1000 V(I) = 0.0
      END

C**********************************************************************C
      SUBROUTINE INVER (A,B)
C----INVERT A THREE BY THREE MATRIX
      REAL X(9)
      DIMENSION A(9),B(9)
      DIMENSION IDATA1(3),IDATA2(3)
      DATA IDATA1 /4,7,1/
      DATA IDATA2 /7,1,4/
C----------------------------------------------------------------------C
      J = 0
      DO 1000 I=1,9,3
      J = J + 1
      ID1 = IDATA1(J)
      ID2 = IDATA2(J)
      CALL CROSS(A(ID1),A(ID2),X(I))
 1000 CONTINUE

      DETA = DET(A)
      IF (ABS(DETA) .LE. 1.0E-20) THEN
         DETA = SIGN(1.0E-20,DETA)
      ELSE
         DETA = 1.0 / DETA
      ENDIF

      CALL CONMAT (X,DETA,X)
      CALL TRNSPZ(X,B)
      END

C**********************************************************************C
      SUBROUTINE VECSUM (X,Y,Z)
C----ADD TWO VECTORS AND RETURN THE SUM IN Z
      DIMENSION X(3), Y(3), Z(3)
C----------------------------------------------------------------------C
      DO 1000 I=1,3
      Z(I) = X(I) + Y(I)
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE VECDIF (X,Y,Z)
C----SUBTRACT TWO VECTORS AND RETURN THE RESULT IN Z
      DIMENSION X(3), Y(3), Z(3)
C----------------------------------------------------------------------C
      DO 1000 I=1,3
      Z(I) = X(I) - Y(I)
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE CPYVEC (X,Y)
C----COPY A VECTOR X INTO A VECTOR Y
      DIMENSION X(3), Y(3)
C----------------------------------------------------------------------C
      DO 1000 I=1,3
      Y(I) = X(I)
 1000 CONTINUE
      END

C**********************************************************************C
      FUNCTION DOT (X,Y)
C----COMPUTE AND RETURN THE DOT PRODUCT OF X AND Y
      DIMENSION X(3),Y(3)
C----------------------------------------------------------------------C
      DOT = 0.0
      DO 1000 I=1,3
      DOT = DOT + X(I) * Y(I)
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE TRNSPZ (A,B)
C----PUT THE TRANSPOSE OF A INTO B
      DIMENSION A(9), B(9)
C----------------------------------------------------------------------C
      J = 0
      DO 1000 I=1,9,3
      J = J + 1
      CALL UNVEC (A(I),B(J),B(J+3),B(J+6))
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE MATMUL (A,B,C)
C----MULTIPLY TWO MATRICIES
      REAL X(9)
      DIMENSION A(9), B(9), C(9)
C----GET THE TRANSPOSE OF B INTO X
C----------------------------------------------------------------------C
      CALL TRNSPZ (B,X)
      IJ = 0
      DO 2000 I=1,9,3
      DO 2000 J=1,9,3
      IJ = IJ + 1
      C(IJ) = DOT(A(I),X(J))
 2000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE CROSS (X,Y,Z)
C----COMPUTE Z = X CROSS Y
      DIMENSION X(3),Y(3),Z(3)
C----------------------------------------------------------------------C
      Z(1) = X(2)*Y(3) - Y(2)*X(3)
      Z(2) =-X(1)*Y(3) + Y(1)*X(3)
      Z(3) = X(1)*Y(2) - Y(1)*X(2)
      END

C**********************************************************************C
      FUNCTION DET(A)
C----RETURN THE VALUE OF THE DETERMINANT OF A MATRIX
      REAL X(3)
      DIMENSION A(9)
C----------------------------------------------------------------------C
      CALL CROSS (A(1),A(4),X)
      DET = DOT (X,A(7))
      END

C**********************************************************************C
      SUBROUTINE CONMAT(AMAT,X,BMAT)
      DIMENSION AMAT(9),BMAT(9)
C----------------------------------------------------------------------C
      DO 1000 I=1,9
      BMAT(I) = X * AMAT(I)
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE UNVEC (X,F,G,H)
C----RETURN THE VECTOR COMPONENTS AS SCALARS
      DIMENSION X(3)
C----------------------------------------------------------------------C
      F = X(1)
      G = X(2)
      H = X(3)
      END

C**********************************************************************C
      SUBROUTINE WRCLV6 (C,CE,V,VE,TEST)
      include 'ITERATE.cmn'
      logical pnumb_
      CHARACTER *6 TEST
      REAL C(6),CE(6),V(6),VE(6)
C----------------------------------------------------------------------C
      if (test .ne. 'WRCLV6') then
         write (*,*) ' '//hm//
     2     ' test string was not WRCLV6 in that routine'
         stop
      endif
      if (ostyle.ne.'CIF ') then

         DO 3000 I=1,6
            WRITE (*,'(1x,a3,1X,F10.3,3X,F8.3,8X,f10.3,3x,f8.2))')
     2         ' '//hm//' ',C(I),CE(I),V(I),VE(I)
 3000    CONTINUE

      WRITE (*,*)
      else
      do ii = 1,6
        cifres = pnumb_(' ',c(ii),0.)
      enddo
      endif
      END
