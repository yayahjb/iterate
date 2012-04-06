C     Bravais General Analysis of Lattices (BGAOL)
C     The Program formerly known as ITERATE
C
C     Lawrence C Andrews[1] and Herbert J. Bernstein[2,*]
C
C     [1] Micro Encoder Inc., 11533 NE 118th St, #200,
C         Kirkland, WA 98034-7111 USA
C     [2] Dowling College, 1300 William Floyd Parkway,
C         Shirley, NY 11967 USA
C     [*] To whom correspondence should be addressed.
C         Email: yaya@dowling.edu
C
C     Copyright 1996, 2012, all rights reserved
C
C*******************************************************
C    You may redistribute this program under the terms
C    of the GPL.
C
C    ALternatively you may redistribute this functions
C    and subroutines of this program as an API under the
C    terms of the LGPL
C*******************************************************
C*************************** GPL NOTICES ******************************
C*                                                                    *
C* This program is free software; you can redistribute it and/or      *
C* modify it under the terms of the GNU General Public License as     *
C* published by the Free Software Foundation; either version 2 of     *
C* (the License, or (at your option) any later version.               *
C*                                                                    *
C* This program is distributed in the hope that it will be useful,    *
C* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      *
C* GNU General Public License for more details.                       *
C*                                                                    *
C* You should have received a copy of the GNU General Public License  *
C* along with this program; if not, write to the Free Software        *
C* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA           *
C* 02111-1307  USA                                                    *
C*                                                                    *
C**********************************************************************/

C************************* LGPL NOTICES *******************************
C*                                                                    *
C* This library is free software; you can redistribute it and/or      *
C* modify it under the terms of the GNU Lesser General Public         *
C* License as published by the Free Software Foundation; either       *
C* version 2.1 of the License, or (at your option) any later version. *
C*                                                                    *
C* This library is distributed in the hope that it will be useful,    *
C* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C* Lesser General Public License for more details.                    *
C*                                                                    *
C* You should have received a copy of the GNU Lesser General Public   *
C* License along with this library; if not, write to the Free         *
C* Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
C* MA  02110-1301  USA                                                *
C*                                                                    *
C**********************************************************************/

C
C In simple terms, what this program does is to find
C the cells which are "close" to the cell given, in
C order to help find the Bravais lattice of highest
C symmetry consistent with the cell.
C
C A central problem in the solution of every crystal
C structure is to determine the correct Bravais lattice
C of the crystal. The Bravais lattices as they are
C usually listed are:
C
C aP  triclinic (anorthic) primitive
C mP  monoclinic primitive
C mS  monoclinic side-centered (usually C-centered)
C oP  orthorhombic primitive
C oS  orthorhombic side-centered
C oF  orthorhombic face-centered
C oI  orthorhombic body-centered
C hP  hexagonal primitive
C hR  hexagonal rhombohedrally-centered
C tP  tetragonal primitive
C tI  tetragonal body-centered
C cP  cubic primitive
C cF  cubic face-centered
C cI  cubic body-centered
C
C Failure to find the highest correct symmetry has several
C consequences, the worst of which is that the structure
C may not be solved. The least of the consequences is that
C some successor to Richard Marsh may publish a paper that
C points out the error, corrects it, and finds a better
C solution to the structure. Many methods have been
C described for finding the correct Bravais lattice. A
C summary of the published methods was published in the
C paper that described the G6 formalism (which is used
C in the program on this web page).
C
C "Lattices and Reduced Cells as Points in 6-Space and
C Selection of Bravais Lattice Type by Projections."
C Lawrence C. Andrews and Herbert J. Bernstein, Acta
C Crystallographica, A44, 1009-1018 (1988).
C
C The program on BGAOL implements a search in G6 for the
C various Bravais lattices that the user's cell may fit.
C For each lattice type, the best metric match is reported.
C If the higher symmetry type is actually correct, then
C that is likely to be the best cell from which to start
C further refinement. However, the possibility exists
C that one of the rejected cells (which did not match as
C well) was actually the correct one to use. The reason
C for this ambiguity is experimental error and its
C propagation in the transformations of the lattices
C in the program. Fortunately, the rejected cells are
C usually quite similar to the accepted one.
C
C A note on standard deviations: First, even in the best
C of circumstances, standard deviations of unit cell
C dimensions from 4-circle diffractometer data are always
C underestimated (by at least a factor of 2). In addition,
C the points chosen for the determination are often not
C well distributed (for example all in the first octant
C of orthorhombic lattices). These less than optimal
C choices cause substantial systematic error. The
C experimental errors are amplified in the mathematical
C conversions between various lattices that any lattice
C search program must perform. It is not a rare occurrence
C for angles to be incorrect by 0.5 degrees in initial
C unit cell determinations.
C
C Note: Even in most well determined unit cells, the actual
C errors in the edge lengths are 0.2 to 0.5 parts per
C thousand. (Note that reproducibility of the measurements
C is substantially better, leading to the illusion that
C diffractometers produce excellent unit cell parameters).
C Use of standard deviations that are too small is a
C common reason for failure of Bravais lattice searches.
C For small molecules, 0.1 Angstroms is a reasonable error
C for the edge lengths, for proteins, 0.4 to 0.5 (or even
C more for preliminary measurements). Accurate unit cell
C parameters must by determined by a number of more complex
C methods and must include extrapolation to remove systematic
C effects. For an excellent summary, see "Xray Structure
C Determination", G.H.Stout and L.H.Jensen, Wiley, 1989.
C
C Note on the name BGAOL -- gaols have lots of cells
C
      include "mkrefl.for"
      include "near6.for"
      include "e3tog6.for"

C**********************************************************************C
      PROGRAM RED
      implicit none

      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      INTEGER MAXPRJ
      PARAMETER (MAXPRJ=42)
      INTEGER ITDESG(MAXPRJ)
      INTEGER IORD(MAXPRJ)
      LOGICAL STAND(MAXPRJ)
      CHARACTER *2 CHRLAT(MAXPRJ)
      real*8 PJNORM(MAXPRJ)
      INTEGER PRJ(36,MAXPRJ)
      real*8 P(36),AP(36)
      real*8 CV(6),CE(6),G(6),GE(6),TG(6),AG(6),COUT(6)
      real*8 VRED(6,MAXPRJ)
      real*8 VOUT(6,MAXPRJ)
      real*8 MPRIM(36)
      INTEGER MRED(36),MRTEST(36),IUM(36),MTEMP(36)
      real*8 M3RED(3,3)
      real*8 GOUT(6),GRED(6),CRED(6),VTEST(6)
      real*8 ZEROS(6)
      LOGICAL INPCEL
      EXTERNAL INPCEL
      CHARACTER LATSYM
      integer nvmax, mxtree
      integer iv, nv
      real*8 band
      integer bandcount
      integer INRMSQ
      logical DEBUG, DEBUGD

      real*8 dbest, dcurr

      PARAMETER (NVMAX=100000)
      PARAMETER (MXTREE=11*NVMAX)
      real*8 TREE(MXTREE)
      real*8 V(6,NVMAX)
      real*8 VBEST(6),AVBEST(6)
      real*8 NewCell(6)
      real*8 CenteredCell(6)

      real*8 DOTVN
      EXTERNAL DOTVN
      INTEGER NPROJ,I,J
      real*8 RATIO,SIZE,ERRSIZ
      real*8 anorm

      REAL*8 BEST(42)
      real*8 distances2Boundaries(21)
      character*7 bnames(21)
      integer ib
      character*3 sprefix

      CHARACTER *100 TITLE

      DATA ZEROS /6*0.0/
      DATA bnames/ '1','2','3','4','5','6','7','8','9',
     *  'A','B','C','D','E','F',
     *  '678X','9ABCDEX','FX','67','9A','CD'/
C----------------------------------------------------------------------C

      iunit0 = 0
      iunit1 = 1
      iunit2 = 2
      iunit3 = 3
      iunit10 = 10
      DEBUG = .true.
      DEBUGD = .false.
      call IUNTMN(6,IUM)
      write(*,*) ("*+",i=1,35)

      CALL BLDPRJ (MAXPRJ,NPROJ,ITDESG,CHRLAT,PJNORM,PRJ,'BLDPRJ')
      write (*,*) ' nproj = ',nproj

      IF (INPCEL(LATSYM,TITLE,CV,CE)) THEN
         CALL CTOG6(CV,CE,G,GE,SIZE,ERRSIZ,RATIO,'CTOG6 ')
         WRITE (*,'('' INPUT CELL AND ERRORS '')')
         CALL WRCELL(CV,CE,'WRCELL')
         WRITE (*,*)
         WRITE (*,*) ' INPUT VECTOR AND ERRORS'
         CALL WRVEC6(G,GE,'WRVEC6')
         WRITE (*,*)
         CALL WRSIZE (SIZE,ERRSIZ,RATIO,'WRSIZE')
         CALL MKPRIM (LATSYM,G,MPRIM,GOUT,'MKPRIM')
         WRITE (*,*) ' AFTER MKPRIM'
         CALL dg6toe3(MPRIM,M3RED)
         WRITE (*,*) "g6 Centering Matrix:"
         WRITE (*,'(12x,6d10.3)') MPRIM
         WRITE (*,*) "E3 Centering Matrix:"
         WRITE (*,'(12x,3f10.3)') ((M3RED(I,J),i=1,3),j=1,3)

         CALL CHKVEC(GOUT)
         CALL IUNTMN(6,MRED)
         CALL REDUCE (GOUT,MRED,GRED,0.D0,'REDUCE')
         CALL CHKVEC(GRED)
         CALL G6TOC (GRED,CRED,'G6TOC ')
         SIZE = DSQRT(DOTVN(6,GRED,GRED))

         ERRSIZ = RATIO * SIZE
         WRITE (*,*)
         WRITE (*,*) TITLE
         WRITE (*,'('' REDUCED CELL  '')')
         WRITE (*,'('' Red. Cell   '',6F11.4)')CRED
         WRITE (*,'('' Red. Vector '',6F10.5)') GRED
         WRITE (*,*) "G6 Reduction Matrix:"
         WRITE (*,'(12x,6I6)') MRED 
         CALL g6toe3(MRED,M3RED)
         WRITE (*,*) "E3 Reduction Matrix:"
         WRITE (*,'(12x,3f10.3)') ((M3RED(I,J),i=1,3),j=1,3)
         

      CALL MKREFL (RATIO,MXTREE,TREE,NVMAX,V,NV,GRED,'MKREFL')





         DO 4000 I=1,NPROJ
            BEST(I) = -19191
            STAND(I) = .false.
            DBEST = 1.0E20
            DO 1000 J=1,36
               P(J) = PRJ(J,I)/PJNORM(I)
 1000          AP(J) = -P(J)
            DO 2000 J=1,36,7
               AP(J) = 1.0D0 + AP(J)
 2000       CONTINUE
            DO 3000  IV=1,NV
            CALL REDUCE (V(1,IV),MRTEST,VTEST,
     2           ERRSIZ/2.D0,'REDUCE')
            CALL vdifn(6,V(1,IV),VTEST,VTEST)
            CALL RMV6 (V(1,IV),AP,AG)
            CALL RMV6 (V(1,IV),P,TG)
            DCURR = DOTVN(6,AG,AG)
            IF ( ((DCURR .LT. DBEST-1.D-10) 
     2        .OR. ((.NOT.STAND(I)).AND.
     3        anorm(6,VTEST).le.ERRSIZ/2.D0)
     4        .AND. SQRT(DCURR).LE.ERRSIZ/2.D0)
     5        .AND. (TG(1).GT.1.0D0 .AND. TG(2).GT.1.0D0
     2        .AND. TG(3).GT.1.0D0
     3        .AND. (
     4           (TG(4)*TG(5)*TG(6) .LE. 1.E-4 
     5            .AND. ITDESG(I).LT.0)
     6           .OR. (TG(4)*TG(5)*TG(6) .GE. -1.E-4
     7            .AND. ITDESG(I).GT.0)))) THEN
               DBEST = DOTVN(6,AG,AG)
               CALL CPYVN(6,TG,VBEST)
               CALL CPYVN(6,AG,AVBEST)
               CALL CPYVN(6,V(1,IV),VRED(1,I))
               IF(anorm(6,VTEST).gt.ERRSIZ/2.D0) THEN
                 STAND(I) = .false.
               ELSE
                 STAND(I) = .true.
               ENDIF
            ENDIF

 3000       CONTINUE
            IF (DSQRT(DOTVN(6,AVBEST,AVBEST)) .GT. ERRSIZ*3.5D0) 
     *        GO TO 4000
            IF (DSQRT(DBEST) .LE. 999.0) THEN
               sprefix = "NON"
               if (STAND(I)) sprefix = "   "
               WRITE (*,*)
               IF (ITDESG(I).GT.0) THEN
               WRITE (*,'(1X,I3,2X,''IT('',I2'')'',
     2            4H +++,1X,A3,A,2X,A2,2X,A100)')
     3            I,IABS(ITDESG(I)),
     4            sprefix,"STANDARD ",CHRLAT(I), TITLE
               ELSE
               WRITE (*,'(1X,I3,2X,''IT('',I2'')'',
     2            4H ---,1X,A3,A,2X,A2,2X,A100)')
     3            I,IABS(ITDESG(I)),
     4            sprefix,"STANDARD ",CHRLAT(I), TITLE
               ENDIF
               WRITE (*,'(6X,F5.1,1X,F5.2,3H V=,6F8.1)')
     2             DSQRT(DOTVN(6,AVBEST,AVBEST)),
     3             DSQRT(DOTVN(6,AVBEST,AVBEST))/ERRSIZ,VBEST
               CALL G6TOC (VBEST,COUT,'G6TOC ')
               WRITE (*,'(6X,''cell'',8X,2(3F8.2,'' ''))') COUT

               BEST(I) = DSQRT(DOTVN(6,AVBEST,AVBEST))
               call cpyvn(6,TG,VOUT(1,I))
            ELSE
               BEST(I) = -19191
               STAND(I) = .false.
            ENDIF
            WRITE (*,*)
            call UncenterVector( VBEST, ITDESG(I), CHRLAT(I), NewCell )
            write(*,"(8H newcell, 6F9.3)" ) NewCell
            CenteredCell(1) = DSQRT( NewCell(1) );
            CenteredCell(2) = DSQRT( NewCell(2) );
            CenteredCell(3) = DSQRT( NewCell(3) );
            CenteredCell(4) =
     *          DACOS(0.5*NewCell(4)/CenteredCell(2)/CenteredCell(3)) *
     *          180.0D0/DATAN(1.0D0)/4.0D0
            CenteredCell(5) =
     *          DACOS(0.5*NewCell(5)/CenteredCell(1)/CenteredCell(3)) *
     *          180.0D0/DATAN(1.0D0)/4.0D0
            CenteredCell(6) =
     *          DACOS(0.5*NewCell(6)/CenteredCell(1)/CenteredCell(2)) *
     *          180.0D0/DATAN(1.0D0)/4.0D0
            write(*,"(8H **cell*, 6F9.3)" ) CenteredCell

         IF (DEBUGD) THEN
            WRITE (*,'(1X,6F10.6)') P
            WRITE (*,*)
         call BoundaryDistances( 1, VBEST,
     2                   distances2Boundaries, ERRSIZ, DEBUG )
         do 3500 ib=1,21
            write(*,*)
     2            " returned boundary ",bnames(ib),
     3            distances2Boundaries(ib)
 3500    continue
         ENDIF
         WRITE (*,*)
         

 4000    CONTINUE

      call SortDists(NPROJ, BEST, IORD, ERRSIZ)
      band = ERRSIZ/2
      bandcount = 0
      write (*,*) "Final List"
      DO iv = 1,2
      if (iv.eq.1) write(*,*)"STANDARD PRESENTATIONS"
      if (iv.eq.2) write(*,*)"NONSTANDARD PRESENTATIONS"
      DO 5000 i=1,NPROJ
         IF ( BEST(IORD(I)) .ge. 0 
     *     .and. ((STAND(IORD(I)).and.iv.eq.1)
     *        .or. ((.NOT.STAND(IORD(I))).and.iv.eq.2))) THEN
           IF ( BEST(IORD(I)) .gt. band) THEN
           IF (bandcount .gt. 0) THEN
             WRITE(*, "(1X,A)"),"-----------------------------"
           ENDIF
           band = ERRSIZ/2.D0
     *       *dble(INT(2.0*BEST(IORD(I))/ERRSIZ)+1)
           bandcount = 0
           if (band .gt. 3.25D0*ERRSIZ) go to 6000
           ENDIF
           sprefix = "NON"
           if (STAND(IORD(I))) sprefix = "   "   
           if (ITDESG(IORD(I)).GT.0) then
           WRITE(*, "(4H IT(, I3, 2H) ,4H +++,1X,A3,
     2        9HSTANDARD ,A2, 2x, F9.3)" )
     3        IABS(ITDESG(IORD(I))), sprefix,
     4        CHRLAT(IORD(I)), BEST(IORD(I))
           else
           WRITE(*, "(4H IT(, I3, 2H) ,4H ---,1x,A3,
     2        9HSTANDARD ,A2, 2x, F9.3)" )
     3        IABS(ITDESG(IORD(I))), sprefix,
     4        CHRLAT(IORD(I)), BEST(IORD(I))
           endif
           bandcount = bandcount+1
         ENDIF
 5000 CONTINUE
      enddo

 6000 CONTINUE
      ENDIF
      END

C**********************************************************************C
      subroutine UncenterVector( VBEST, ITDESG, CHRLAT, NewCell )
      implicit none
      real*8 VBEST(6), NewCell(6)
      integer ITDESG
      integer iCellType
      Character*2 CHRLAT

      real*8 matrices(36,44)
      integer i

      character*1 latticeSymbol

      data (matrices(i,01),i=1,36)/
     *   1.00,  1.00,  1.00, -1.00,  1.00, -1.00,
     *   1.00,  1.00,  1.00, -1.00, -1.00,  1.00,
     *   1.00,  1.00,  1.00,  1.00, -1.00, -1.00,
     *  -2.00,  2.00, -2.00,  0.00,  2.00,  0.00,
     *  -2.00, -2.00,  2.00,  0.00,  0.00,  2.00,
     *   2.00, -2.00, -2.00,  2.00,  0.00,  0.00
     * /
      data (matrices(i,02),i=1,36)/
     *   1.00,  1.00,  0.00, -0.00,  0.00, -1.00,
     *   1.00,  0.00,  1.00,  0.00, -1.00, -0.00,
     *   1.00,  1.00,  1.00,  1.00,  1.00,  1.00,
     *   2.00, -0.00, -2.00, -1.00,  0.00,  1.00,
     *  -2.00,  2.00, -0.00,  1.00, -1.00,  0.00,

     *  -2.00, -0.00,  0.00, -1.00,  1.00,  1.00
     * /
      data (matrices(i,03),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,04),i=1,36)/
     *   1.00,  1.00,  0.00, -0.00,  0.00, -1.00,
     *   1.00,  0.00,  1.00,  0.00, -1.00, -0.00,
     *   1.00,  1.00,  1.00,  1.00,  1.00,  1.00,
     *   2.00, -0.00, -2.00, -1.00,  0.00,  1.00,
     *  -2.00,  2.00, -0.00,  1.00, -1.00,  0.00,
     *  -2.00, -0.00,  0.00, -1.00,  1.00,  1.00
     * /
      data (matrices(i,05),i=1,36)/
     *   1.00,  0.00,  1.00,  0.00,  1.00,  0.00,
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   0.00,  2.00,  0.00,  1.00,  1.00,  1.00,
     *   0.00,  0.00,  2.00,  1.00,  1.00,  1.00,
     *   2.00,  0.00,  0.00,  1.00,  1.00,  1.00
     * /
      data (matrices(i,06),i=1,36)/
     *   1.00,  0.00,  1.00,  0.00,  1.00,  0.00,
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  2.00,  0.00,  1.00,  1.00,  1.00,
     *   2.00,  0.00,  0.00,  1.00,  1.00,  1.00,
     *   0.00,  0.00,  2.00,  1.00,  1.00,  1.00
     * /
      data (matrices(i,07),i=1,36)/
     *   1.00,  0.00,  1.00,  0.00,  1.00,  0.00,
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  2.00,  0.00,  1.00,  1.00,  1.00,
     *   2.00,  0.00,  0.00,  1.00,  1.00,  1.00,
     *   0.00,  0.00,  2.00,  1.00,  1.00,  1.00
     * /
      data (matrices(i,08),i=1,36)/
     *   1.00,  0.00,  1.00,  0.00,  1.00,  0.00,
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   0.00,  2.00,  0.00,  1.00,  1.00,  1.00,
     *   0.00,  0.00,  2.00,  1.00,  1.00,  1.00,
     *   2.00,  0.00,  0.00,  1.00,  1.00,  1.00
     * /
      data (matrices(i,09),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  1.00,  0.00,  0.00, -0.00, -1.00,
     *   1.00,  1.00,  9.00, -3.00, -3.00,  1.00,
     *   2.00, -2.00,  0.00,  3.00, -3.00,  0.00,
     *  -2.00, -0.00,  0.00,  0.00,  3.00, -1.00,
     *  -2.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,10),i=1,36)/
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   1.00,  1.00,  0.00,  0.00, -0.00, -1.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *  -0.00,  0.00,  0.00,  1.00, -1.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  1.00,  0.00,
     *  -2.00,  2.00,  0.00,  0.00,  0.00,  0.00
     * /
      data (matrices(i,11),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,12),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,13),i=1,36)/
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   1.00,  1.00,  0.00,  0.00, -0.00, -1.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *  -0.00,  0.00,  0.00,  1.00, -1.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  1.00,  0.00,
     *  -2.00,  2.00,  0.00,  0.00,  0.00,  0.00
     * /
      data (matrices(i,14),i=1,36)/
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   1.00,  1.00,  0.00,  0.00, -0.00, -1.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *  -0.00,  0.00,  0.00,  1.00, -1.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  1.00,  0.00,
     *  -2.00,  2.00,  0.00,  0.00,  0.00,  0.00
     * /
      data (matrices(i,15),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  1.00,  4.00,  2.00,  2.00,  1.00,
     *   0.00,  2.00,  0.00,  2.00,  0.00,  1.00,
     *   2.00,  0.00,  0.00,  0.00,  2.00,  1.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,16),i=1,36)/
     *   1.00,  1.00,  0.00, -0.00,  0.00, -1.00,
     *   1.00,  1.00,  4.00,  2.00,  2.00,  1.00,
     *   1.00,  1.00,  0.00, -0.00, -0.00,  1.00,
     *  -2.00, -2.00,  0.00, -2.00, -2.00, -2.00,
     *  -2.00,  2.00,  0.00, -0.00,  0.00,  0.00,
     *   2.00, -2.00,  0.00, -2.00,  2.00,  0.00
     * /
      data (matrices(i,17),i=1,36)/
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   1.00,  1.00,  0.00,  0.00,  0.00,  1.00,
     *   1.00,  0.00,  1.00,  0.00,  1.00,  0.00,
     *   2.00,  0.00,  0.00,  1.00,  1.00,  1.00,
     *   0.00,  0.00,  2.00,  1.00,  1.00,  1.00,
     *   0.00,  2.00,  0.00,  1.00,  1.00,  1.00
     * /
      data (matrices(i,18),i=1,36)/
     *   0.00,  1.00,  1.00, -1.00,  0.00, -0.00,
     *   1.00,  1.00,  1.00,  1.00, -1.00, -1.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   2.00, -0.00, -0.00, -0.00, -1.00, -1.00,
     *   0.00, -0.00,  0.00,  0.00,  1.00, -1.00,
     *   0.00,  2.00, -2.00,  0.00,  1.00, -1.00
     * /
      data (matrices(i,19),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00, -0.00, -0.00,
     *   1.00,  1.00,  1.00,  1.00, -1.00, -1.00,
     *   0.00,  1.00,  1.00, -1.00,  0.00, -0.00,
     *  -0.00, -2.00,  2.00,  0.00, -1.00,  1.00,
     *  -0.00, -0.00,  0.00,  0.00, -1.00,  1.00,
     *   2.00,  0.00,  0.00,  0.00, -1.00, -1.00
     * /
      data (matrices(i,20),i=1,36)/
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   0.00,  1.00,  1.00, -1.00,  0.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00, -0.00,  0.00,  0.00,  1.00, -1.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  1.00,
     *   0.00, -2.00,  2.00,  0.00,  0.00,  0.00
     * /
      data (matrices(i,21),i=1,36)/
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00
     * /
      data (matrices(i,22),i=1,36)/
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00
     * /
      data (matrices(i,23),i=1,36)/
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   0.00,  1.00,  1.00, -1.00,  0.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00, -0.00,  0.00,  0.00,  1.00, -1.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  1.00,
     *   0.00, -2.00,  2.00,  0.00,  0.00,  0.00
     * /
      data (matrices(i,24),i=1,36)/
     *   1.00,  4.00,  1.00,  2.00,  1.00,  2.00,
     *   0.00,  1.00,  1.00, -1.00,  0.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00, -0.00,  0.00,  0.00,  1.00, -1.00,
     *   2.00,  0.00,  0.00,  0.00,  1.00,  2.00,
     *   0.00, -4.00,  2.00,  1.00,  1.00, -1.00
     * /
      data (matrices(i,25),i=1,36)/
     *   0.00,  1.00,  1.00,  1.00,  0.00,  0.00,
     *   0.00,  1.00,  1.00, -1.00,  0.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00, -0.00,  0.00,  0.00,  1.00, -1.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  1.00,
     *   0.00, -2.00,  2.00,  0.00,  0.00,  0.00
     * /
      data (matrices(i,26),i=1,36)/
     *   1.00,  4.00,  0.00,  0.00, -0.00, -2.00,
     *   1.00,  0.00,  4.00,  0.00, -2.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *  -2.00,  0.00,  0.00,  0.00,  2.00,  0.00,
     *  -2.00,  0.00,  0.00,  0.00,  0.00,  2.00,
     *   2.00,  0.00,  0.00,  4.00, -2.00, -2.00
     * /
      data (matrices(i,27),i=1,36)/
     *   0.00,  1.00,  1.00, -1.00, -0.00,  0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  1.00,  1.00,  1.00, -1.00, -1.00,
     *  -2.00,  0.00,  0.00,  0.00,  1.00,  1.00,
     *  -0.00,  2.00, -2.00,  0.00,  1.00, -1.00,
     *   0.00,  0.00, -0.00,  0.00, -1.00,  1.00
     * /
      data (matrices(i,28),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00, -0.00, -0.00,
     *   1.00,  0.00,  4.00,  0.00, -2.00, -0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *  -0.00,  0.00,  0.00,  2.00,  0.00, -1.00,
     *  -0.00,  0.00,  0.00,  0.00,  0.00, -1.00,
     *   2.00,  0.00,  0.00,  0.00, -2.00, -0.00
     * /
      data (matrices(i,29),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00, -0.00, -0.00,
     *   1.00,  4.00,  0.00,  0.00, -0.00, -2.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *  -0.00,  0.00,  0.00,  2.00, -1.00,  0.00,
     *  -0.00,  0.00,  0.00,  0.00, -1.00,  0.00,
     *   2.00,  0.00,  0.00,  0.00, -0.00, -2.00
     * /
      data (matrices(i,30),i=1,36)/
     *   0.00,  1.00,  0.00, -0.00,  0.00, -0.00,
     *   0.00,  1.00,  4.00, -2.00,  0.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00, -0.00,  0.00,  0.00,  2.00, -1.00,
     *   0.00, -0.00,  0.00,  0.00,  0.00, -1.00,
     *   0.00,  2.00,  0.00, -2.00,  0.00, -0.00
     * /
      data (matrices(i,31),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,32),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,33),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,34),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00
     * /
      data (matrices(i,35),i=1,36)/
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,36),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  0.00,  4.00, -0.00,  2.00, -0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *  -0.00,  0.00, -0.00, -2.00, -0.00, -1.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *  -2.00,  0.00, -0.00,  0.00, -2.00,  0.00
     * /
      data (matrices(i,37),i=1,36)/
     *   1.00,  0.00,  4.00,  0.00,  2.00,  0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  0.00,  0.00,  2.00,  0.00,  1.00,
     *   2.00,  0.00,  0.00,  0.00,  2.00,  0.00
     * /
      data (matrices(i,38),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  4.00,  0.00, -0.00, -0.00,  2.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *  -0.00, -0.00,  0.00, -2.00, -1.00, -0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *  -2.00, -0.00,  0.00,  0.00,  0.00, -2.00
     * /
      data (matrices(i,39),i=1,36)/
     *   1.00,  4.00,  0.00,  0.00,  0.00,  2.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  2.00,  1.00,  0.00,
     *   2.00,  0.00,  0.00,  0.00,  0.00,  2.00
     * /
      data (matrices(i,40),i=1,36)/
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  4.00,  2.00, -0.00, -0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00, -0.00, -0.00, -0.00, -2.00, -1.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00, -2.00, -0.00, -2.00,  0.00,  0.00
     * /
      data (matrices(i,41),i=1,36)/
     *   0.00,  1.00,  4.00,  2.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   0.00,  0.00,  0.00,  0.00,  2.00,  1.00,
     *   0.00,  2.00,  0.00,  2.00,  0.00,  0.00
     * /
      data (matrices(i,42),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  1.00,  4.00,  2.00,  2.00,  1.00,
     *  -0.00, -2.00, -0.00, -2.00, -0.00, -1.00,
     *  -2.00, -0.00, -0.00, -0.00, -2.00, -1.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /
      data (matrices(i,43),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   1.00,  1.00,  4.00,  2.00,  2.00,  1.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  2.00,  0.00,  2.00,  0.00,  1.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00,
     *   2.00,  0.00,  0.00,  0.00,  2.00,  1.00
     * /
      data (matrices(i,44),i=1,36)/
     *   1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  1.00,  0.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  1.00,  0.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  1.00,  0.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
     *   0.00,  0.00,  0.00,  0.00,  0.00,  1.00
     * /

      iCellType = iabs(ITDESG)

      latticeSymbol(1:1) = CHRLAT(2:2)

      if (iCellType .LE. 43 ) then
         call RMV6( VBEST, matrices(1,iCellType), NewCell )
      elseif (iCellType .EQ. 15 ) then
         call RMV6( VBEST, matrices(1,iCellType), NewCell )
      else
         do i=1,6
            NewCell(i) = VBEST(i)
         enddo
      endif

      end


C**********************************************************************C
      subroutine SortDists(NDIST, DIST, IORD, ERRSIZ )
      implicit none
      integer NDIST, IORD(NDIST), II, IGAP, ITEMP
      real*8 ERRSIZ
      logical INORD
      real*8 DIST(NDIST)
      do II = 1, NDIST
        IORD(II) = II
      enddo
      if (NDIST .eq. 1) return

C     bring all elements .le. ERRSIZ/2 to the front
 100  continue
      INORD = .true.
      do ii = 1, NDIST-1
      if (DIST(IORD(II)) .gt. ERRSIZ/2
     *  .AND. DIST(IORD(II+1)).le.ERRSIZ/2) then
        ITEMP = IORD(II)
        IORD(II) = IORD(II+1)
        IORD(II+1) = ITEMP
        INORD = .false.
      endif
      enddo
      if (.not. INORD) go to 100

C     sort all elements .gt. ERRSIZ/2
      IGAP = NDIST+1
 1000 IGAP = IGAP/2
      if (IGAP .le. 0) return
 2000 continue
      INORD = .true.
      do II = 1, NDIST-IGAP, IGAP
      if (DIST(IORD(II)) .gt. DIST(IORD(II+IGAP))
     *  .AND. DIST(IORD(II)) .gt. ERRSIZ/2) then
        ITEMP = IORD(II)
        IORD(II) = IORD(II+IGAP)
        IORD(II+IGAP) = ITEMP
        INORD = .false.
      endif
      enddo
      if (.not. INORD) go to 2000
      go to 1000
      return
      end

C**********************************************************************C
      subroutine BoundaryDistances( NV, V, distances2Boundaries,
     2             ERRSIZ, DEBUG )
      implicit none
      real*8 V(6,NV)
      real*8 basicBoundaries(36,21)
      real*8 vtemp(6)
      integer FillPrjList
      real*8 boundary67
      real*8 boundary9A
      real*8 boundaryCD
      real*8 boundary678X
      real*8 boundary9ABCDEX
      real*8 boundaryFX
      real*8 a
      real*8 ERRSIZ
      integer nv, i, ip, iv, n
      real*8 anorm
      logical DEBUG
      character*1 XCONE
      real*8 distances2Boundaries(21)
C----------------------------------------------------------------------C
C get the projectors
      n = FillPrjList( basicBoundaries )
C convert the projectors to perp
      call MakePerpList( 21, basicBoundaries )

      do 2000 ip=1,21
C         write(*,*) "ip ",ip
C         write(*,"(6f7.3)" ) (basicBoundaries(i,ip),i=1,36)
 2000 continue


C now loop thru all the vectors
      do iv=1,nv
         if ( nv .eq. 1 ) then
             write(*,"(5x,'G6 VECTOR ' 5(f12.4,',' ),f12.4)" )
     2              (v(i,iv),i=1,6)
         else
             write(*,*) "VECTOR", iv
             write(*,"(5x,1Hg,i1, f12.7)" ) (i, v(i,iv),i=1,6)
         endif

         call rmv6( v(1,iv), basicBoundaries(1,16), vtemp )
         boundary678X = anorm(6,vtemp)
         distances2Boundaries(16) = boundary678X;
         call rmv6( v(1,iv), basicBoundaries(1,17), vtemp )
         boundary9ABCDEX = anorm(6,vtemp)
         distances2Boundaries(17) = boundary9ABCDEX;
         call rmv6( v(1,iv), basicBoundaries(1,18), vtemp )
         boundaryFX = anorm(6,vtemp)
         distances2Boundaries(18) = boundaryFX;
         call rmv6( v(1,iv), basicBoundaries(1,19), vtemp )
         boundary67 = anorm(6,vtemp)
         distances2Boundaries(19) = boundary67;
         call rmv6( v(1,iv), basicBoundaries(1,20), vtemp )
         boundary9A = anorm(6,vtemp)
         distances2Boundaries(20) = boundary9A;
         call rmv6( v(1,iv), basicBoundaries(1,21), vtemp )
         boundaryCD = anorm(6,vtemp)
         distances2Boundaries(21) = boundaryCD;


         do ip=1,15
            call rmv6( v(1,iv), basicBoundaries(1,ip), vtemp )
            XCONE = ' '
            if ((v(4,iv)*v(5,iv)*v(6,iv) .le. 0.D0
     *        .and. (ip .eq. 6 .or. ip .eq. 7
     *          .or. ip .eq. 9 .or. ip .eq. 10
     *          .or. ip .eq. 12.or. ip .eq. 13))
     *        .or. (v(4,iv)*v(5,iv)*v(6,iv) .gt. 0.D0
     *        .and. (ip .eq. 8 .or. ip .eq. 11
     *          .or. ip .eq. 14. or. ip .eq. 15))) then
            XCONE = 'X'
            if ( ip .eq. 6 .or. ip .eq. 7 .or. ip .eq. 8 )
     *        call rmv6( v(1,iv), basicBoundaries(1,16), vtemp )
            if ( ip .ge. 9 .and. ip .le. 14)
     *        call rmv6( v(1,iv), basicBoundaries(1,17), vtemp )
            if ( ip .eq. 15)
     *        call rmv6( v(1,iv), basicBoundaries(1,18), vtemp )
            endif
            a = anorm(6,vtemp)

            if ( ip .eq.  6 .and. v(5,iv) .le. v(6,iv)
     *         .and. v(5,iv) .gt. 0.D0) then
               if (boundary67 .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,Z2,4x,g10.5)") "BOUNDARY67",ip,boundary67
               endif
               distances2Boundaries(ip) = boundary67
            elseif ( ip .eq.  7 .and. v(5,iv) .ge. v(6,iv)
     *         .and. v(6,iv) .gt. 0.D0) then
               if (boundary67 .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,Z2,4x,g10.5)") "BOUNDARY67",ip,boundary67
               endif
               distances2Boundaries(ip) = boundary67
            elseif ( ip .eq.  9 .and. v(4,iv) .le. v(6,iv)
     *         .and. v(4,iv) .gt. 0.D0 ) then
               if (boundary9A .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,Z2,4x,g10.5)") "BOUNDARY9A",ip,boundary9A
               endif
               distances2Boundaries(ip) = boundary9A
            elseif ( ip .eq. 10 .and. v(4,iv) .ge. v(6,iv)
     *          .and. v(6,iv) .gt. 0.D0 ) then
               if (boundary9A .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,Z2,4x,g10.5)") "BOUNDARY9A",ip,boundary9A
               endif
               distances2Boundaries(ip) = boundary9A
            elseif ( ip .eq. 12 .and. v(4,iv) .le. v(5,iv)
     *         .and. v(4,iv) .gt. 0.D0 ) then
               if (boundaryCD .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,Z2,4x,g10.5)") "BOUNDARYCD",ip,boundaryCD
               endif
               distances2Boundaries(ip) = boundaryCD
            elseif ( ip .eq. 13 .and. v(4,iv) .ge. v(5,iv)
     *         .and. v(5,iv) .gt. 0.D0 ) then
               if (boundaryCD .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,Z2,4x,g10.5)") "BOUNDARYCD",ip,boundaryCD
               endif
               distances2Boundaries(ip) = boundaryCD
            else
               if (a .lt. ERRSIZ*3.5D0 .or. DEBUG) then
               write(*,"(a,a,1X,Z2,4x,g10.5)")
     *         "BOUNDARY", XCONE, ip, a
               endif
               distances2Boundaries(ip) = a
            endif
         enddo

         write(*,*)

      enddo
      end

C**********************************************************************C
      subroutine MakePerpList( NV, basicBoundaries )
      implicit none
      integer iprj, i, nv
      real*8 basicBoundaries(36,NV)
C----------------------------------------------------------------------C
      do 3000 iprj=1,NV
         do 1000 i=1,36
 1000    basicBoundaries(i,iprj) = -basicBoundaries(i,iprj)
         do 2000 i=1,36,7
 2000    basicBoundaries(i,iprj) = 1.0D0 + basicBoundaries(i,iprj)
 3000 continue
      end

C**********************************************************************C
      SUBROUTINE CHKVEC(V)
      implicit none
      integer i, nbad, iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      real*8 V(6)
C----------------------------------------------------------------------C
      NBAD = 0
      DO 1000 I=1,3
         IF(V(I).LE. 0.0) THEN
            WRITE (*,*) ' BAD VECTOR, I=',I,' ',V(I)
            NBAD = NBAD + 1
         ENDIF
 1000 CONTINUE
      IF (NBAD .GT. 0) STOP
      END

C**********************************************************************C
      REAL*8 FUNCTION DOTVN (N,V1,V2)
      implicit none
      real*8 V1(N),V2(N)
      integer i,n
C----------------------------------------------------------------------C
      DOTVN = 0.0
      DO 1000 I=1,N
         DOTVN = DOTVN + V1(I)*V2(I)
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE IMV6 (V1,M,V2)
      implicit none
      real*8 V1(6),V2(6)
      INTEGER M(36)
      integer i,j
      real*8 sum
C----------------------------------------------------------------------C
      DO 3000 I=1,6
      SUM = 0.0
      DO 2000 J=1,6
         SUM = SUM + M(6*(I-1)+J)*V1(J)
 2000 CONTINUE
      V2(I) = SUM
 3000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE RMV6 (V1,M,V2)
      implicit none
      real*8 V1(6),V2(6)
      real*8 M(36)
      integer i, j
      real*8 sum
C----------------------------------------------------------------------C
      DO 3000 I=1,6
      SUM = 0.0
      DO 2000 J=1,6
         SUM = SUM + M(6*(I-1)+J)*V1(J)
 2000 CONTINUE
      V2(I) = SUM
 3000 CONTINUE
      END

C**********************************************************************C
      real*8 function anorm( n, v )
      implicit none
      integer n, i
      real*8 v(n)
      real*8 sum
C----------------------------------------------------------------------C
      sum = 0
      do 1000 i=1,n
         sum = sum + v(i)*v(i)
 1000 continue
      anorm = DSQRT(sum)
      end

C**********************************************************************C
      SUBROUTINE WRCELL (C,CE,TEST)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      real*8 C(6),CE(6)
      CHARACTER *6 TEST
      integer i
      real*8 sum
C----------------------------------------------------------------------C
      IF (TEST .NE. 'WRCELL') THEN
         WRITE (*,*) ' TEST WAS WRONG IN WRCELL'
         STOP
      ENDIF

      SUM = 0.0
      DO 1000 I=1,6
 1000 SUM = SUM + CE(I)

      IF (SUM .EQ. 0) THEN
         WRITE (*,'(1X,6F10.3)') C
      ELSE
         DO 3000 I=1,6
            IF (CE(I) .GT. 0.0) THEN
               WRITE (*,'(1X,F10.3,2X,F7.3)') C(I),CE(I)
            ELSE
               WRITE (*,'(1X,F10.3,2X,F7.3)') C(I)
            ENDIF
 3000    CONTINUE
      ENDIF
      WRITE (*,*)
      END

C**********************************************************************C
      LOGICAL FUNCTION OKCELL (LATSYM,CV)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      integer nlatt, i
      PARAMETER (NLATT=16)
      CHARACTER *1 SYMLST(NLATT)
      CHARACTER *1 LATSYM
      real*8 CV(6), V(6)
      DATA SYMLST /'P','A','B','C','I','F','R','V',
     2             'p','a','b','c','i','f','r','v'/
C----------------------------------------------------------------------C
      OKCELL = .TRUE.

C     input lattice type was "V", expect that it's a vector
      if ( LATSYM .eq. 'V' ) then
         call cpyvn( 6, cv, v )
         call g6toc( v, cv, "G6TOC " )
      endif

      IF (CV(4).GT.175.0D0 .OR. CV(5).GT.175.0D0
     2              .OR. CV(6).GT.175.0D0) THEN
         WRITE (*,*) ' THE LATTICE ANGLES MUST BE LESS THAN 175 DEGREES'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(4).LT.5.0D0 .OR. CV(5).LT.5.0D0 .OR. CV(6).LT.5.0D0) THEN
         WRITE (*,*) ' THE LATTICE ANGLES MUST EXCEED 5.0 DEGREES'
         OKCELL = .FALSE.
      ENDIF

      DO 1000 I=1,NLATT
         IF (LATSYM.EQ.SYMLST(I)) THEN
            GO TO 1100
         ENDIF
 1000 CONTINUE
      WRITE (*,*) ' XTAL CLASS SYMBOL ',LATSYM,' IS NOT IMPLEMENTED'
      OKCELL = .FALSE.
 1100 CONTINUE



      IF (CV(4).GE.CV(5)+CV(6)) THEN
         WRITE (*,*) ' ERROR, ALPHA EXCEEDS BETA PLUS GAMMA'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(5).GE.CV(4)+CV(6)) THEN
         WRITE (*,*) ' ERROR, BETA EXCEEDS ALPHA PLUS GAMMA'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(6).GE.CV(4)+CV(5)) THEN
         WRITE (*,*) ' ERROR, GAMMA EXCEEDS ALPHA PLUS GAMMA'
         OKCELL = .FALSE.
      ENDIF

      IF (CV(4)+CV(5)+CV(6) .GT. 355.0D0) THEN
         WRITE (*,*) ' THE SUM OF THE CELL ANGLES MUST BE LESS THAN'
         WRITE (*,*) ' 355.0 DEGREES'
         OKCELL = .FALSE.
      ENDIF

      DO 2000 I=1,3
         IF (CV(I) .LT. 1.0D0) THEN
            WRITE (*,*) ' THE AXIAL LENGTHS MUST EXCEED 1.0)'
            OKCELL = .FALSE.
         ENDIF
 2000 CONTINUE
      END

C**********************************************************************C
      LOGICAL FUNCTION INPCEL (LATSYM,TITLE,CV,CE)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      EXTERNAL OKCELL
      LOGICAL OKCELL
      CHARACTER *1 LATSYM
      CHARACTER *100 TITLE
      CHARACTER *102 INPUTLINE
      real*8 CV(6),CE(6)
C----------------------------------------------------------------------C

 1000 CONTINUE
      WRITE (*,*) ' Input Xtal Class '
      READ (*,"(A)", END=9000 ) INPUTLINE
      LATSYM = INPUTLINE(1:1)
      TITLE = INPUTLINE(3:101)

      IF (LATSYM .GE. 'a' .AND. LATSYM .LE. 'z')
     2 LATSYM = CHAR(ICHAR(LATSYM)-ICHAR('a')+ICHAR('A'))
      WRITE (*,*) ' Input Cell Parameters'
      READ (*,*, END=9000) CV
      INPCEL = OKCELL(LATSYM,CV)
      IF (.NOT. INPCEL) GO TO 1000

      WRITE (*,*) ' Input Standard Deviations of Cell Parameters'
      READ (*,*, END=9000) CE
      RETURN
 9000 CONTINUE
      INPCEL = .FALSE.
      RETURN
      END

C**********************************************************************C
      SUBROUTINE CTOG6 (CV,CVE,G,GE,SIZE,ERRSIZ,RATIO,TEST)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      CHARACTER *6 TEST
      real*8 COSI(4:6)
      real*8 C(6),CV(6),CE(6),CVE(6),G(6),GE(6)
      integer i,j,i1,i2
      real*8 size, errsiz,ratio
      real*8 rad
      parameter (RAD = DATAN2(0.0D0,-1.0D0) / 180.0D0)
C----------------------------------------------------------------------C

      IF (TEST .NE. 'CTOG6 ') THEN
         WRITE (*,*) ' TEST IS WRONG IN CTOG6'
         STOP
      ENDIF

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
         COSI(J) = DCOS(C(J))
         IF (C(I).NE.0.0D0) THEN
            G(J)=2.0*C(1)*C(2)*C(3)*COSI(J)/C(I)
         ELSE
            G(J) = 0.0D0
         ENDIF
         GE(I) = 2.D0*ABS(C(I)*CE(I))
C         WRITE (*,*) ' GE(I) ',I,' ',GE(I)
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
         GE(J) = 2.0*DSQRT(G(I1)*(COSI(J)*CE(I2))**2 +
     2                    G(I2)*(COSI(J)*CE(I1))**2 +
     3                    G(I1)*G(I2)*(SIN(C(J))*CE(J))**2)

C
C   NOTE THE UNITS IMBALANCE ABOVE
C

C      WRITE (*,*) ' GE(J) ',J,' ',GE(J)
 2000 CONTINUE

      SIZE = 0.0D0
      ERRSIZ = 0.0D0
      DO 3000 I=1,6
         SIZE = SIZE + G(I)*G(I)
         ERRSIZ = ERRSIZ + GE(I)*GE(I)
 3000 CONTINUE

      SIZE = DSQRT(SIZE)
      ERRSIZ = DSQRT(ERRSIZ)
      RATIO = ERRSIZ / SIZE
      END

C**********************************************************************C
      SUBROUTINE WRVEC6(V,VE,TEST)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      real*8 V(6),VE(6)
      CHARACTER *6 TEST
      integer i
      real*8 sum
C----------------------------------------------------------------------C

      IF (TEST .NE. 'WRVEC6') THEN
         WRITE (*,*) ' TEST WAS WRONG IN WRVEC6'
         STOP
      ENDIF

      SUM = 0.0
      DO 1000 I=1,6
 1000 SUM = SUM + ABS(VE(I))

      IF (SUM .EQ. 0D0) THEN
         WRITE (*,'(1X,6F10.2)') V
      ELSE
         DO 3000 I=1,6
            WRITE (*,'(1X,F10.3,3X,F6.3)') V(I),VE(I)
 3000    CONTINUE
      ENDIF
      WRITE (*,*)
      END

C**********************************************************************C
      SUBROUTINE WRSIZE (SIZE,ERRSIZ,RATIO,TEST)
      implicit none
      real*8 size, errsiz, ratio
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      CHARACTER *6 TEST
C----------------------------------------------------------------------C
      IF (TEST .NE. 'WRSIZE') THEN
         WRITE (*,*) ' TEST WAS WRONG IN WRSIZE'
         STOP
      ENDIF
      WRITE (*,'('' INPUT VECTOR SIZE AND ERROR AND RATIO '',3F10.2)')
     2  SIZE,ERRSIZ,RATIO
      END

C**********************************************************************C
      SUBROUTINE ICPYVN(N,M1,M2)
      implicit none
      integer n
      INTEGER M1(N),M2(N)
      integer i
C----------------------------------------------------------------------C
      DO 1000 I=1,N
 1000 M2(I) = M1(I)
      END

C**********************************************************************C
      SUBROUTINE IUNTMN (N,M)
      implicit none
      INTEGER N, i, j
      INTEGER M(N,N)
C----------------------------------------------------------------------C
      DO 1000 I=1,N
      DO 1000 J=1,N
 1000 M(I,J) = 0
      DO 2000 I=1,N
 2000 M(I,I) = 1
      END

C**********************************************************************C
      INTEGER FUNCTION INRMSQ (N,M)
      implicit none
      INTEGER N, i
      INTEGER M(N)
C----------------------------------------------------------------------C
      INRMSQ = 0
      DO 1000 I=1,N
      INRMSQ = INRMSQ+M(I)*M(I)
 1000 CONTINUE
      RETURN
      END


C**********************************************************************C
      SUBROUTINE RUNTMN (N,M)
      implicit none
      INTEGER N,i,j
      real*8 M(N,N)
C----------------------------------------------------------------------C
      DO 1000 I=1,N
      DO 1000 J=1,N
 1000 M(I,J) = 0.0
      DO 2000 I=1,N
 2000 M(I,I) = 1.0
      END

C**********************************************************************C
      SUBROUTINE ZEROS (N,V)
      implicit none
      integer i,n
      real*8 V(N)
C----------------------------------------------------------------------C
      DO 1000 I=1,N
 1000 V(I) = 0.0
      END

C**********************************************************************C
      SUBROUTINE MKPRIM (LATSYM,GIN,M1,GOUT,TEST)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      CHARACTER *6 TEST
      CHARACTER LATSYM
      real*8 GIN(6),GOUT(6)
      real*8 M(36),M1(36)
      integer i
C----------------------------------------------------------------------C

      IF (TEST .NE. 'MKPRIM') THEN
         WRITE (*,*) ' TEST WAS WRONG IN MKPRIM'
         STOP
      ENDIF
      CALL ZEROS (36,M)

      IF (LATSYM .EQ. 'P' .OR. LATSYM .EQ. 'V') THEN
         CALL RUNTMN(6,M)

      ELSEIF (LATSYM .EQ. 'I') THEN
         M(1) = 1D0
         M(8) = 1D0
         DO 1000 I=13,18
 1000    M(I) = 0.25D0
         M(20) = 1D0
         M(22) = 0.5D0
         M(24) = 0.5D0
         M(25) = 1D0
         M(29) = 0.5D0
         M(30) = 0.5D0
         M(36) = 1D0

      ELSEIF (LATSYM .EQ. 'F') THEN
         M(1) = 0.25D0
         M(2) = 0.25D0
         M(6) = 0.25D0
         M(7) = 0.25D0
         M(9) = 0.25D0
         M(11) = 0.25D0
         M(14) = 0.25D0
         M(15) = 0.25D0
         M(16) = 0.25D0
         M(21) = 0.5D0
         DO 2000 I=22,24
 2000    M(I)= 0.25D0
         M(26) = 0.5D0
         DO 2100 I=28,30
 2100    M(I) = 0.25D0
         M(31) = 0.5D0
         DO 2200 I=34,36
 2200    M(I) = 0.25D0
      ELSEIF (LATSYM .EQ. 'A') THEN
         M(1) = 1D0
         M(8) = 1D0
         M(14) = 0.25D0
         M(15) = 0.25D0
         M(16) = 0.25D0
         M(20) = 1D0
         M(22) = 0.5D0
         M(29) = 0.5D0
         M(30) = 0.5D0
         M(36) = 1D0
      ELSEIF (LATSYM .EQ. 'B') THEN
         M(1) = 1D0
         M(8) = 1D0
         M(13) = 0.25D0
         M(15) = 0.25D0
         M(17) = 0.25D0
         M(22) = 0.5D0
         M(24) = 0.5D0
         M(25) = 1D0
         M(29) = 0.5D0
         M(36) = 1D0
      ELSEIF (LATSYM .EQ. 'C') THEN
         M(1) = 1D0
         M(7) = 0.25D0
         M(8) = 0.25D0
         M(12) = 0.25D0
         M(15) = 1D0
         M(22) = 0.5D0
         M(23) = 0.5D0
         M(29) = 1D0
         M(31) = 1D0
         M(36) = 0.5D0

      ELSE
         WRITE (*,*) ' DIDN-T FIND LATTICE SYMBOL ',LATSYM
         STOP
      ENDIF
      CALL CPYVN(36,M,M1)
      CALL RMV6(GIN,M1,GOUT)
      END

C**********************************************************************C
      SUBROUTINE IZEROS (N,A)
      implicit none
      integer i,n
      INTEGER A(N)
C----------------------------------------------------------------------C
      DO 1000 I=1,N
 1000 A(I) = 0
      END

C**********************************************************************C
      SUBROUTINE MKNORM (VI,M,VOUT,TEST)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      integer i, ncycle, numneg, nzero
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      CHARACTER *6 TEST
      real*8 VIN(6),VI(6),VOUT(6)
      logical PPPMMM
      INTEGER M(36),MAT(36),M1(36)
      LOGICAL AGAIN
C----------------------------------------------------------------------C

      IF (TEST .NE. 'MKNORM') THEN
         WRITE (*,*) ' TEST WAS WRONG IN MKNORM'
         STOP
      ENDIF
C      CALL CPYVN(6,VI,VIN)
      DO 500 I=1,6
  500 VIN(I) = VI(I)
      CALL IUNTMN(6,M)
      NCYCLE = 0
 1000 CONTINUE
      NCYCLE = NCYCLE + 1
      AGAIN =.FALSE.
      CALL IZEROS (36,MAT)
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
         CALL IMIM6(MAT,M,M1)
         CALL ICPYVN(36,M1,M)
         CALL IMV6(VIN,MAT,VOUT)
         CALL CPYVN(6,VOUT,VIN)
      ENDIF
      IF (AGAIN .AND. NCYCLE.LT.4) GO TO 1000

      NUMNEG = 0
      NZERO = 0
      PPPMMM = .true.
      DO 2000 I=4,6
         IF (VIN(I).LT.1.0D-10) THEN
           PPPMMM= .NOT. PPPMMM
           NUMNEG = NUMNEG + 1
         ENDIF
         IF (DABS(VIN(I)).LT.1.0D-10) NZERO = NZERO+1
 2000 CONTINUE
      IF (NZERO .GT. 0) PPPMMM = .false.
      CALL IUNTMN(6,MAT)
      IF (PPPMMM .AND. NUMNEG.NE.0) THEN
        DO I = 4,6
          IF (VIN(I).LT.1.0D-10) THEN
            MAT(6*(I-1)+I) = -1
          ENDIF
        ENDDO
      ENDIF
      IF ((.NOT.PPPMMM) .AND. NUMNEG.NE.3) THEN
        DO I = 4,6
          IF (VIN(I).GE.1.0D-10) THEN
            MAT(6*(I-1)+I) = -1
          ENDIF
        ENDDO
      ENDIF
      CALL IMIM6(MAT,M,M1)
      CALL ICPYVN(36,M1,M)
      CALL IMV6(VIN,MAT,VOUT)
      CALL CPYVN(6,VOUT,VIN)
      END

C**********************************************************************C
      SUBROUTINE REDUCE (VI,M,VOUT,DELTA,TEST)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      CHARACTER *6 TEST

      real*8 VI(6),VIN(6), VOUT(6)
      real*8 DELTA
      LOGICAL AGAIN
      INTEGER M(36),M1(36),M2(36)
      real*8 ZEROS(6)
      integer i, last, ncycle
      DATA ZEROS /6*0.0/
C----------------------------------------------------------------------C

      IF (TEST .NE. 'REDUCE') THEN
         WRITE (*,*) ' TEST WAS WRONG IN REDUCE'
         STOP
      ENDIF
      CALL CPYVN(6,VI,VIN)
      NCYCLE = 0
 1000 CONTINUE
      LAST = 0
      CALL MKNORM (VIN,M,VOUT,'MKNORM')
      CALL CPYVN(6,VOUT,VIN)
      CALL IUNTMN (6,M1)
      IF (ABS(VIN(4)) .GT. ABS(VIN(2)) + DELTA) THEN
          M1(14) = 1
          M1(16) = INT(-dsign (1.0D0,VIN(4)))
          M1(20) = INT(-2.0*dsign(1.0D0,VIN(4)))
          M1(30) = M1(16)
          AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
          LAST = 5
      ELSEIF (ABS(VIN(5)) .GT. ABS(VIN(1)) + DELTA) THEN
         M1(13) = 1
         M1(17) = INT(-dsign (1.0D0,VIN(5)))
         M1(24) = M1(17)
         M1(25) = 2*M1(17)
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 6
      ELSEIF (ABS(VIN(6)) .GT. ABS(VIN(1))+ DELTA) THEN
         M1(7) = 1
         M1(12) = INT(-dsign(1.0D0,VIN(6)))
         M1(23) = M1(12)
         M1(31) = 2*M1(12)
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 7
      ELSEIF (VIN(4)+VIN(5)+VIN(6)+ABS(VIN(1))+ABS(VIN(2))+DELTA
     2  .LT. 0.0D0) THEN
         DO 2000 I=13,18
 2000    M1(I) = 1
         M1(20) = 2
         M1(24) = 1
         M1(25) = 2
        M1(30) = 1
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 8
      ELSEIF ( (ABS(VIN(4)-VIN(2)).LE.DELTA
     2            .AND. 2.0D0*VIN(5)+DELTA.LT.VIN(6)) .OR.
     3         (ABS(VIN(4)+VIN(2)).LE.DELTA
     4            .AND. VIN(6)+DELTA.LT.0.0D0)) THEN
         M1(14) = 1
         M1(16) = INT(-dsign(1.0D0,VIN(4)))
         M1(20) = 2*M1(16)
         M1(30) = M1(16)
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 15
      ELSEIF ( (ABS(VIN(5)-VIN(1)).LE.DELTA
     2            .AND. 2.0D0*VIN(4)+DELTA.LT.VIN(6)) .OR.
     3         (ABS(VIN(5)+VIN(1)).LE.DELTA
     4            .AND. VIN(6)+DELTA.LT.0.0D0)) THEN
         M1(13) = 1
         M1(17) = INT(-dsign(1.0D0,VIN(5)))
         M1(24) = M1(17)
         M1(25) = 2*M1(17)
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 16
      ELSEIF ( (ABS(VIN(6)-VIN(1)).LE.DELTA
     2            .AND. 2.0D0*VIN(4)+DELTA.LT.VIN(5)) .OR.
     3         (ABS(VIN(6)+VIN(1)).LE.DELTA
     4            .AND. VIN(5)+DELTA.LT.0.0D0)) THEN
         M1(7) = 1
         M1(12) = INT(-dsign(1.0D0,VIN(6)))
         M1(23) = M1(12)
         M1(31) = 2*M1(12)
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 17
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 123
      ELSEIF ( ABS((VIN(4)+VIN(5)+VIN(6)
     2 +ABS(VIN(1))+ABS(VIN(2))))
     3 .LE.DELTA .AND.
     4 ( 2.0D0*(ABS(VIN(1))+VIN(5))+VIN(6).GT.DELTA) ) THEN
         DO 3000 I=13,18
 3000    M1(I) = 1
         M1(20) = 2
         M1(24) = 1
         M1(25) = 2
         M1(30) = 1
         AGAIN = .TRUE.
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN (36,M2,M)
         CALL IMV6(VIN,M1,VOUT)
         LAST = 18
      ELSE
         AGAIN = .FALSE.
         CALL CPYVN (6,VIN,VOUT)
      ENDIF

      CALL MKNORM (VOUT,M1,VIN,'MKNORM')
         CALL IMIM6(M1,M,M2)
         CALL ICPYVN(36,M2,M)
         CALL CPYVN(6,VIN,VOUT)

      IF (VIN(1).LT. 0.0D0 .OR. VIN(2).LT.0.0D0
     2        .OR. VIN(3).LT.0.0D0) THEN
         WRITE (*,*) ' NEG. SQ. AXIS ',NCYCLE
         CALL WRVEC6(VIN,ZEROS,'WRVEC6')
         CALL WRVEC6(VOUT,ZEROS,'WRVEC6')
         READ (*,*)
      ENDIF

      NCYCLE = NCYCLE + 1
      IF (NCYCLE .LT. 25 .AND. AGAIN) GO TO 1000
      END


C**********************************************************************C
      SUBROUTINE G6TOC (G,C,TEST)
      implicit none
      integer  iunit0,iunit1,iunit2,iunit3,iunit10
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      CHARACTER *6 TEST
      real*8 G(6),C(6)
      real*8 rad
      parameter (rad=180.0D0/4.0D0/datan(1.0D0))
      real*8 AC
      integer i
C----------------------------------------------------------------------C
      IF (TEST .NE. 'G6TOC ') THEN
         WRITE (*,*) ' TEST WAS WRONG IN G6TOC'
         STOP
      ENDIF
      DO 900 I=1,3
         IF (G(I) .LE. 0.0) THEN
            WRITE (*,*) ' G(I)<=0, I=',I,'  ',G(I)
         ENDIF
  900 CONTINUE

      DO 1000 I=1,3
 1000 C(I) = DSQRT(G(I))
      AC = 0.5D0*G(4)/C(2)/C(3)
      IF (ABS(AC) .LE. 1.0) THEN
         C(4) = rad*DACOS(AC)
      ELSE
         C(4) = 0.0
         WRITE (*,*) ' ARG>1.0 C(4) ',AC,1.0D0-ABS(AC)
      ENDIF

      AC = 0.5D0*G(5)/C(1)/C(3)
      IF (ABS(AC) .LE. 1.0D0) THEN
         C(5) = rad*DACOS(AC)
      ELSE
         C(5) = 0.0
         WRITE (*,*) ' ARG>1.0 C(5) ',AC,1.0D0-ABS(AC)
      ENDIF

      AC = 0.5D0*G(6)/C(1)/C(2)
      IF (ABS(AC) .LE. 1.0D0) THEN
         C(6) = rad*DACOS(AC)
      ELSE
         C(6) = 0.0D0
         WRITE (*,*) ' ARG>1.0 C(6) ',AC,1.0D0-ABS(AC)
      ENDIF
      END

C**********************************************************************C
      SUBROUTINE IMIM6(M1,M2,M3)
      implicit none
      INTEGER M1(6,6),M2(6,6),M3(6,6)
      integer i,j,k
C----------------------------------------------------------------------C
      CALL IZEROS (36,M3)
      DO 3000 I=1,6
      DO 2000 J=1,6
        M3(I,J) = 0
        DO 1000 K=1,6
           M3(I,J) = M3(I,J) + M1(I,K)*M2(K,J)
 1000   CONTINUE
 2000 CONTINUE
 3000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE CPYVN (N,V1,V2)
      implicit none
      integer i,n
      real*8 V1(N),V2(N)
C----------------------------------------------------------------------C
      DO 1000 I=1,N
 1000 V2(I) = V1(I)
      END

C**********************************************************************C
      subroutine bldprj (maxprj,nproj,itdesg,chrlat,pjn,prj,test)
      implicit none
      integer iunit0,iunit1,iunit2,iunit3,iunit10
      integer maxprj, i, ii, iproj, itemp, j, nprob
      common /files/ iunit0,iunit1,iunit2,iunit3,iunit10
      real*8 pjn(maxprj)
      integer prj(36,MAXPRJ)
      integer itdesg(maxprj)
      character *2 chrlat(maxprj)
      integer ngtype(42)
      integer projct(36,42)
      integer iprj(6,6,42)
      integer nproj
      equivalence (projct,iprj)
      real*8 pjnorm(42)
      real*8 sum
      character *2 lattyp(42)
      character *6 test

      data ngtype(1) /-3/
      data lattyp(1) /'cP'/
      data pjnorm(1) /3D0/
      data( projct(i,1),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(2) /-5/
      data lattyp(2) /'cI'/
      data pjnorm(2) /39D0/
      data (projct(i,2),i=1,36) /
     1   9, 9, 9,-6,-6,-6,
     2   9, 9, 9,-6,-6,-6,
     3   9, 9, 9,-6,-6,-6,
     4  -6,-6,-6, 4, 4, 4,
     5  -6,-6,-6, 4, 4, 4,
     6  -6,-6,-6, 4, 4, 4 /


      data ngtype(3) /1/
      data lattyp(3) /'cF'/
      data pjnorm(3) /6D0/
      data (projct(i,3),i=1,36) /
     1  1,1,1,1,1,1,
     2  1,1,1,1,1,1,
     3  1,1,1,1,1,1,
     4  1,1,1,1,1,1,
     5  1,1,1,1,1,1,
     6  1,1,1,1,1,1 /



      data ngtype(4) /-11/
      data lattyp(4) /'tP'/
      data pjnorm(4) /2D0/
      data (projct(i,4),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(5) /-21/
      data lattyp(5) /'tP'/
      data pjnorm(5) /2D0/
      data (projct(i,5),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(6) /-15/
      data lattyp(6) /'tI'/
      data pjnorm(6) /4D0/
      data (projct(i,6),i=1,36) /
     1  1, 1, 0,-1,-1, 0,
     2  1, 1, 0,-1,-1, 0,
     3  0, 0, 4, 0, 0, 0,
     4 -1,-1, 0, 1, 1, 0,
     5 -1,-1, 0, 1, 1, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(7) /-6/
      data lattyp(7) /'tI'/
      data pjnorm(7) /26D0/
      data (projct(i,7),i=1,36) /
     1  6, 6, 6,-4,-4,-4,
     2  6, 6, 6,-4,-4,-4,
     3  6, 6, 6,-4,-4,-4,
     4 -4,-4,-4, 7, 7,-6,
     5 -4,-4,-4, 7, 7,-6,
     6 -4,-4,-4,-6,-6,20 /


      data ngtype(8) /-7/
      data lattyp(8) /'tI'/
      data pjnorm(8) /26D0/
      data (projct(i,8),i=1,36) /
     1  6, 6, 6,-4,-4,-4,
     2  6, 6, 6,-4,-4,-4,
     3  6, 6, 6,-4,-4,-4,
     4 -4,-4,-4,20,-6,-6,
     5 -4,-4,-4,-6, 7, 7,
     6 -4,-4,-4,-6, 7, 7 /


      data ngtype(9) /18/
      data lattyp(9) /'tI'/
      data pjnorm(9) /26D0/
      data (projct(i,9),i=1,36) /
     1  8, 0, 0, 4, 8, 8,
     2  0,13,13, 0, 0, 0,
     3  0,13,13, 0, 0, 0,
     4  4, 0, 0, 2, 4, 4,
     5  8, 0, 0, 4, 8, 8,
     6  8, 0, 0, 4, 8, 8 /


      data ngtype(10) /-12/
      data lattyp(10) /'hP'/
      data pjnorm(10) /3D0/
      data (projct(i,10),i=1,36) /
     1  1, 1, 0, 0, 0,-1,
     2  1, 1, 0, 0, 0,-1,
     3  0, 0, 3, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1,-1, 0, 0, 0, 1 /


      data ngtype(11) /-22/
      data lattyp(11) /'hP'/
      data pjnorm(11) /3D0/
      data (projct(i,11),i=1,36) /
     1  3, 0, 0, 0, 0, 0,
     2  0, 1, 1,-1, 0, 0,
     3  0, 1, 1,-1, 0, 0,
     4  0,-1,-1, 1, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(12) /9/
      data lattyp(12) /'hR'/
      data pjnorm(12) /5D0/
      data (projct(i,12),i=1,36) /
     1  1,1,0,1,1,1,
     2  1,1,0,1,1,1,
     3  0,0,5,0,0,0,
     4  1,1,0,1,1,1,
     5  1,1,0,1,1,1,
     6  1,1,0,1,1,1  /


      data ngtype(13) /2/
      data lattyp(13) /'hR'/
      data pjnorm(13) /3D0/
      data (projct(i,13),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,1,1,1,
     5  0,0,0,1,1,1,
     6  0,0,0,1,1,1 /


      data ngtype(14) /-4/
      data lattyp(14) /'hR'/
      data pjnorm(14) /3D0/
      data (projct(i,14),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,1,1,1,
     5  0,0,0,1,1,1,
     6  0,0,0,1,1,1 /


      data ngtype(15) /-24/
      data lattyp(15) /'hR'/
      data pjnorm(15) /53D0/
      data (projct(i,15),i=1,36) /
     1  27,  3,  3,  6,-18,-18,
     2   3, 18, 18,-17, -2, -2,
     3   3, 18, 18,-17, -2, -2,
     4   6,-17,-17, 19, -4, -4,
     5 -18, -2, -2, -4, 12, 12,
     6 -18, -2, -2, -4, 12, 12 /


      data ngtype(16) /-32/
      data lattyp(16) /'oP'/
      data pjnorm(16) /1D0/
      data (projct(i,16),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(17) /-36/
      data lattyp(17) /'oS'/
      data pjnorm(17) /2D0/
      data (projct(i,17),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(18) /-38/
      data lattyp(18) /'oS'/
      data pjnorm(18) /2D0/
      data (projct(i,18),i=1,36) /
     1  1, 0, 0, 0, 0,-1,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1, 0, 0, 0, 0, 1 /


      data ngtype(19) /-13/
      data lattyp(19) /'oS'/
      data pjnorm(19) /2D0/
      data (projct(i,19),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,2 /


      data ngtype(20) /-23/
      data lattyp(20) /'oS'/
      data pjnorm(20) /2D0/
      data (projct(i,20),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(21) /-40/
      data lattyp(21) /'oS'/
      data pjnorm(21) /2D0/
      data (projct(i,21),i=1,36) /
     1  2, 0, 0, 0, 0, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(22) /-16/
      data lattyp(22) /'oF'/
      data pjnorm(22) /10D0/
      data (projct(i,22),i=1,36) /
     1  3, 3, 0,-2,-2,-2,
     2  3, 3, 0,-2,-2,-2,
     3  0, 0,10, 0, 0, 0,
     4 -2,-2, 0, 3, 3,-2,
     5 -2,-2, 0, 3, 3,-2,
     6 -2,-2, 0,-2,-2, 8  /


      data ngtype(23) /26/
      data lattyp(23) /'oF'/
      data pjnorm(23) /13D0/
      data (projct(i,23),i=1,36) /
     1  4, 0, 0, 2, 4, 4,
     2  0,13, 0, 0, 0, 0,
     3  0, 0,13, 0, 0, 0,
     4  2, 0, 0, 1, 2, 2,
     5  4, 0, 0, 2, 4, 4,
     6  4, 0, 0, 2, 4, 4  /


      data ngtype(24) /-8/
      data lattyp(24) /'oF'/
      data pjnorm(24) /13D0/
      data (projct(i,24),i=1,36) /
     1  3, 3, 3,-2,-2,-2,
     2  3, 3, 3,-2,-2,-2,
     3  3, 3, 3,-2,-2,-2,
     4 -2,-2,-2,10,-3,-3,
     5 -2,-2,-2,-3,10,-3,
     6 -2,-2,-2,-3,-3,10  /


      data ngtype(25) /19/
      data lattyp(25) /'oI'/
      data pjnorm(25) /6D0/
      data (projct(i,25),i=1,36) /
     1  2,0,0,0,2,2,
     2  0,3,3,0,0,0,
     3  0,3,3,0,0,0,
     4  0,0,0,6,0,0,
     5  2,0,0,0,2,2,
     6  2,0,0,0,2,2  /


      data ngtype(26) /-42/
      data lattyp(26) /'oI'/
      data pjnorm(26) /2D0/
      data (projct(i,26),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(27) /-33/
      data lattyp(27) /'mP'/
      data pjnorm(27) /1D0/
      data (projct(i,27),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,1,0,
     6  0,0,0,0,0,0  /


      data ngtype(28) /-35/
      data lattyp(28) /'mP'/
      data pjnorm(28) /1D0/
      data (projct(i,28),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,1,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0  /


      data ngtype(29) /-34/
      data lattyp(29) /'mP'/
      data pjnorm(29) /1D0/
      data (projct(i,29),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,1  /


      data ngtype(30) /-39/
      data lattyp(30) /'mS'/
      data pjnorm(30) /2D0/
      data (projct(i,30),i=1,36) /
     1  1, 0, 0, 0, 0,-1,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 2, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1, 0, 0, 0, 0, 1  /


      data ngtype(31) /-41/
      data lattyp(31) /'mS'/
      data pjnorm(31) /2D0/
      data (projct(i,31),i=1,36) /
     1  2, 0, 0, 0, 0, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5  0, 0, 0, 0, 2, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(32) /-37/
      data lattyp(32) /'mS'/
      data pjnorm(32) /2D0/
      data (projct(i,32),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 2, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(33) /10/
      data lattyp(33) /'mS'/
      data pjnorm(33) /2D0/
      data (projct(i,33),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,1,1,0,
     5  0,0,0,1,1,0,
     6  0,0,0,0,0,2  /


      data ngtype(34) /-14/
      data lattyp(34) /'mS'/
      data pjnorm(34) /2D0/
      data (projct(i,34),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,1,1,0,
     5  0,0,0,1,1,0,
     6  0,0,0,0,0,2  /


      data ngtype(35) /20/
      data lattyp(35) /'mS'/
      data pjnorm(35) /2D0/
      data (projct(i,35),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,1,1,
     6  0,0,0,0,1,1  /


      data ngtype(36) /-25/
      data lattyp(36) /'mS'/
      data pjnorm(36) /2D0/
      data (projct(i,36),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,1,1,
     6  0,0,0,0,1,1  /


      data ngtype(37) /28/
      data lattyp(37) /'mS'/
      data pjnorm(37) /10D0/
      data (projct(i,37),i=1,36) /
     1  5, 0, 0, 0, 5, 0,
     2  0,10, 0, 0, 0, 0,
     3  0, 0,10, 0, 0, 0,
     4  0, 0, 0, 2, 0, 4,
     5  5, 0, 0, 0, 5, 0,
     6  0, 0, 0, 4, 0, 8  /


      data ngtype(38) /30/
      data lattyp(38) /'mS'/
      data pjnorm(38) /10D0/
      data (projct(i,38),i=1,36) /
     1  10, 0, 0, 0, 0, 0,
     2   0, 5, 0, 5, 0, 0,
     3   0, 0,10, 0, 0, 0,
     4   0, 5, 0, 5, 0, 0,
     5   0, 0, 0, 0, 2, 4,
     6   0, 0, 0, 0, 4, 8  /


      data ngtype(39) /29/
      data lattyp(39) /'mS'/
      data pjnorm(39) /10D0/
      data (projct(i,39),i=1,36) /
     1  5, 0, 0, 0, 0, 5,
     2  0,10, 0, 0, 0, 0,
     3  0, 0,10, 0, 0, 0,
     4  0, 0, 0, 2, 4, 0,
     5  0, 0, 0, 4, 8, 0,
     6  5, 0, 0, 0, 0, 5  /


      data ngtype(40) /-43/
      data lattyp(40) /'mS'/
      data pjnorm(40) /20D0/
      data (projct(i,40),i=1,36) /
     1  11, 1, 0, 1,-9,-4,
     2   1,11, 0,-9, 1,-4,
     3   0, 0,20, 0, 0, 0,
     4   1,-9, 0,11, 1,-4,
     5  -9, 1, 0, 1,11,-4,
     6  -4,-4, 0,-4,-4,16  /


      data ngtype(41) /-17/
      data lattyp(41) /'mS'/
      data pjnorm(41) /10D0/
      data (projct(i,41),i=1,36) /
     1  3, 3, 0,-2,-2,-2,
     2  3, 3, 0,-2,-2,-2,
     3  0, 0,10, 0, 0, 0,
     4 -2,-2, 0, 8,-2,-2,
     5 -2,-2, 0,-2, 8,-2,
     6 -2,-2, 0,-2,-2, 8  /


      data ngtype(42) /27/
      data lattyp(42) /'mS'/
      data pjnorm(42) /3D0/
      data (projct(i,42),i=1,36) /
     1  1,0,0,0,1,1,
     2  0,3,0,0,0,0,
     3  0,0,3,0,0,0,
     4  0,0,0,3,0,0,
     5  1,0,0,0,1,1,
     6  1,0,0,0,1,1  /
C----------------------------------------------------------------------C

      nprob = 0
      if (test .ne. 'BLDPRJ') then
         write (*,*) ' test string was not BLDPRJ in that routine'
         stop
      endif
      do 4000 iproj=1,42
C        if (ngtype(iproj).le. 0 .or. ngtype(iproj).gt.43) then
C            write (*,*) ' bad ngtype ',iproj,ngtype(iproj)
C            nprob = nprob + 1
C        endif
        if (lattyp(iproj).eq.' ') then
           write (*,*) ' blank lattyp ',iproj
           nprob = nprob + 1
        endif
        if (pjnorm(iproj) .le. 0D0 .or. pjnorm(iproj).gt. 60D0) then
           write (*,*) ' bad pjnorm ',iproj,' ',pjnorm(iproj)
           nprob = nprob + 1
        endif
        do 1000 i=1,5
        do 1000 j=i+1,6
          if (iprj(i,j,iproj) .ne. iprj(i,j,iproj)) then
             write (*,*) ' bad projector ',iproj,' ',i,' ',j
             write (*,*) iprj(i,j,iproj),' ',iprj(i,j,iproj)
             nprob = nprob + 1
          endif
 1000   continue
        do 2000 i=1,36
          if (projct(i,iproj)/pjnorm(iproj).gt. 1.00001D0) then
             write (*,*) ' bad projector, value > 1.0 '
             write (*,*) iproj,' pjnorm ',pjnorm(iproj),
     2    ' i ',i,' ',projct(i,iproj)
             write (*,*)
             nprob = nprob + 1
          endif
 2000    continue
        do 3500 i=1,3
        sum = 0.0
        do 3000 j=1,6
          sum = sum + abs(iprj(i,j,iproj))
 3000   continue
        if (sum .le. 0) then
           write (*,*) ' zero xyz row, ',iproj
           nprob = nprob + 1
        endif
 3500   continue
 4000 continue

      NPROJ = 0
      do 5100 i=1,MIN(MAXPRJ,42)
         IF (IABS(NGTYPE(I)) .GT. 0) THEN
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
         write (*,*) ' ',nprob,' problems were found with projectors'
         stop
      endif

      itemp = projct(1,1)
      projct(1,1) = itemp

      end

C**********************************************************************C
      integer function FillPrjList( eightfoldProjectors )
      implicit none
      real*8 eightfoldProjectors(36,*)
      character*8 level8
      real*8 prj(36,21)
      integer i, ip

      integer mxtree
      parameter (mxtree = 1000)
C----------------------------------------------------------------------C


C The projectors for the 15 base types (5-D boundaries
C in G6), plus a few extra for internal boundaries
C     case1
      data (prj(i,1),i=1,36) /
     * .5,.5,0,0,0,0,
     * .5,.5,0,0,0,0,
     * 0,0,1,0,0,0,
     * 0,0,0,1,0,0,
     * 0,0,0,0,1,0,
     * 0,0,0,0,0,1 /

C     case2
      data (prj(i,2),i=1,36) /
     * 1,0,0,0,0,0,
     * 0,.5,.5,0,0,0,
     * 0,.5,.5,0,0,0,
     * 0,0,0,1,0,0,
     * 0,0,0,0,1,0,
     * 0,0,0,0,0,1 /

C     case3
      data (prj(i,3),i=1,36) /
     * 1,0,0,0,0,0,
     * 0,1,0,0,0,0,
     * 0,0,1,0,0,0,
     * 0,0,0,0,0,0,
     * 0,0,0,0,1,0,
     * 0,0,0,0,0,1 /

C     case4
      data (prj(i,4),i=1,36) /
     * 1,0,0,0,0,0,
     * 0,1,0,0,0,0,
     * 0,0,1,0,0,0,
     * 0,0,0,1,0,0,
     * 0,0,0,0,0,0,
     * 0,0,0,0,0,1 /

C     case5
      data (prj(i,5),i=1,36) /
     * 1,0,0,0,0,0,
     * 0,1,0,0,0,0,
     * 0,0,1,0,0,0,
     * 0,0,0,1,0,0,
     * 0,0,0,0,1,0,
     * 0,0,0,0,0,0 /

C     case6
      data (prj(i,6),i=1,36) /
     * 1,  0, 0,  0, 0, 0,
     * 0, .5, 0, .5, 0, 0,
     * 0,  0, 1,  0, 0, 0,
     * 0, .5, 0, .5, 0, 0,
     * 0,  0, 0,  0, 1, 0,
     * 0,  0, 0,  0, 0, 1 /

C     case7
      data (prj(i,7),i=1,36) /
     * 1, 0, 0, 0, 0, 0,
     * 0,.5, 0,.5, 0, 0,
     * 0, 0, 1, 0, 0, 0,
     * 0,.5, 0,.5, 0, 0,
     * 0, 0, 0, 0, 1, 0,
     * 0, 0, 0, 0, 0, 1 /

C     case8
      data (prj(i,8),i=1,36) /
     * 1,  0, 0,  0, 0, 0,
     * 0, .5, 0,-.5, 0, 0,
     * 0,  0, 1,  0, 0, 0,
     * 0,-.5, 0, .5, 0, 0,
     * 0,  0, 0,  0, 1, 0,
     * 0,  0, 0,  0, 0, 1 /

C     case9
      data (prj(i,9),i=1,36) /
     * .5, 0, 0, 0,.5, 0,
     *  0, 1, 0, 0, 0, 0,
     *  0, 0, 1, 0, 0, 0,
     *  0, 0, 0, 1, 0, 0,
     * .5, 0, 0, 0,.5, 0,
     *  0, 0, 0, 0, 0, 1 /

C     caseA
      data (prj(i,10),i=1,36) /
     * .5, 0, 0, 0,.5, 0,
     *  0, 1, 0, 0, 0, 0,
     *  0, 0, 1, 0, 0, 0,
     *  0, 0, 0, 1, 0, 0,
     * .5, 0, 0, 0,.5, 0,
     *  0, 0, 0, 0, 0, 1 /

C     caseB
      data (prj(i,11),i=1,36) /
     *  .5, 0, 0, 0,-.5, 0,
     *   0, 1, 0, 0,  0, 0,
     *   0, 0, 1, 0,  0, 0,
     *   0, 0, 0, 1,  0, 0,
     * -.5, 0, 0, 0, .5, 0,
     *   0, 0, 0, 0,  0, 1 /

C     caseC
      data (prj(i,12),i=1,36) /
     * .5, 0, 0, 0, 0,.5,
     *  0, 1, 0, 0, 0, 0,
     *  0, 0, 1, 0, 0, 0,
     *  0, 0, 0, 1, 0, 0,
     *  0, 0, 0, 0, 1, 0,
     * .5, 0, 0, 0, 0,.5 /

C     caseD
      data (prj(i,13),i=1,36) /
     * .5, 0, 0, 0, 0,.5,
     *  0, 1, 0, 0, 0, 0,
     *  0, 0, 1, 0, 0, 0,
     *  0, 0, 0, 1, 0, 0,
     *  0, 0, 0, 0, 1, 0,
     * .5, 0, 0, 0, 0,.5 /

C     caseE
      data (prj(i,14),i=1,36) /
     *  .5, 0, 0, 0, 0,-.5,
     *   0, 1, 0, 0, 0,  0,
     *   0, 0, 1, 0, 0,  0,
     *   0, 0, 0, 1, 0,  0,
     *   0, 0, 0, 0, 1,  0,
     * -.5, 0, 0, 0, 0, .5 /

C     caseF
      data (prj(i,15),i=1,36) /
     *  .8,-.2, 0,-.2,-.2,-.2,
     * -.2, .8, 0,-.2,-.2,-.2,
     *   0,  0, 1,  0,  0,  0,
     * -.2,-.2, 0, .8,-.2,-.2,
     * -.2,-.2, 0,-.2, .8,-.2,
     * -.2,-.2, 0,-.2,-.2, .8 /

C     case678X g2=g4=g5=g6=0
      data (prj(i,16),i=1,36) /
     * 1, 0, 0, 0,  0,  0,
     * 0, 0, 0, 0,  0,  0,
     * 0, 0, 1, 0,  0,  0,
     * 0, 0, 0, 0,  0,  0,
     * 0, 0, 0, 0,  0,  0,
     * 0, 0, 0, 0,  0,  0 /

C     case9ABCDEX g1=g4=g5=g6=0
      data (prj(i,17),i=1,36) /
     * 0, 0, 0,  0,  0,  0,
     * 0, 1, 0,  0,  0,  0,
     * 0, 0, 1,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0 /

C     caseFX g1=g2=g4=g5=g6=0
      data (prj(i,18),i=1,36) /
     * 0, 0, 0,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0,
     * 0, 0, 1,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0,
     * 0, 0, 0,  0,  0,  0 /

C     DATA boundary67 /
      data (prj(i,19),i=1,36) /
     * 1,  0, 0,  0,  0,  0,
     * 0, .5, 0, .5,  0,  0,
     * 0,  0, 1,  0,  0,  0,
     * 0, .5, 0, .5,  0,  0,
     * 0,  0, 0,  0, .5, .5,
     * 0,  0, 0,  0, .5, .5 /

C     DATA boundary9A /
      data (prj(i,20),i=1,36) /
     * .5, 0, 0,  0, .5,  0,
     *  0, 1, 0,  0,  0,  0,
     *  0, 0, 1,  0,  0,  0,
     *  0, 0, 0, .5,  0, .5,
     * .5, 0, 0,  0, .5,  0,
     *  0, 0, 0, .5,  0, .5 /

C     DATA boundaryCD /
      data (prj(i,21),i=1,36) /
     * .5, 0, 0,  0,  0, .5,
     *  0, 1, 0,  0,  0,  0,
     *  0, 0, 1,  0,  0,  0,
     *  0, 0, 0, .5, .5,  0,
     *  0, 0, 0, .5, .5,  0,
     * .5, 0, 0,  0,  0, .5 /

      data level8 /"12679ACD"/


      do 1000 ip=1,21
      do 1000 i=1,36
 1000 eightfoldProjectors(i,ip) = prj(i,ip)

      FillPrjList = 21
      end

C**********************************************************************C
      subroutine FillList( nlevels, list, vp, stack )
      implicit none
      integer nlevels, list(nlevels), vp(nlevels), stack(nlevels)
      integer i
C----------------------------------------------------------------------C
      do 1000 i=1,nlevels
 1000 stack(i) = list(vp(i))
      end

C**********************************************************************C
      logical function NextIndex( nlevels, maxprj, stack )
      implicit none
C after a projector has been completed, look for the next (lexically)
C one to treat
      integer stack(8)
      integer nstack
      integer i, k, nlevels, maxprj
C----------------------------------------------------------------------C
      nstack = nlevels
c      write(*,*) "in NEXTPRJ"
C  search down to find the next incrementable element
      do 5000 i=1,nlevels
c            write(*,*) "NEXTPRJ 5000"
         if ( nstack .gt. 1 .and. stack(nstack) .gt. maxprj-i ) then
            nstack = nstack-1
            if ( nstack .lt. 0 ) then
               NextIndex = .FALSE.
               return
            else
               go to 5000
            endif
         else
c            write(*,*) "NextIndex A", i,nstack,stack(nstack)
            stack(nstack) = stack(nstack) + 1
            if( nstack .eq. 1 .and. stack(1) .gt. maxprj-nlevels+1 )then
               NextIndex = .FALSE.
c            write(*,*) "NextIndex B", i,nstack,stack(nstack)
               return
            endif
c            write(*,*) "NextIndex C", i,nstack, "   stack  ", stack
            do 3000 k=nstack+1,nlevels
               stack(k) = stack(k-1)+1
c            write(*,*) "NextIndex D", i,nstack, "   stack  ", stack
 3000       continue
            NextIndex = .TRUE.
            return
         endif
 5000 continue
      NextIndex = .TRUE.
      end

C**********************************************************************C
      real*8 function G6SIGN( a, s )
      implicit none
      real*8 a, s
C----------------------------------------------------------------------C
      if ( s .gt. 1.e-20 ) then
         G6SIGN = a
      else
         G6SIGN = -a
      endif
      end

C**********************************************************************C
      subroutine vsumn( n, v1,v2,v3 )
      implicit none
      integer n, i
      real*8 v1(n), v2(n), v3(n)
C----------------------------------------------------------------------C
      do 1000 i=1,n
         v3(i) = v1(i)+v2(i)
 1000 continue
      end

C**********************************************************************C
      subroutine vdifn( n, v1,v2,v3 )
      implicit none
      integer n, i
      real*8 v1(n), v2(n), v3(n)
C----------------------------------------------------------------------C
      do 1000 i=1,n
         v3(i) = v1(i)-v2(i)
 1000 continue
      end

C**********************************************************************C
      subroutine ivdifn( n, v1,v2,v3 )
      implicit none
      integer n, i
      integer v1(n), v2(n), v3(n)
C----------------------------------------------------------------------C
      do 1000 i=1,n
          v3(i) = v1(i)-v2(i)
 1000 continue
      end

C**********************************************************************C
      subroutine  v2Cart(v ,v1,v2,v3 )
      implicit none
      real*8 v(6),v1(3),v2(3),v3(3),cell(6)
      real*8 mat8, SinAl, SinBe, SinGA, Root, Sqr 
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


C**********************************************************************C
      SUBROUTINE CONMAT(AMAT,X,BMAT)
      implicit none
      real*8  X,AMAT(9),BMAT(9)
      integer I
C----------------------------------------------------------------------C
      DO 1000 I=1,9
      BMAT(I) = X * AMAT(I)
 1000 CONTINUE
      END


C**********************************************************************C
      SUBROUTINE CROSS (X,Y,Z)
      implicit none
C----COMPUTE Z = X CROSS Y
      real*8 X(3),Y(3),Z(3)
C----------------------------------------------------------------------C
      Z(1) = X(2)*Y(3) - Y(2)*X(3)
      Z(2) =-X(1)*Y(3) + Y(1)*X(3)
      Z(3) = X(1)*Y(2) - Y(1)*X(2)
      END

C**********************************************************************C
      FUNCTION DET(A)
      implicit none
C----RETURN THE VALUE OF THE DETERMINANT OF A MATRIX
      REAL*8 X(3)
      REAL*8 DET, DOT, A(9)
C----------------------------------------------------------------------C
      CALL CROSS (A(1),A(4),X)
      DET = DOT (X,A(7))
      END

C**********************************************************************C
      FUNCTION DOT (X,Y)
      implicit none
C----COMPUTE AND RETURN THE DOT PRODUCT OF X AND Y
      real*8 DOT, X(3),Y(3)
      integer i
C----------------------------------------------------------------------C
      DOT = 0.0
      DO 1000 I=1,3
      DOT = DOT + X(I) * Y(I)
 1000 CONTINUE
      END

C**********************************************************************C
      SUBROUTINE INVER (A,B)
      implicit none
C----INVERT A THREE BY THREE MATRIX
      REAL*8 X(9)
      REAL*8 A(9),B(9), DET, DETA
      INTEGER IDATA1(3),IDATA2(3), I, J, ID1, ID2
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
         DETA = DSIGN(1.0D-20,DETA)
      ELSE
         DETA = 1.0D0 / DETA
      ENDIF

      CALL CONMAT (X,DETA,X)
      CALL TRNSPZ(X,B)
      END


C**********************************************************************C
      SUBROUTINE TRNSPZ (A,B)
      implicit none
C----PUT THE TRANSPOSE OF A INTO B
      real*8 A(9), B(9)
      integer I, J
C----------------------------------------------------------------------C
      J = 0
      DO 1000 I=1,9,3
      J = J + 1
      CALL UNVEC (A(I),B(J),B(J+3),B(J+6))
 1000 CONTINUE
      END


C**********************************************************************C
      subroutine SetStack( vp, nstack, stack )
      implicit none
      integer vp(8), stack(*)
      integer i, nstack
C----------------------------------------------------------------------C
      nstack = 0

      do 1000 i=1,8
         if ( vp(i) .le. 0 ) go to 9000
         nstack = nstack + 1
         stack(nstack) = vp(i)
 1000 continue
 9000 continue
      end

C**********************************************************************C
      subroutine Squaring( n, prjout )
      implicit none
      real*8 prjout(36)
      real*8 m1(36)
      integer i, n
C perform n squarings of the input matrix
C This uses the squared power method
C----------------------------------------------------------------------C
      do 1000 i=1,n
         call mm6( prjout, prjout, m1 )
         call cpyvn( 36, m1, prjout )
 1000 continue
      end

C**********************************************************************C
      subroutine unitmx(n,a)
      implicit none
      real*8 a(n,n)
      integer i,j,n
C----------------------------------------------------------------------C
      do 1000 i=1,n
      do 1000 j=1,n
 1000 a(i,j) = 0.0D0
      do 2000 i=1,n
 2000 a(i,i) = 1.0D0
      end



C***********************************************************************
      SUBROUTINE mm6(M1,M2,M3)
      implicit none

      real*8 M1(36),M2(36),M3(36)
      integer i36,j,k,irow,icol
C-----------------------------------------------------------------------
      do 900 i36=1,36
  900 m3(i36) = 0.0

      DO 3000 I36=1,36
      IROW = (I36+5)/6
      ICOL = MOD(I36-1,6)+1
      K = ICOL-6
      DO 2000 J=6*IROW-5,6*IROW
        K = K + 6
        M3(I36) = M3(I36) + M1(J)*M2(K)
 2000 CONTINUE
 3000 CONTINUE
      END


C***********************************************************************
      SUBROUTINE mm3(M1,M2,M3)
      implicit none

      real*8 M1(9),M2(9),M3(9)
      integer i9,j,k,irow,icol
C-----------------------------------------------------------------------
      do 900 i9=1,9
  900 m3(i9) = 0.0

      DO 3000 I9=1,9
      IROW = (I9+2)/3
      ICOL = MOD(I9-1,3)+1
      K = ICOL-3
      DO 2000 J=3*IROW-2,3*IROW
        K = K + 3
        M3(I9) = M3(I9) + M1(J)*M2(K)
 2000 CONTINUE
 3000 CONTINUE
      END


C***********************************************************************
      SUBROUTINE rowofm3(M,IROW,R)
      implicit none

      real*8 M(9),R(3)
      integer IROW,ICOL
C-----------------------------------------------------------------------
      DO ICOL=1,3
        R(ICOL) = M(ICOL+3*IROW-3)
      ENDDO
      END


C***********************************************************************
      SUBROUTINE colofm3(M,ICOL,C)
      implicit none

      real*8 M(9),C(3)
      integer IROW,ICOL
C-----------------------------------------------------------------------
      DO IROW=1,3
        C(IROW) = M(ICOL+3*IROW-3)
      ENDDO
      END


C***********************************************************************
      SUBROUTINE m3rows(M,R1,R2,R3)
      implicit none

      real*8 M(9),R1(3),R2(3),R3(3)
      integer ICOL
C-----------------------------------------------------------------------
      DO ICOL=1,3
        M(ICOL)   = R1(ICOL)
        M(ICOL+3) = R2(ICOL)
        M(ICOL+6) = R3(ICOL)
      ENDDO
      END

C***********************************************************************
      SUBROUTINE m3cols(M,C1,C2,C3)
      implicit none

      real*8 M(9),C1(3),C2(3),C3(3)
      integer IROW
C-----------------------------------------------------------------------
      DO IROW=1,3
        M(3*IROW-2) = C1(IROW)
        M(3*IROW-1) = C2(IROW)
        M(3*IROW  ) = C3(IROW)
      ENDDO
      END

C**********************************************************************C
      SUBROUTINE UNVEC (X,F,G,H)
C----RETURN THE VECTOR COMPONENTS AS SCALARS
      implicit none

      real*8 X(3),F,G,H
C----------------------------------------------------------------------C
      F = X(1)
      G = X(2)
      H = X(3)
      END

C**********************************************************************C
      function root (a)
      implicit none

      real*8 root, a
C----------------------------------------------------------------------C
      root = DSQRT( MAX( 0.D0,a))
      end

C**********************************************************************C
      function sqr (a)
      implicit none

      real*8 sqr,a
C----------------------------------------------------------------------C
      sqr = a*a
      end



