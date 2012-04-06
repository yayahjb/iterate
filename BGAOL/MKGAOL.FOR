
C     MKGAOL -- Make GAOL (make lots of cells)
C     a subroutine of
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


C      CALL MKGAOL (MXTREE,TREE,NVMAX,V,NV,GRED,GE,
C    *  RATIO, VDIST, IVB,'MKGAOL')

C**********************************************************************C
      SUBROUTINE MKGAOL
     * (MXTREE,TREE,NVMAX,V,NV,GRED,GE,RATIO,
     *  VDIST,IVB,TEST)
      implicit none
      integer mxswtc, mxtree, nvmax, nv
      real*8 ratio
      PARAMETER (MXSWTC=3)
      CHARACTER *6 TEST
      real*8 TREE(MXTREE), V(6,NVMAX), GRED(6), GE(6)
      real*8 VDIST(MXTREE)
      real*8 NEWGE(6)
      real*8 VT(6)
      real*8 dmin, grmin
      real*8 dotvn
      integer IVB(MXTREE)
      integer i, icycle, irefl, nv1, id
      integer irefl2
      integer nearst
      INTEGER MS(36,15)
      INTEGER ISEED(15)
      real*8 SEEDS(6,15)
      real*8 NEWSEEDS(6,15,15)
      real*8 SEEDDIST(15),NEWSEEDDIST(15)
      LOGICAL GOODSEED(15),NEWGOODSEED(15,15)
      real*8 prj(36,21), prjperp(36,21)
      
C    M_1 (g1 = g2, a -> b, b -> a)
      DATA (MS(I,1), I=1,36) /
     2 0,1,0,0,0,0,
     3 1,0,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,0,1,0,
     6 0,0,0,1,0,0,
     7 0,0,0,0,0,1 /

C    M_2 (g2 = g3, b -> c, c -> b)
      DATA (MS(I,2), I=1,36) /
     2 1,0,0,0,0,0,
     3 0,0,1,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,1,0,0,
     6 0,0,0,0,0,1,
     7 0,0,0,0,1,0 /

C    M_3 (g4 = 0, a -> -a)
      DATA (MS(I,3),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,1,0,0,
     6 0,0,0,0,-1,0,
     7 0,0,0,0,0,-1 /

C    M_4 (g5 = 0, b -> -b)
      DATA (MS(I,4),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,-1,0,0,
     6 0,0,0,0,1,0,
     7 0,0,0,0,0,-1 /

C    M_5 (g6 = 0, c -> -c)
      DATA (MS(I,5),I=1,36) /
     2 1,0,0,0,0,0,
     3 0,1,0,0,0,0,
     4 0,0,1,0,0,0,
     5 0,0,0,-1,0,0,
     6 0,0,0,0,-1,0,
     7 0,0,0,0,0,1 /

C    M_6 (g2 = g4, g5 >= g6, b -> -b, c -> b - c)
      DATA (MS(I,6),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 0, 1, 1,-1, 0, 0,
     5 0,-2, 0,-1, 0, 0,
     6 0, 0, 0, 0,-1, 1,
     7 0, 0, 0, 0, 0,-1 /

C    M_7 (g2 = g4, g5 < g6, c -> b - c)
      DATA (MS(I,7),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 0, 1, 1,-1, 0, 0,
     5 0, 2, 0,-1, 0, 0,
     6 0, 0, 0, 0,-1, 1,
     7 0, 0, 0, 0, 0, 1 /

C    M_8 (g2 = -g4, a -> -a, c -> b + c)
      DATA (MS(I,8),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 0, 1, 1, 1, 0, 0,
     5 0, 2, 0, 1, 0, 0,
     6 0, 0, 0, 0,-1,-1,
     7 0, 0, 0, 0, 0,-1 /

C    M_9 (g1 = g5, g4 >= g6, b -> -b, c -> c - a)
      DATA (MS(I,9),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 1, 0, 1, 0,-1, 0,
     5 0, 0, 0,-1, 0, 1,
     6-2, 0, 0, 0, 1, 0,
     7 0, 0, 0, 0, 0,-1 /

C    M_A (g1 = g5, g4 < g6, c -> a - c)
      DATA (MS(I,10),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 1, 0, 1, 0,-1, 0,
     5 0, 0, 0,-1, 0, 1,
     6 2, 0, 0, 0,-1, 0,
     7 0, 0, 0, 0, 0, 1 /

C    M_B (g1 = -g5, b -> -b, c -> a + c)
      DATA (MS(I,11),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 1, 0, 1, 0, 1, 0,
     5 0, 0, 0,-1, 0,-1,
     6 2, 0, 0, 0, 1, 0,
     7 0, 0, 0, 0, 0,-1 /

C    M_C (g1 = g6, g4 >= g6, b -> -b, b -> b - a)
      DATA (MS(I,12),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 1, 1, 0, 0, 0,-1,
     4 0, 0, 1, 0, 0, 0,
     5 0, 0, 0,-1, 1, 0,
     6 0, 0, 0, 0,-1, 0,
     7-2, 0, 0, 0, 0,-1 /

C    M_D (g1 = g6, g4 < g5, b -> a - b)
      DATA (MS(I,13),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 1, 1, 0, 0, 0,-1,
     4 0, 0, 1, 0, 0, 0,
     5 0, 0, 0,-1, 1, 0,
     6 0, 0, 0, 0, 1, 0,
     7 2, 0, 0, 0, 0,-1 /

C    M_E (g1 = -g6, b -> a + b, c -> -c )
      DATA (MS(I,14),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 1, 1, 0, 0, 0, 1,
     4 0, 0, 1, 0, 0, 0,
     5 0, 0, 0,-1,-1, 0,
     6 0, 0, 0, 0,-1, 0,
     7 2, 0, 0, 0, 0, 1 /

C    M_F (g1+g2+g3+g4+g5+g6 = g3, c -> -(a+b+c))
      DATA (MS(I,15),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 0, 1, 0, 0, 0, 0,
     4 1, 1, 1, 1, 1, 1,
     5 0,-2, 0,-1, 0,-1,
     6-2, 0, 0, 0,-1,-1,
     7 0, 0, 0, 0, 0, 1 /
     
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


     
C----------------------------------------------------------------------C

      IF (TEST .NE. 'MKGAOL') THEN
         WRITE (*,*) ' TEST WAS WRONG IN MKGAOL'
         STOP
      ENDIF
      
C     Make the perpendicular projectors
      
      do i = 1,21
        call MAKEPERP(PRJ(1,i),PRJPERP(1,i))
      enddo
      
C     Generate seed points on each of the 15 boundaries      
      
      call GENSEEDS(GRED,GE,SEEDS,SEEDDIST,
     * GOODSEED, PRJ, PRJPERP)
     
      TREE(1) = 0
      NV = 0
      DMIN = 0.1*DSQRT(DOTVN(6,GRED,GRED))*RATIO
      grmin = gred(1)**2+gred(2)**2+gred(3)**2
      WRITE (*,*) ' DMIN IN MKGOAL ',DMIN
      
C     Populate the tree with the original cell and
C     the nearby neighbors on the boundary
      NV = NV+1
      CALL BLDTRE (MXTREE,GRED,NV*16,TREE,'BLDTRE')
      do i = 1,15
        if (GOODSEED(i)) then
          NV1 = NV + 1
          IF (NEARST(MXTREE,SEEDS(1,i),DMIN,TREE,NV1,ID,
     *      'NEARST') .EQ. 0) THEN
            NV = NV + 1
            CALL BLDTRE (MXTREE,SEEDS(1,i),NV*16+i,TREE,'BLDTRE')
          ENDIF
        endif
      enddo
      
C     For each seed point that is good, we apply the
C     matching matrix, yielding a possible new
C     seed point on a different boundary.  Those
C     become the new seed points.
C
C     In each case it is important to keep the original
C     seed distances to original cell and to keep with
C     each derived cell the index of the original
C     seed for that cell, because the final distance
C     from the original cell to the projection of the
C     drived cell into a lattice charcacter is just the
C     hypoteneuse of the triangle.  The needed index is
C     held in the bottom 4 bits of the tree entry ID.
C
      DO i = 1,6
        NEWGE(i) = DMIN*1.D-4+GE(i)/1.D4
      ENDDO
      DO i = 1,15
        ISEED(i) = i
      ENDDO
      DO ICYCLE=1,MXSWTC+1
      DO IREFL=1,15
         DO i = 1,15
           DO IREFL2 = 1,15
             NEWGOODSEED(IREFL2,IREFL) = .false.
           ENDDO
         ENDDO
         IF (GOODSEED(IREFL)) THEN
C
C          If the seed is good, apply the matrix
C          for that boundary and map the boundaries
C          on which the image of the seed lands
C
           call imv6(SEEDS(1,i),MS(1,i),vt)
           call GENSEEDS(vt,NEWGE,NEWSEEDS(1,1,IREFL),
     *       NEWSEEDDIST, NEWGOODSEED(1,IREFL), 
     *       PRJ, PRJPERP)
           DO IREFL2 = 1,15
             IF (NEWGOODSEED(IREFL2,IREFL)) THEN
               IF (NEWSEEDS(1,IREFL2,IREFL)**2
     *            + NEWSEEDS(2,IREFL2,IREFL)**2
     *            + NEWSEEDS(3,IREFL2,IREFL)**2 
     *            .lt. 10.0*grmin) THEN
                 NV1 = NV + 1
                 IF(NEARST(MXTREE,NEWSEEDS(1,IREFL2,IREFL),
     2             DMIN,TREE,NV1,ID,'NEARST') .EQ. 0)
     3             THEN
                     NV = NV + 1
                     CALL BLDTRE (MXTREE,
     2                 NEWSEEDS(1,IREFL2,IREFL),
     3                 NV*16+ISEED(IREFL),TREE,'BLDTRE')
                       CALL CPYVN(6,NEWSEEDS(1,IREFL2,IREFL),
     2                   V(1,NV))
                       IVP(1,NV) = ISEED(IREFL)
                       VDIST(NV) = SEEDDIST(IVP(1,NV))
                       IF(
                   ELSE
                     NEWGOODSEED(IREFL2,IREFL) = .false.
                   ENDIF
                 ELSE
                   NEWGOODSEED(IREFL2,IREFL) = .false.
                 ENDIF
               ELSE
                 NEWGOODSEED(IREFL2,IREFL) = .false.
               ENDIF
C             ENDIF
           ENDDO
         ENDIF
      ENDDO
      DO i = 1,15
         GOODSEED(i) = .false.
         DO IREFL2 = 1,15
           IF (NEWGOODSEED(i,IREFL2)) then
             CALL CPYVN(6,NEWSEEDS(1,i,IREFL2),SEEDS(1,i))
             ISEED(i) = IREFL2
           ENDIF
         ENDDO
      ENDDO
      ENDDO
      WRITE (*,*) ' NV IN MKREFL ',NV
      END
      
      SUBROUTINE MKPERP( PRJ, PRJPERP)
      implicit none
      real*8 PRJ(36),PRJPERP(36)
      integer i
      do i = 1,36
        PRJPERP(i) = - PRJ(i)
      enddo
      do i = 1,36,7
        PRJPERP(i) = 1.D0+PRJPERP(i)
      enddo
      RETURN
      END
      
      SUBROUTINE GENSEEDS( GVEC, GE, SEEDS, SEEDDIST,
     * GOODSEED, PRJ, PRJPERP)
      implicit none
      real*8 GVEC(6), GE(6)
      real*8 PRJ(36,21),PRJPERP(36,21)
      real*8 SEEDS(6,15), SEEDDIST(15)
      logical GOODSEED(15)
      real*8 vtemp(6)
      real*8 boundary67
      real*8 boundary9A
      real*8 boundaryCD
      real*8 boundary678X
      real*8 boundary9ABCDEX
      real*8 boundaryFX
      real*8 a
      integer i, ip
      real*8 anorm
C----------------------------------------------------------------------C

C     Compute the seeds and the seed distances

      call rmv6( GVEC, PRJPERP(1,16), vtemp )
      boundary678X = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,17), vtemp )
      boundary9ABCDEX = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,18), vtemp )
      boundaryFX = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,19), vtemp )
      boundary67 = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,20), vtemp )
      boundary9A = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,21), vtemp )
      boundaryCD = anorm(6,vtemp)


      do ip=1,15
        call rmv6( GVEC, PRJ(1,ip), SEEDS(1,ip))
        call rmv6( GVEC, PRJPERP(1,ip), vtemp )
        a = anorm(6,vtemp)
        if ((GVEC(4)*GVEC(5)*GVEC(6) .le. 0.D0
     *        .and. (ip .eq. 6 .or. ip .eq. 7
     *          .or. ip .eq. 9 .or. ip .eq. 10
     *          .or. ip .eq. 12.or. ip .eq. 13))
     *        .or. (GVEC(4)*GVEC(5)*GVEC(6) .gt. 0.D0
     *        .and. (ip .eq. 8 .or. ip .eq. 11
     *          .or. ip .eq. 14. or. ip .eq. 15))) then
            if ( ip .eq. 6 .or. ip .eq. 7 .or. ip .eq. 8 ) then
              call rmv6( GVEC, PRJ(1,16), SEEDS(1,ip) )
              call rmv6( GVEC, PRJPERP(1,16), vtemp )
            endif
            if ( ip .ge. 9 .and. ip .le. 14) then
              call rmv6( GVEC, PRJ(1,17), SEEDS(1,ip) )
              call rmv6( GVEC, PRJPERP(1,17), vtemp )
            endif
            if ( ip .eq. 15) then
              call rmv6( GVEC, PRJ(1,18), SEEDS(1,ip) )
              call rmv6( GVEC, PRJPERP(1,18), vtemp )
            endif
            a = anorm(6,vtemp)
        endif
        if ( ip .eq.  6 .and. GVEC(5) .le. GVEC(6)
     *     .and. GVEC(5) .gt. 0.D0) then
           a = boundary67
        elseif ( ip .eq.  7 .and. GVEC(5) .ge. GVEC(6)
     *     .and. GVEC(6) .gt. 0.D0) then
           a = boundary67
        elseif ( ip .eq.  9 .and. GVEC(4) .le. GVEC(6)
     *     .and. GVEC(4) .gt. 0.D0 ) then
           a = boundary9A
        elseif ( ip .eq. 10 .and. GVEC(4) .ge. GVEC(6)
     *     .and. GVEC(6) .gt. 0.D0 ) then
           a = boundary9A
        elseif ( ip .eq. 12 .and. GVEC(4) .le. GVEC(5)
     *     .and. GVEC(4) .gt. 0.D0 ) then
           a = boundaryCD
        elseif ( ip .eq. 13 .and. GVEC(4) .ge. GVEC(5)
     *     .and. GVEC(5) .gt. 0.D0 ) then
           a = boundaryCD
        endif
        
C       get distance and mark as bad if outside
C       3.5 times the errorbox

        SEEDDIST(ip) = a
        GOODSEED(ip) = .true.
        do i = 1,6
          if (abs(GVEC(i)-SEEDS(i,ip)) .gt.
     *      3.5D0*GE(i)) GOODSEED(ip) = .false.
        enddo
      enddo
      
      return
      end

