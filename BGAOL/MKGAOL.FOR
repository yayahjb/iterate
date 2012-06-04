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
      BLOCK DATA PRJS
      real*8 prj(36,21), prjperp(36,21)
      common /projectors/ prj,prjperp

C The projectors for the 15 base types (5-D boundaries
C in G6), plus a few extra for internal boundaries
C     case1
      data (prj(i,1),i=1,36) /
     * .5D0,.5D0,0.D0,0.D0,0.D0,0.D0,
     * .5D0,.5D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

C     case2
      data (prj(i,2),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,.5D0,.5D0,0.D0,0.D0,0.D0,
     * 0.D0,.5D0,.5D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

C     case3
      data (prj(i,3),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

C     case4
      data (prj(i,4),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

C     case5
      data (prj(i,5),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /

C     case6
      data (prj(i,6),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 1.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 1.D0 /

C     case7
      data (prj(i,7),i=1,36) /
     * 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0, .5D0, 0.D0, .5D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0, .5D0, 0.D0, .5D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /

C     case8
      data (prj(i,8),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0, -.5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0, -.5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 1.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 1.D0 /

C     case9
      data (prj(i,9),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /

C     caseA
      data (prj(i,10),i=1,36) /
     * .5D0, 0.D0, 0.D0, 0.D0,.5D0, 0.D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     * .5D0, 0.D0, 0.D0, 0.D0,.5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /

C     caseB
      data (prj(i,11),i=1,36) /
     *   .5D0, 0.D0, 0.D0, 0.D0, -.5D0, 0.D0,
     *   0.D0, 1.D0, 0.D0, 0.D0,  0.D0, 0.D0,
     *   0.D0, 0.D0, 1.D0, 0.D0,  0.D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 1.D0,  0.D0, 0.D0,
     *  -.5D0, 0.D0, 0.D0, 0.D0,  .5D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0,  0.D0, 1.D0 /

C     caseC
      data (prj(i,12),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /

C     caseD
      data (prj(i,13),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /

C     caseE
      data (prj(i,14),i=1,36) /
     *   .5D0, 0.D0, 0.D0, 0.D0, 0.D0, -.5D0,
     *   0.D0, 1.D0, 0.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 1.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 1.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0, 1.D0,  0.D0,
     *  -.5D0, 0.D0, 0.D0, 0.D0, 0.D0,  .5D0 /

C     caseF
      data (prj(i,15),i=1,36) /
     *  .8D0,-.2D0, 0.D0,-.2D0,-.2D0,-.2D0,
     * -.2D0, .8D0, 0.D0,-.2D0,-.2D0,-.2D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     * -.2D0,-.2D0, 0.D0, .8D0,-.2D0,-.2D0,
     * -.2D0,-.2D0, 0.D0,-.2D0, .8D0,-.2D0,
     * -.2D0,-.2D0, 0.D0,-.2D0,-.2D0, .8D0 /

C     case678X g2=g4=g5=g6=0
      data (prj(i,16),i=1,36) /
     * 1.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 1.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0 /

C     case9ABCDEX g1=g4=g5=g6=0
      data (prj(i,17),i=1,36) /
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0 /

C     caseFX g1=g2=g4=g5=g6=0
      data (prj(i,18),i=1,36) /
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0 /

C     DATA boundary67 /
      data (prj(i,19),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0,  0.D0,  0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0,  0.D0,  0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0,  .5D0,  .5D0,
     * 0.D0,  0.D0, 0.D0,  0.D0,  .5D0,  .5D0 /

C     DATA boundary9A /
      data (prj(i,20),i=1,36) /
     *  .5D0, 0.D0, 0.D0,  0.D0,  .5D0,  0.D0,
     *  0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  0.D0,  .5D0,
     *  .5D0, 0.D0, 0.D0,  0.D0,  .5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  0.D0,  .5D0 /

C     DATA boundaryCD /
      data (prj(i,21),i=1,36) /
     *  .5D0, 0.D0, 0.D0,  0.D0,  0.D0,  .5D0,
     *  0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  .5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  .5D0,  0.D0,
     *  .5D0, 0.D0, 0.D0,  0.D0,  0.D0, .5D0 /



      END

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
      real*8 VDIST(NVMAX)
      real*8 NEWGE(6)
      real*8 VT(6)
      real*8 dmin, gerr, grmin
      real*8 XDOTVN, NCDIST
      integer IVB(NVMAX)
      integer i, icycle, irefl, nv1, nexamined, id, jj
      integer irefl2
      integer nearst
      INTEGER MS(36,15)
      INTEGER ISEED(15)
      real*8 SEEDS(6,15)
      real*8 NEWSEEDS(6,15,15)
      real*8 SEEDDIST(15),NEWSEEDDIST(15,15)
      LOGICAL GOODSEED(15),NEWGOODSEED(15,15)
      LOGICAL NEARRED

      real*8 prj(36,21), prjperp(36,21)
      common /projectors/ prj,prjperp
      
C
C    The following matrices are the transofrmation
C    matrices that may be applied at the associated
C    boundaries
C
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
     5 0,-2, 0, 1, 0, 0,
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

C    M_C (g1 = g6, g4 >= g5, b -> -b, b -> b - a)
      DATA (MS(I,12),I=1,36) /
     2 1, 0, 0, 0, 0, 0,
     3 1, 1, 0, 0, 0,-1,
     4 0, 0, 1, 0, 0, 0,
     5 0, 0, 0,-1, 1, 0,
     6 0, 0, 0, 0,-1, 0,
     7-2, 0, 0, 0, 0, 1 /

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
     

     
C----------------------------------------------------------------------C

      IF (TEST .NE. 'MKGAOL') THEN
         WRITE (*,*) ' TEST WAS WRONG IN MKGAOL'
         STOP
      ENDIF
      
C     Make the perpendicular projectors
C     Each of the PERP projectors takes a g6
C     vector to the space othgonal to the
C     corresponding PRJ.  The length of that
C     vector is therefore the distance to that
C     polytope
      
      do i = 1,21
        call MKPERP(PRJ(1,i),PRJPERP(1,i))
        if (i.eq.2) THEN
C          write(*,'(A,6(6f8.3/))')
C     *      "PRJPERP(2)",(PRJPERP(jj,i),jj=1,36)
        endif
      enddo
      
C     Generate seed points on each of the 15 boundaries      
      
      call GENSEEDS(GRED,GE,SEEDS,SEEDDIST,
     * GOODSEED, PRJ, PRJPERP)
     
      do i=1,MXTREE
        Tree(i) = 0
      enddo

      NV = 0
      gerr = DSQRT(XDOTVN(6,GRED,GRED))*RATIO
      DMIN = min(0.05,0.1*DSQRT(XDOTVN(6,GRED,GRED))*RATIO)
      grmin = gred(1)**2+gred(2)**2+gred(3)**2
      WRITE (*,*) ' DMIN IN MKGOAL ',DMIN
      
C     Populate the tree with the original cell and
C     the nearby neighbors on the boundary
      NV = NV+1
C      write(*,*) "calling bldtre"
      CALL BLDTRE (MXTREE,GRED,NV*16,TREE,'BLDTRE')
      CALL CPYVN(6,GRED,V(1,NV))
      nexamined = NV
      IVB(NV) = 0
      VDIST(NV) = 0.D0
C      write(*,*)"SEED DIST",NV,SEEDDIST(NV)
C      call printg6('Storing ',GRED)
C
C     Add the good seeds that are not duplicates to
C     the tree, and for each one apply the matrix
C     for that boundary and if the image is nearly
C     reduced, add that as well.
C
      do i = 1,15
        if (GOODSEED(i)) then
C        write(*,'(a,Z2,6f9.2)')
C     *   "I, SEEDS(I) =",I, SEEDS(1,i),SEEDS(2,i),
C     *   SEEDS(3,i), SEEDS(4,i), SEEDS(5,i), SEEDS(6,i)

          NV1 = NV + 1
          IF (NEARST(MXTREE,SEEDS(1,i),DMIN,TREE,NV1,ID,
     *      'NEARST') .EQ. 0) THEN
            NV = NV + 1
            CALL BLDTRE (MXTREE,SEEDS(1,i),NV*16+i,TREE,'BLDTRE')
            CALL CPYVN(6,SEEDS(1,i),V(1,NV))
            IVB(NV) = 1
            VDIST(NV) = SEEDDIST(i)
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, SEED DIST",NV,VDIST(NV)
C            call printg6('Store SEED ',SEEDS(1,i))
          ENDIF
          call imv6(SEEDS(1,i),MS(1,i),vt)
          NV1 = NV+1
          if (NCDIST(vt,GRED).LT.gerr .and.
     *       NEARRED(vt,1.0D-6,'NEARRED') .AND.
     *       NEARST(MXTREE,VT,DMIN,TREE,NV1,ID,
     *      'NEARST') .EQ. 0) THEN
            NV = NV+1
            CALL BLDTRE (MXTREE,VT,NV*16+i,TREE,'BLDTRE')
            CALL CPYVN(6,VT,V(1,NV))
            IVB(NV) = 1
            VDIST(NV) = SEEDDIST(i)
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, M*SEED DIST",NV,VDIST(NV)
C            call printg6('Store M*SEED ',VT)
          ENDIF
C        else
C           call printg6('Reject SEED ',SEEDS(1,i))
        endif
      enddo

      DO I = 1,6
        NEWGE(I) = max(GE(I),DMIN)
      enddo

C
C     Now we have seeds to examine from nexamined+1
C     through NV.  In the course of doing so, NV may
C     increase.

1000  nexamined = nexamined+1
      if (nexamined .gt. NV .OR. NV.GT.NVMAX-31) THEN
         WRITE (*,*) ' NV IN MKREFL ',NV
         return
      ENDIF
      if (ivb(nexamined).gt.10) go to 1000

      DO I = 1,6
        NEWGE(I) = max(GE(I)/(1.414D0**IVB(nexamined)),DMIN)
      ENDDO

C     Generate seed points on each of the 15 boundaries      

      call GENSEEDS(V(1,nexamined),NEWGE,SEEDS,SEEDDIST,
     * GOODSEED, PRJ, PRJPERP)
C
C     Add the good seeds that are not duplicates to
C     the tree, and for each one apply the matrix
C     for that boundary and if the image is nearly
C     reduced, add that as well.
C
C     Differs from the prior loop only in the distance
C     calculation
C
      do i = 1,15
        if (GOODSEED(i) .and.
     *    SEEDS(1,i)**2 + SEEDS(2,i)**2 + SEEDS(3,i)**2
     *    .lt. 10.0D0*grmin) then
          NV1 = NV + 1
          IF (NEARST(MXTREE,SEEDS(1,i),DMIN,TREE,NV1,ID,
     *      'NEARST') .EQ. 0) THEN
            NV = NV + 1
            CALL BLDTRE (MXTREE,SEEDS(1,i),NV*16+i,TREE,'BLDTRE')
            CALL CPYVN(6,SEEDS(1,i),V(1,NV))
            IVB(NV) = IVB(nexamined)+1
            VDIST(NV) = dsqrt(SEEDDIST(i)**2+
     *        VDIST(nexamined)**2)
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, SEED DIST",NV,VDIST(NV)
C            call printg6('Store SEED ',SEEDS(1,i))
          ENDIF
          call imv6(SEEDS(1,i),MS(1,i),vt)
          NV1 = NV+1
          if (NCDIST(vt,GRED).LT.gerr .and.
     *       NEARRED(vt,1.0D-6,'NEARRED') .AND.
     *       NEARST(MXTREE,VT,DMIN,TREE,NV1,ID,
     *      'NEARST') .EQ. 0) THEN
            NV = NV+1
            CALL BLDTRE (MXTREE,VT,NV*16+i,TREE,'BLDTRE')
            CALL CPYVN(6,VT,V(1,NV))
            IVB(NV) = IVB(nexamined)+1
            VDIST(NV) = dsqrt(SEEDDIST(i)**2+
     *        VDIST(nexamined)**2)
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, M*SEED DIST",NV,VDIST(NV)
C            call printg6('Store M*SEED ',VT)
          ENDIF
C        else
C           call printg6('Reject SEED ',SEEDS(1,i))
        endif
      enddo
      go to 1000
      end
      
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
C      WRITE(*,*) "PRJ",PRJ
C      WRITE(*,*) "PRJPERP",PRJPERP
      RETURN
      END
      
      SUBROUTINE PRINTG6(TEXT,VECTOR)
      implicit none
      character*(*) TEXT
      real*8 VECTOR(6),FVECTOR(6)
      WRITE(*,'(A,1x,6F9.2)')TEXT,VECTOR

      CALL FOLDBED(VECTOR,FVECTOR,0.05D1)
      WRITE(*,'(A,1x,6F9.2)')"FOLDED: ",FVECTOR
      RETURN
      END


      SUBROUTINE GENSEEDS( GVEC, GE, SEEDS, SEEDDIST,
     * GOODSEED, PRJ, PRJPERP)
C
C     GENSEEDS TAKES AN ARBiTRARY G6 VECTOR AND PROJECTS
C     IT ONTO EACH OF THE 15 BOUNDARIES.  EVEN IF THE
C     ORIGINAL VECTOR WAS REDUCED THE PROJECTED VECTOR
C     MAY FAIL TO BE REDUCED, AND MAY EVEN FAIL TO BE
C     NEARLY REDUCED.  IT WILL BE MARKED AS A GOOD SEED
C     PROVIDED IT IS BOTH NEARLY REDUCED AND IS WITHIN
C     3.5 TIMES THE ERROR BOX GE OF GVEC.
C
      implicit none
      real*8 GVEC(6), GE(6)
      real*8 PRJ(36,21),PRJPERP(36,21)
      real*8 SEEDS(6,15), SEEDDIST(15)
      logical GOODSEED(15)
      character*20 SEEDLAB
      real*8 vtemp(6)
      real*8 boundary67
      real*8 boundary9A
      real*8 boundaryCD
      real*8 a
      integer i, ip
      real*8 anorm
      logical NEARRED
C----------------------------------------------------------------------C

C     Compute the seeds and the seed distances

C      write(*,*) "GENSEEDS GVEC, GE"
C      write(*,*) GVEC
C      write(*,*) GE

      call rmv6( GVEC, PRJPERP(1,19), vtemp )
      boundary67 = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,20), vtemp )
      boundary9A = anorm(6,vtemp)
      call rmv6( GVEC, PRJPERP(1,21), vtemp )
      boundaryCD = anorm(6,vtemp)

C      call printg6('GENSEEDS GVEC ',GVEC)

      do ip=1,15
        call rmv6( GVEC, PRJ(1,ip), SEEDS(1,ip))
        call rmv6( GVEC, PRJPERP(1,ip), vtemp )
C        write(*,'(A,i2/,6(6F9.3/))'),"PRJ",
C     *  ip,(PRJ(i,ip),i=1,36)
C        write(*,'(A,i2/,6(6F9.3/))'),"PRJPERP",
C     *  ip,(PRJPERP(i,ip),i=1,36)
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
        GOODSEED(ip) = NEARRED(SEEDS(1,ip),1.0D-6,'NEARRED')
        if (.NOT.GOODSEED(ip)) then
C         call printg6("Non-reduced projection",seeds(1,ip))
C         return
        endif
        do i = 1,6
          if (abs(GVEC(i)-SEEDS(i,ip)) .gt.
     *      3.5D0*GE(i)) GOODSEED(ip) = .false.
        enddo
C        write(seedlab,'(''BDRY'', I2)')IP
C        call printg6(seedlab,seeds(1,ip))
      enddo
      
      return
      end

C
C     Compute XS (distances to boundary sets) and
C             YS (signed distances along boundaries)
C

      SUBROUTINE BDCOORD(GVEC,XS,YS)
      implicit none
      real*8 sqrt2, sqrt5
      real*8 XS(6),YS(6)
      real*8 GVEC(6)
C
C     The boundaries are in the order
C       1, 2, 678, 9AB, CDE, F
      sqrt2 = SQRT(2.D0)
      sqrt5 = SQRT(5.D0)
      XS(1) = (GVEC(2) - GVEC(1))/sqrt2
      XS(2) = (GVEC(3) - GVEC(2))/sqrt2
      XS(3) = (GVEC(2) - ABS(GVEC(4)))/sqrt2
      XS(4) = (GVEC(1) - ABS(GVEC(5)))/sqrt2
      XS(5) = (GVEC(1) - ABS(GVEC(6)))/sqrt2
      XS(6) = (GVEC(1)+GVEC(2) +
     *  GVEC(4) + GVEC(5) + GVEC(6))/sqrt5
      YS(1) = (ABS(GVEC(4)) - ABS(GVEC(5)))/sqrt2
      YS(2) = (ABS(GVEC(5)) - ABS(GVEC(6)))/sqrt2
      YS(3) = (GVEC(5) - ABS(GVEC(6)/2.D0))/sqrt2
      YS(4) = (GVEC(4) - ABS(GVEC(6)/2.D0))/sqrt2
      YS(5) = (GVEC(4) - ABS(GVEC(5)/2.D0))/sqrt2
      YS(6) = (GVEC(2) - GVEC(1) + GVEC(4) - GVEC(5))
     *  /2.D0
      END


C
C     Compute the distance reflecting from a boundary
C
      real*8 function mirdist(h1,h2,db)
      implicit none
      real*8 h1, h2, db
      mirdist = sqrt((h1+h2)**2+db**2)
      end

C
C     Compute the minimal distance between two Niggli-reduced
C     vectors in the Niggli Cone
C
C
      REAL*8 FUNCTION NCDIST(gv1,gv2)
      implicit none
      real*8 gv1(6), gv2(6), pgv1(6,9), pgv2(6,9)
      real*8 fpgv1(6,9), fpgv2(6,9)
      real*8 anorm, vtemp(6)
      real*8 prj(36,21), prjperp(36,21)
      integer ind, jnd, ip, bindex(9)
      real*8 vdist(9,3),h1,h2,db, mirdist
      logical valid1(9),valid2(9),isneg(9),ispos(9)
      common /projectors/ prj,prjperp
      data bindex/1,2,6,8,9,11,12,14,15/
      data isneg/.true.,.true.,.false.,.true.,
     *  .false.,.true.,.false.,.true.,.true./
      data ispos/.true.,.true.,.true.,.false.,
     *  .true.,.false.,.true.,.false.,.false./

C      write(*,*) "gv1 ",gv1
C      write(*,*) "gv2 ",gv2
      call vdifn(6,gv1,gv2,vtemp)
      NCDIST = anorm(6,vtemp)
      do ind=1,9
C        write(*,*) "ind ", ind
        ip = bindex(ind)
        valid1(ind) = .false.
        valid2(ind) = .false.
        if (isneg(ind)) then
          if (gv1(4)*gv1(5)*gv1(6).le.1.e-6) 
     *      valid1(ind) = .true.
          if (gv2(4)*gv2(5)*gv2(6).le.1.e-6)
     *      valid2(ind) = .true.
        endif
        if (ispos(ind)) then
          if (gv1(4)*gv1(5)*gv1(6).gt.1.e-6)
     *      valid1(ind) = .true.
          if (gv2(4)*gv2(5)*gv2(6).gt.1.e-6)
     *      valid2(ind) = .true.
        endif
        if (valid2(ind)) then
          call rmv6(gv2,PRJ(1,ip),pgv2(1,ind))
          call vdifn(6,gv2,pgv2(1,ind),vtemp)
          vdist(ind,2) = anorm(6,vtemp)
          call foldbed(pgv2(1,ind),fpgv2(1,ind),1.0D0)
          if (ind.eq.8 .or. ind.eq.11 .or. ind.eq.14) then
            jnd = ind-2
            if (valid1(jnd)) then
            call vdifn(6,fpgv1(1,jnd),fpgv2(1,ind),vtemp)
            vdist(ind,3) = anorm(6,vtemp)**2
            h1 = vdist(jnd,1)
            h2 = vdist(ind,2)
            db = vdist(ind,3)
            if (mirdist(h1,h2,db).lt.NCDIST)
     *        NCDIST = mirdist(h1,h2,db)
            endif
          endif
        endif
        if (valid1(ind)) then
          call rmv6(gv1,PRJ(1,ip),pgv1(1,ind))
          call vdifn(6,gv1,pgv1(1,ind),vtemp)
          vdist(ind,1) = anorm(6,vtemp)
          call foldbed(pgv1(1,ind),fpgv1(1,ind),1.0D0)
          if (valid2(ind)) then
            call vdifn(6,fpgv1(1,ind),fpgv2(1,ind),vtemp)
            vdist(ind,3) = anorm(6,vtemp)**2
            h1 = vdist(ind,1)
            h2 = vdist(ind,2)
            db = vdist(ind,3)
            if (mirdist(h1,h2,db).lt.NCDIST) 
     *        NCDIST=mirdist(h1,h2,db)
          endif
          if (ind.eq.8 .or. ind.eq.11 .or. ind.eq.14) then
            jnd = ind-2
            if (valid2(jnd)) then
            call vdifn(6,fpgv1(1,ind),fpgv2(1,jnd),vtemp)
            vdist(ind,3) = anorm(6,vtemp)**2
            h1 = vdist(ind,1)
            h2 = vdist(jnd,2)
            db = vdist(ind,3)
            if (mirdist(h1,h2,db).lt.NCDIST) 
     *        NCDIST=mirdist(h1,h2,db)
            endif
          endif
        endif
      enddo
      end



C
C     FOLDBED -- fold the boundaries where appropriate
C     taking a G6 vector GVEC to a folded G6 vector
C     FGVEC.  
C     All folds try to go as much as possible to --- to
C     allow for interactions with the F boundary, and
C     to |g4|>=|g5|>=|g6| when possible to avoid infinite
C     loops.
C
C
      SUBROUTINE FOLDBED(GVEC,FGVEC,DELTA)
      implicit none
      real*8 GVEC(6), FGVEC(6), DELTA 
      real*8 sqrt2, sqrt5, temp

      integer i
      logical folded
      integer folds
      
      do i = 1,6
      FGVEC(i) = GVEC(i)
      enddo

      sqrt2 = sqrt(2.D0)
      sqrt5 = sqrt(5.D0)
      folds = 0

 100  continue
      folded = .false.
C      write(*,*)"folding",FGVEC
C      print *,"folds ",folds

C
C     fold the 678 boundary along g4 around |g6|/2
C
      if ((FGVEC(2) - ABS(FGVEC(4))).LE.sqrt2*DELTA
     *  .AND. 2.D0*FGVEC(5).GT.ABS(FGVEC(6))) THEN
        FGVEC(5) = MIN(FGVEC(5),ABS(FGVEC(6))-FGVEC(5))
        IF(FGVEC(5).LE.0.D0) THEN
          FGVEC(4) = -ABS(FGVEC(4))
          FGVEC(6) = -ABS(FGVEC(6))
        ENDIF
        folded = .true.
C        print *,"folding on 678"

      ENDIF


C
C     fold the 9AB boundary along g4 around |g6|/2
C
      if ((FGVEC(1) - ABS(FGVEC(5))).LE.sqrt2*DELTA
     *  .AND. 2.D0*FGVEC(4).GT.ABS(FGVEC(6))) THEN
        FGVEC(4) = MIN(FGVEC(4),ABS(FGVEC(6))-FGVEC(4))
        IF(FGVEC(4).LE.0.D0) THEN
          FGVEC(5) = -ABS(FGVEC(5))
          FGVEC(6) = -ABS(FGVEC(6))
        ENDIF
        folded = .true.
C        print *,"folding on 9AB"

      ENDIF

C
C     fold the CDE boundary along g4 around |g5|/2
C
      if ((FGVEC(1) - ABS(FGVEC(6))).LE.sqrt2*DELTA
     *  .AND. 2.D0*FGVEC(4).GT.ABS(FGVEC(5))) THEN
        FGVEC(4) = MIN(FGVEC(4),ABS(FGVEC(5))-FGVEC(4))
        IF(FGVEC(4).LE.0.D0) THEN
          FGVEC(5) = -ABS(FGVEC(5))
          FGVEC(6) = -ABS(FGVEC(6))
        ENDIF
        folded = .true.
C        print *,"folding on CDE"

      ENDIF


C
C     fold the F boundary along 
C        [0,0,0,1/sqrt(2),-1/sqrt(2),0] around
C        [g1,g2,g3,
C           g4+(g1-g2-g4+g5)/2,g5-(g1-g2-g4+g5)/2,g6]
C     by swapping g4 and g5
C
      if (FGVEC(1)+FGVEC(1)
     *  +FGVEC(4)+FGVEC(5)+FGVEC(6).LE.sqrt5*DELTA
     *  .AND. FGVEC(1)-FGVEC(2)-FGVEC(4)+FGVEC(5)
     *        .LT.0.D0) THEN
        
        temp = FGVEC(4)
        FGVEC(4) = FGVEC(5)
        FGVEC(5) = temp
        folded = .true.
C        print *,"folding on F"
      ENDIF

C
C     fold the 1 boundary along 
C        [0,0,0,1/sqrt(2),-1/sqrt(2),0] around
C        [g1,g2,g3,
C           (g4+g5)/2,(g4+g5)/2,g6]
C     by swapping g4 and g5
C
      if (ABS(FGVEC(2)-FGVEC(1)).LE.DELTA
     *  .AND. ABS(FGVEC(4)).LT.ABS(FGVEC(5))) THEN
        
C        print *,"folding on 1", FGVEC(4), FGVEC(5)
        temp = FGVEC(4)
        FGVEC(4) = FGVEC(5)
        FGVEC(5) = temp
        folded = .true.
C        print *,"folded on 1", FGVEC(4), FGVEC(5)

      ENDIF

C
C     fold the 2 boundary along 
C        [0,0,0,0,1/sqrt(2),-1/sqrt(2)] around
C        [g1,g2,g3,g4
C           (g5+g6)/2,(g5+g6)/2]
C     by swapping g5 and g6
C
      if (ABS(FGVEC(3)-FGVEC(2)).LE.DELTA
     *  .AND. ABS(FGVEC(5)).LT.ABS(FGVEC(6))) THEN
        
        temp = FGVEC(5)
        FGVEC(5) = FGVEC(6)
        FGVEC(6) = temp
        folded = .true.
C        print *,"folding on 2"


      ENDIF

      folds = folds+1
      if (folded.and.folds.LT.5) go to 100

      END




