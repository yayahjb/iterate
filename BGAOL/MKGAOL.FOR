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
      real*8 prj(36,22), prjperp(36,22)
      real*8 prjhat(36,15)
      integer ms(36,18)
      common /projectors/ prj,prjperp,prjhat,ms

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

       data (prjhat(i,1),i=1,36) /
     * .5D0,.5D0,0.D0,0.D0,0.D0,0.D0,
     * .5D0,.5D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,.5D0,.5D0,0.D0,
     * 0.D0,0.D0,0.D0,.5D0,.5D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

      data (prjperp(i,1),i=1,36) /
     *  .5D0,-.5D0,0.D0,0.D0,0.D0,0.D0,
     * -.5D0, .5D0,0.D0,0.D0,0.D0,0.D0,
     *  0.D0, 0.D0,0.D0,0.D0,0.D0,0.D0,
     *  0.D0, 0.D0,0.D0,0.D0,0.D0,0.D0,
     *  0.D0, 0.D0,0.D0,0.D0,0.D0,0.D0,
     *  0.D0, 0.D0,0.D0,0.D0,0.D0,0.D0 /


C     case2
      data (prj(i,2),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,.5D0,.5D0,0.D0,0.D0,0.D0,
     * 0.D0,.5D0,.5D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

      data (prjhat(i,2),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,.5D0,.5D0,0.D0,0.D0,0.D0,
     * 0.D0,.5D0,.5D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,.5D0,.5D0,
     * 0.D0,0.D0,0.D0,0.D0,.5D0,.5D0 /

      data (prjperp(i,2),i=1,36) /
     * 0.D0, 0.D0, 0.D0,0.D0,0.D0,0.D0,
     * 0.D0, .5D0,-.5D0,0.D0,0.D0,0.D0,
     * 0.D0,-.5D0, .5D0,0.D0,0.D0,0.D0,
     * 0.D0, 0.D0, 0.D0,0.D0,0.D0,0.D0,
     * 0.D0, 0.D0, 0.D0,0.D0,0.D0,0.D0,
     * 0.D0, 0.D0, 0.D0,0.D0,0.D0,0.D0 /

C     case3
      data (prj(i,3),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

      data (prjhat(i,3),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /

      data (prjperp(i,3),i=1,36) /
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /


C     case4
      data (prj(i,4),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

      data (prjhat(i,4),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /

      data (prjperp(i,4),i=1,36) /
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /

C     case5
      data (prj(i,5),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /

      data (prjhat(i,5),i=1,36) /
     * 1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,1.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0 /

      data (prjperp(i,5),i=1,36) /
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     * 0.D0,0.D0,0.D0,0.D0,0.D0,1.D0 /

C     case6
      data (prj(i,6),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 1.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 1.D0 /

      data (prjhat(i,6),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, .2D0, .4D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, .4D0, .8D0 /

      data (prjperp(i,6),i=1,36) /
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0, -.5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0, -.5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0 /

C     case7
      data (prj(i,7),i=1,36) /
     * 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0, .5D0, 0.D0, .5D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0, .5D0, 0.D0, .5D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /

      data (prjhat(i,7),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, .2D0, .4D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, .4D0, .8D0 /

      data (prjperp(i,7),i=1,36) /
     * 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0, .5D0, 0.D0,-.5D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0,-.5D0, 0.D0, .5D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /


C     case8
      data (prj(i,8),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0, -.5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0, -.5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 1.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 1.D0 /

      data (prjperp(i,8),i=1,36) /
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, 0.D0, 0.D0 /

C     case9
      data (prj(i,9),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /

      data (prjhat(i,9),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .2D0, 0.D0, .4D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .4D0, 0.D0, .8D0 /

      data (prjperp(i,9),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0,-.5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * -.5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /


C     caseA
      data (prj(i,10),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /


      data (prjhat(i,10),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .2D0, 0.D0, .4D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .4D0, 0.D0, .8D0 /

      data (prjperp(i,10),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0,-.5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * -.5D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /


C     caseB
      data (prj(i,11),i=1,36) /
     *   .5D0, 0.D0, 0.D0, 0.D0, -.5D0, 0.D0,
     *   0.D0, 1.D0, 0.D0, 0.D0,  0.D0, 0.D0,
     *   0.D0, 0.D0, 1.D0, 0.D0,  0.D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 1.D0,  0.D0, 0.D0,
     *  -.5D0, 0.D0, 0.D0, 0.D0,  .5D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0,  0.D0, 1.D0 /

      data (prjperp(i,11),i=1,36) /
     *   .5D0, 0.D0, 0.D0, 0.D0,  .5D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0,  0.D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0,  0.D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0,  0.D0, 0.D0,
     *   .5D0, 0.D0, 0.D0, 0.D0,  .5D0, 0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0,  0.D0, 0.D0 /


C     caseC
      data (prj(i,12),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /

      data (prjhat(i,12),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .2D0, .4D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .4D0, .8D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /

      data (prjperp(i,12),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0,-.5D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /


C     caseD
      data (prj(i,13),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /

      data (prjhat(i,13),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0,
     *  0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .2D0, .4D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, .4D0, .8D0, 0.D0,
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /

      data (prjperp(i,13),i=1,36) /
     *  .5D0, 0.D0, 0.D0, 0.D0, 0.D0,-.5D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     * -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0 /


C     caseE
      data (prj(i,14),i=1,36) /
     *   .5D0, 0.D0, 0.D0, 0.D0, 0.D0, -.5D0,
     *   0.D0, 1.D0, 0.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 1.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 1.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0, 1.D0,  0.D0,
     *  -.5D0, 0.D0, 0.D0, 0.D0, 0.D0,  .5D0 /

      data (prjperp(i,14),i=1,36) /
     *   .5D0, 0.D0, 0.D0, 0.D0, 0.D0,  .5D0,
     *   0.D0, 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,
     *   0.D0, 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,
     *   .5D0, 0.D0, 0.D0, 0.D0, 0.D0,  .5D0 /

C     caseF
      data (prj(i,15),i=1,36) /
     *  .8D0,-.2D0, 0.D0,-.2D0,-.2D0,-.2D0,
     * -.2D0, .8D0, 0.D0,-.2D0,-.2D0,-.2D0,
     *  0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,
     * -.2D0,-.2D0, 0.D0, .8D0,-.2D0,-.2D0,
     * -.2D0,-.2D0, 0.D0,-.2D0, .8D0,-.2D0,
     * -.2D0,-.2D0, 0.D0,-.2D0,-.2D0, .8D0 /

      data (prjhat(i,15),i=1,36) /
     * .55D0, .05D0, 0.D0, .05D0,-.45D0,-.2D0,
     * .05D0, .55D0, 0.D0,-.45D0, .05D0,-.2D0,
     *  0.D0,  0.D0, 1.D0,  0.D0,  0.D0, 0.D0,
     * .05D0,-.45D0, 0.D0, .55D0, .05D0,-.2D0,
     *-.45D0, .05D0, 0.D0, .05D0, .55D0,-.2D0,
     * -.2D0, -.2D0, 0.D0, -.2D0, -.2D0, .8D0 /

      data (prjperp(i,15),i=1,36) /
     *  .2D0, .2D0, 0.D0, .2D0, .2D0, .2D0,
     *  .2D0, .2D0, 0.D0, .2D0, .2D0, .2D0,
     *  0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,
     *  .2D0, .2D0, 0.D0, .2D0, .2D0, .2D0,
     *  .2D0, .2D0, 0.D0, .2D0, .2D0, .2D0,
     *  .2D0, .2D0, 0.D0, .2D0, .2D0, .2D0 /


C     case678X g2=g4=g5=g6=0
      data (prj(i,16),i=1,36) /
     * 1.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 1.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0 /

      data (prjperp(i,16),i=1,36) /
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 1.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 1.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  1.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0, 0.D0,  0.D0,  1.D0 /


C     case9ABCDEX g1=g4=g5=g6=0
      data (prj(i,17),i=1,36) /
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0 /

      data (prjperp(i,17),i=1,36) /
     * 1.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  1.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  1.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  1.D0 /

C     caseFX g1=g2=g4=g5=g6=0
      data (prj(i,18),i=1,36) /
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0 /

      data (prjperp(1,18),i=1,36) /
     * 1.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  1.D0,  0.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  1.D0,  0.D0,
     * 0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  1.D0 /


C     DATA boundary67 /
      data (prj(i,19),i=1,36) /
     * 1.D0,  0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0,  0.D0,  0.D0,
     * 0.D0,  0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0,  .5D0, 0.D0,  .5D0,  0.D0,  0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0,  .5D0,  .5D0,
     * 0.D0,  0.D0, 0.D0,  0.D0,  .5D0,  .5D0 /

      data (prjperp(i,19),i=1,36) /
     * 0.D0,  0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0,  .5D0, 0.D0, -.5D0,  0.D0,  0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     * 0.D0, -.5D0, 0.D0,  .5D0,  0.D0,  0.D0,
     * 0.D0,  0.D0, 0.D0,  0.D0,  .5D0, -.5D0,
     * 0.D0,  0.D0, 0.D0,  0.D0, -.5D0,  .5D0 /

C     DATA boundary9A /
      data (prj(i,20),i=1,36) /
     *  .5D0, 0.D0, 0.D0,  0.D0,  .5D0,  0.D0,
     *  0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  0.D0,  .5D0,
     *  .5D0, 0.D0, 0.D0,  0.D0,  .5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  0.D0,  .5D0 /

      data (prjperp(i,20),i=1,36) /
     *  .5D0, 0.D0, 0.D0,  0.D0, -.5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  0.D0, -.5D0,
     * -.5D0, 0.D0, 0.D0,  0.D0, -.5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0, -.5D0,  0.D0,  .5D0 /


C     DATA boundaryCD /
      data (prj(i,21),i=1,36) /
     *  .5D0, 0.D0, 0.D0,  0.D0,  0.D0,  .5D0,
     *  0.D0, 1.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 1.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  .5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0,  .5D0,  0.D0,
     *  .5D0, 0.D0, 0.D0,  0.D0,  0.D0, .5D0 /

      data (prjperp(i,21),i=1,36) /
     *  .5D0, 0.D0, 0.D0,  0.D0,  0.D0, -.5D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  .5D0, -.5D0,  0.D0,
     *  0.D0, 0.D0, 0.D0, -.5D0,  .5D0,  0.D0,
     * -.5D0, 0.D0, 0.D0,  0.D0,  0.D0,  .5D0 /


C      case boundary12 g1=g2=g3
       data (prj(i,22),i=1,36) /
     *  0.3333333333D0, 0.3333333333D0,0.3333333333D0,
     *  0.D0,  0.D0,  0.D0,
     *  0.3333333333D0, 0.3333333333D0,0.3333333333D0,
     *  0.D0,  0.D0,  0.D0,
     *  0.3333333333D0, 0.3333333333D0,0.3333333333D0,
     *  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  1.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  1.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  1.D0 /


       data (prjperp(i,22),i=1,36) /
     *  0.6666666667D0,-0.3333333333D0,-0.3333333333D0,
     *  0.D0,  0.D0,  0.D0,
     * -0.3333333333D0, 0.6666666667D0,-0.3333333333D0,
     *  0.D0,  0.D0,  0.D0,
     * -0.3333333333D0,-0.3333333333D0, 0.6666666667D0,
     *  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0,
     *  0.D0, 0.D0, 0.D0,  0.D0,  0.D0,  0.D0 /

C
C    The following matrices are the transformation
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

C    M_1.M_2 (a -> b, b -> c, c -> a)
      DATA (MS(I,16), I=1,36) /
     2 0,0,1,0,0,0,
     3 1,0,0,0,0,0,
     4 0,1,0,0,0,0,
     5 0,0,0,0,0,1,
     6 0,0,0,1,0,0,
     7 0,0,0,0,1,0 /

C    M_2.M_1 (a -> c, b -> a, c -> b)
      DATA (MS(I,17), I=1,36) /
     2 0,1,0,0,0,0,
     3 0,0,1,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,1,0,
     6 0,0,0,0,0,1,
     7 0,0,0,1,0,0 /

C    M_2.M_1.M_2 (a -> c, c -> a)
      DATA (MS(I,18), I=1,36) /
     2 0,0,1,0,0,0,
     3 0,1,0,0,0,0,
     4 1,0,0,0,0,0,
     5 0,0,0,0,0,1,
     6 0,0,0,0,1,0,
     7 0,0,0,1,0,0/

      END

      Block Data histogram
      integer ihist(22,36)
      common /hist/ihist
      data ihist/792*0/
      end

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
      real*8 VT(6),VTRED(6)
      real*8 dmin, gerr, grmin
      real*8 XDOTVN, NCDIST
      integer IVB(NVMAX)
      integer i, nv1, nexamined, id
      integer nearst_N
      INTEGER ms(36,18)
      real*8 SEEDS(6,15)
      real*8 SEEDDIST(15)
      LOGICAL GOODSEED(15)
      LOGICAL NEARRED,NEAR2RED,NRED
      real*8 dred

      real*8 prj(36,22), prjperp(36,22)
      real*8 prjhat(36,15)
      common /projectors/ prj,prjperp,prjhat,ms





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

      do i = 1,22
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
      if (DMIN.GT.999999999.9989D0) THEN
        WRITE (*,'(a,G13.3)') ' DMIN IN MKGOAL ',DMIN
      else
        WRITE (*,'(a,F13.3)') ' DMIN IN MKGOAL ',DMIN
      endif
C     Populate the tree with the original cell and
C     the nearby neighbors on the boundary
      NV = NV+1
C      write(*,*) "calling bldtre_n"
      CALL BLDTRE_N (.FALSE., MXTREE,6, GRED,NV*16,TREE,'BLDTRE_N')
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
          IF (NEARST_N(.FALSE.,MXTREE,6,SEEDS(1,i),DMIN,TREE,NV1,ID,
     *      'NEARST_N') .EQ. 0) THEN
            NV = NV + 1
            CALL BLDTRE_N (.FALSE., MXTREE,6, SEEDS(1,i),
     *            NV*16+i,TREE,'BLDTRE_N')
            CALL CPYVN(6,SEEDS(1,i),V(1,NV))
            IVB(NV) = 1
            VDIST(NV) = SEEDDIST(i)
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, SEED DIST",NV,VDIST(NV)
C            call printg6('Store SEED ',SEEDS(1,i))
          ENDIF
          call imv6(SEEDS(1,i),MS(1,i),vt)
          NV1 = NV+1
          NRED=NEARRED(vt,1.0D-6,'NEARRED')
          if (NRED) THEN
          VDIST(NV1) = NCDIST(vt,GRED)
          if (VDIST(NV1).LT.gerr .and.
     *       NEARST_N(.false.,MXTREE,6,VT,DMIN,TREE,NV1,ID,
     *      'NEARST_N') .EQ. 0) THEN
            NV = NV+1
            CALL BLDTRE_N (.FALSE.,MXTREE,6,VT,NV*16+i,TREE,'BLDTRE_N')
            CALL CPYVN(6,VT,V(1,NV))
            IVB(NV) = 1
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, M*SEED DIST",NV,VDIST(NV)
C            call printg6('Store M*SEED ',VT)
          ENDIF
          ELSE
          NRED=NEAR2RED(vt,gerr,vtred,dred,'NEAR2RED')
          VDIST(NV1) = NCDIST(vtred,GRED)
          if (VDIST(NV1).LT.gerr .and.
     *       NEARST_N(.false.,MXTREE,6,vtred,DMIN,TREE,NV1,ID,
     *      'NEARST_N') .EQ. 0) THEN
            NV = NV+1
            CALL BLDTRE_N(.FALSE.,MXTREE,6,vtred,
     *          NV*16+i,TREE,'BLDTRE_N')
            CALL CPYVN(6,vtred,V(1,NV))
            IVB(NV) = 1
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, M*SEED ** DIST",NV,VDIST(NV)
C            call printg6('Store M*SEED ** ',VT)
          ENDIF
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

      gerr = DSQRT(XDOTVN(6,NEWGE,NEWGE))

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
          IF (NEARST_N(.false.,MXTREE,6,SEEDS(1,i),DMIN
     *       *(1.85D0**IVB(nexamined)),TREE,NV1,ID,
     *      'NEARST_N') .EQ. 0) THEN
            NV = NV + 1
            CALL BLDTRE_N (.FALSE., MXTREE,6, SEEDS(1,i),
     *           NV*16+i,TREE,'BLDTRE_N')
            CALL CPYVN(6,SEEDS(1,i),V(1,NV))
            IVB(NV) = IVB(nexamined)+1
            VDIST(NV) = dsqrt(SEEDDIST(i)**2+
     *        VDIST(nexamined)**2)
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, SEED DIST,IB",NV,VDIST(NV),IB
C            call printg6('Store SEED ',SEEDS(1,i))
          ENDIF
          call imv6(SEEDS(1,i),MS(1,i),vt)
          NV1 = NV+1
          NRED=NEARRED(vt,1.0D-6,'NEARRED')
          if (NRED) THEN
          VDIST(NV1) = NCDIST(vt,GRED)
          if (VDIST(NV1).LT.gerr .and.
     *       NEARST_N(.false.,MXTREE,6,VT,DMIN
     *       *(1.85D0**IVB(nexamined)),TREE,NV1,ID,
     *      'NEARST_N') .EQ. 0) THEN
            NV = NV+1
            CALL BLDTRE_N (.FALSE.,MXTREE,6,VT,
     *          NV*16+i,TREE,'BLDTRE_N')
            CALL CPYVN(6,VT,V(1,NV))
            IVB(NV) = IVB(nexamined)+1
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, M*SEED DIST,IB",NV,VDIST(NV),IB
C            call printg6('Store M*SEED ',VT)
          ENDIF
          ELSE
          NRED=NEAR2RED(vt,gerr,vtred,dred,'NEAR2RED')
          VDIST(NV1) = NCDIST(vtred,GRED)
          if (VDIST(NV1).LT.gerr .and.
     *       NEARST_N(.false.,MXTREE,6,vtred,DMIN
     *       *(1.85D0**IVB(nexamined)),TREE,NV1,ID,
     *      'NEARST_N') .EQ. 0) THEN
            NV = NV+1
            CALL BLDTRE_N (.FALSE.,MXTREE,6,vtred,
     *           NV*16+i,TREE,'BLDTRE_N')
            CALL CPYVN(6,vtred,V(1,NV))
            IVB(NV) = 1
            if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
C            write(*,*)"NV, M*SEED ** DIST",NV,VDIST(NV)
C            call printg6('Store M*SEED ** ',VT)
          ENDIF
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
      real*8 VECTOR(6)
C      real*8 FVECTOR(6)
      WRITE(*,'(A,1x,6F9.2)')TEXT,VECTOR

      RETURN
      END


      SUBROUTINE GENSEEDS( GVEC, GE, SEEDS, SEEDDIST,
     * GOODSEED, PRJ, PRJPERP)
C
C     GENSEEDS TAKES AN ARBITRARY G6 VECTOR AND PROJECTS
C     IT ONTO EACH OF THE 15 BOUNDARIES.  EVEN IF THE
C     ORIGINAL VECTOR WAS REDUCED THE PROJECTED VECTOR
C     MAY FAIL TO BE REDUCED, AND MAY EVEN FAIL TO BE
C     NEARLY REDUCED.  IT WILL BE MARKED AS A GOOD SEED
C     PROVIDED IT IS BOTH NEARLY REDUCED AND IS WITHIN
C     3.5 TIMES THE ERROR BOX GE OF GVEC.
C
      implicit none
      real*8 GVEC(6), GE(6)
      real*8 PRJ(36,22),PRJPERP(36,22)
      real*8 SEEDS(6,15), SEEDDIST(15)
      logical GOODSEED(15)
C      character*20 SEEDLAB
      real*8 vtemp(6)
      real*8 boundary67
      real*8 boundary9A
      real*8 boundaryCD
      real*8 a,dred
      integer i, ip
      real*8 anorm, gerr
      logical NEARRED,NEAR2RED
C----------------------------------------------------------------------C

C     Compute the seeds and the seed distances

C      write(*,*) "GENSEEDS GVEC, GE"
C      write(*,*) GVEC
C      write(*,*) GE

      gerr= anorm(6,GE)

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
        if (((GVEC(4)-1.d-6*sqrt(GVEC(2)*GVEC(3)))*
     *     (GVEC(5)-1.d-6*sqrt(GVEC(1)*GVEC(3)))*
     *     (GVEC(6)-1.d-6*sqrt(GVEC(1)*GVEC(3)))
     *     .le. 0.D0
     *        .and. (ip .eq. 6 .or. ip .eq. 7
     *          .or. ip .eq. 9 .or. ip .eq. 10
     *          .or. ip .eq. 12.or. ip .eq. 13))
     *        .or. (
     *     (GVEC(4)-1.d-6*sqrt(GVEC(2)*GVEC(3)))*
     *     (GVEC(5)-1.d-6*sqrt(GVEC(1)*GVEC(3)))*
     *     (GVEC(6)-1.d-6*sqrt(GVEC(1)*GVEC(3)))
     *     .gt. 0.D0
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
          GOODSEED(ip) = NEAR2RED(SEEDS(1,ip),gerr,
     *    vtemp,dred,'NEAR2RED')
          call CPYVN(6,vtemp,SEEDS(1,ip))
C         call printg6("Non-reduced projection",seeds(1,ip))

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


C     Compute the best distance between 2 G6 vectors
C     allowing for cell-preserving sign changes in
C     g4,5,6
C
      real*8 function g456dist(v1,v2)
      implicit none
      real*8 v1(6),v2(6),vtemp(6)
      real*8 xdot
      common/xdebug/xdebug
      logical xdebug

      integer ii

      xdot = 0.D0

      do ii = 1, 6
        vtemp(ii) = v1(ii)-v2(ii)
        xdot = xdot+vtemp(ii)*vtemp(ii)
      enddo
      g456dist = sqrt(xdot
     *  +4.D0*min(0.D0,
     *         v1(4)*v2(4)+v1(5)*v2(5),
     *         v1(4)*v2(4)+v1(6)*v2(6),
     *         v1(5)*v2(5)+v1(6)*v2(6)))
C      if (xdebug) then
C        write(7,'(a,1x,6f9.2)'),"g456dist v1",v1
C        write(7,'(a,1x,6f9.2)'),"g456dist v2",v2
C        write(7,'(a,1x,f9.2)'),"g456dist dist", g456dist
C      endif
      end

C
C     Compute the minimal distance between two Niggli-reduced
C     vectors in the Niggli Cone
C
C
      REAL*8 FUNCTION NCDIST(gv1,gv2)
      implicit none
      real*8 gv1(6), gv2(6)
      real*8 prj(36,22), prjperp(36,22)
      real*8 FOLDMDIST
      real*8 g456dist
      common/xdebug/xdebug
      logical xdebug
      integer ms(36,18)
      real*8 prjhat(36,15)
      common /projectors/ prj,prjperp,prjhat,ms

      if (xdebug) then
        write(7,'(a,1x,6f9.2)') "NCDIST gv1 ",gv1
        write(7,'(a,1x,6f9.2)') "NCDIST gv2 ",gv2
      endif
      NCDIST = g456dist(gv1,gv2)
      NCDIST=foldmdist(gv1,gv2,7,NCDIST)
      NCDIST=foldmdist(gv1,gv2,10,NCDIST)
      NCDIST=foldmdist(gv1,gv2,13,NCDIST)
      NCDIST=foldmdist(gv1,gv2,15,NCDIST)
      NCDIST=foldmdist(gv1,gv2,22,NCDIST)
      if (xdebug) then
        write(7,'(a,1x,6f9.2)') "NCDIST dist ",NCDIST
      endif

      end

C     DistF -- distance to the F boundary
      real*8 function DistF(g)
      implicit none
      real*8 g(6)

      DistF = abs(
     *  g(1)+g(2)+g(4)+g(5)+g(6)
     *  + 2.D0*min(0.D0,-g(4)-g(5),
     *  -g(4)-g(6),-g(5)-g(6)))/sqrt(5.D0)
      end


C     FOLDMDIST(GVEC1,GVEC2,IP,CFOLDDIST)
C       return the distance betweem two F6 vectors
C       going through the given boundary
C       (7, A, C, F) after taking the second
C       vector through the 6 exchanges of A, B, C
C
      real*8 function FOLDMDIST(GVEC1,GVEC2,IP,CFOLDDIST)
      implicit none

      REAL*8 GVEC1(6), GVEC2(6), CFOLDDIST
      REAL*8 FOLDxDIST
      integer ip
      REAL*8 m1g1(6),m2g1(6),m1m2g1(6)
      REAL*8 m2m1g1(6),m2m1m2g1(6)
      REAL*8 m1g2(6),m2g2(6),m1m2g2(6)
      REAL*8 m2m1g2(6),m2m1m2g2(6)
      real*8 prj(36,22), prjperp(36,22)
      integer ms(36,18)
      real*8 prjhat(36,15)
      common /projectors/ prj,prjperp,prjhat,ms
      common/xdebug/xdebug
      logical xdebug
      real*8 h1,h1m1,h1m2,h1m1m2,h1m2m1,h1m2m1m2
      real*8 h2,h2m1,h2m2,h2m1m2,h2m2m1,h2m2m1m2
      real*8 distF
    
      call imv6(gvec1,MS(1,1),m1g1)
      call imv6(gvec1,MS(1,2),m2g1)
      call imv6(gvec1,MS(1,16),m1m2g1)
      call imv6(gvec1,MS(1,17),m2m1g1)
      call imv6(gvec1,MS(1,18),m2m1m2g1)
      call imv6(gvec2,MS(1,1),m1g2)
      call imv6(gvec2,MS(1,2),m2g2)
      call imv6(gvec2,MS(1,16),m1m2g2)
      call imv6(gvec2,MS(1,17),m2m1g2)
      call imv6(gvec2,MS(1,18),m2m1m2g2)
      FOLDMDIST = CFOLDDIST
      
 
      if (ip.ge. 6 .and. ip.le. 8) go to 1007
      if (ip.ge. 9 .and. ip.le.11) go to 1010
      if (ip.ge.12 .and. ip.le.14) go to 1013
      if (ip.eq.15) go to 1015
      if (ip.eq.22) go to 1022

 1007 continue
      h1 = abs(gvec1(2)-abs(gvec1(4)))/sqrt(2.D0)
      h2 = abs(gvec2(2)-abs(gvec2(4)))/sqrt(2.D0)
      h1m1 = abs(m1g1(2)-abs(m1g1(4)))/sqrt(2.D0)
      h2m1 = abs(m1g2(2)-abs(m1g2(4)))/sqrt(2.D0)
      h1m2 = abs(m2g1(2)-abs(m2g1(4)))/sqrt(2.D0)
      h2m2 = abs(m2g2(2)-abs(m2g2(4)))/sqrt(2.D0)
      h1m1m2 = abs(m1m2g1(2)-abs(m1m2g1(4)))/sqrt(2.D0)
      h2m1m2 = abs(m1m2g2(2)-abs(m1m2g2(4)))/sqrt(2.D0)
      h1m2m1 = abs(m2m1g1(2)-abs(m2m1g1(4)))/sqrt(2.D0)
      h2m2m1 = abs(m2m1g2(2)-abs(m2m1g2(4)))/sqrt(2.D0)
      h1m2m1m2 = 
     * abs(m2m1m2g1(2)-abs(m2m1m2g1(4)))/sqrt(2.D0)
      h2m2m1m2 =
     * abs(m2m1m2g2(2)-abs(m2m1m2g2(4)))/sqrt(2.D0)
      go to 1100

 1010 continue
      h1 = abs(gvec1(1)-abs(gvec1(5)))/sqrt(2.D0)
      h2 = abs(gvec2(1)-abs(gvec2(5)))/sqrt(2.D0)
      h1m1 = abs(m1g1(1)-abs(m1g1(5)))/sqrt(2.D0)
      h2m1 = abs(m1g2(1)-abs(m1g2(5)))/sqrt(2.D0)
      h1m2 = abs(m2g1(1)-abs(m2g1(5)))/sqrt(2.D0)
      h2m2 = abs(m2g2(1)-abs(m2g2(5)))/sqrt(2.D0)
      h1m1m2 = abs(m1m2g1(1)-abs(m1m2g1(5)))/sqrt(2.D0)
      h2m1m2 = abs(m1m2g2(1)-abs(m1m2g2(5)))/sqrt(2.D0)
      h1m2m1 = abs(m2m1g1(1)-abs(m2m1g1(5)))/sqrt(2.D0)
      h2m2m1 = abs(m2m1g2(1)-abs(m2m1g2(5)))/sqrt(2.D0)
      h1m2m1m2 =
     * abs(m2m1m2g1(1)-abs(m2m1m2g1(5)))/sqrt(2.D0)
      h2m2m1m2 =
     * abs(m2m1m2g2(1)-abs(m2m1m2g2(5)))/sqrt(2.D0)
      go to 1100

 1013 continue
      h1 = abs(gvec1(1)-abs(gvec1(6)))/sqrt(2.D0)
      h2 = abs(gvec2(1)-abs(gvec2(6)))/sqrt(2.D0)
      h1m1 = abs(m1g1(1)-abs(m1g1(6)))/sqrt(2.D0)
      h2m1 = abs(m1g2(1)-abs(m1g2(6)))/sqrt(2.D0)
      h1m2 = abs(m2g1(1)-abs(m2g1(6)))/sqrt(2.D0)
      h2m2 = abs(m2g2(1)-abs(m2g2(6)))/sqrt(2.D0)
      h1m1m2 = abs(m1m2g1(1)-abs(m1m2g1(6)))/sqrt(2.D0)
      h2m1m2 = abs(m1m2g2(1)-abs(m1m2g2(6)))/sqrt(2.D0)
      h1m2m1 = abs(m2m1g1(1)-abs(m2m1g1(6)))/sqrt(2.D0)
      h2m2m1 = abs(m2m1g2(1)-abs(m2m1g2(6)))/sqrt(2.D0)
      h1m2m1m2 =
     * abs(m2m1m2g1(1)-abs(m2m1m2g1(6)))/sqrt(2.D0)
      h2m2m1m2 =
     * abs(m2m1m2g2(1)-abs(m2m1m2g2(6)))/sqrt(2.D0)
      go to 1100


 1015 continue
      h1 = DistF(gvec1)
      h2 = DistF(gvec2)
      h1m1 = DistF(m1g1)
      h2m1 = DistF(m1g2)
      h1m2 = DistF(m2g1)
      h2m2 = DistF(m2g2)
      h1m1m2 = DistF(m1m2g1)
      h2m1m2 = DistF(m1m2g2)
      h1m2m1 = DistF(m2m1g1)
      h2m2m1 = DistF(m2m1g2)
      h1m2m1m2 = DistF(m2m1m2g1)
      h2m2m1m2 = DistF(m2m1m2g2)
      goto 1100

 1022 continue
      h1 = sqrt(2.D0*(gvec1(1)**2+gvec1(2)**2+gvec1(3)**2
     *     - gvec1(2)*gvec1(3)
     *     - gvec1(1)*gvec1(3)
     *     - gvec1(1)*gvec1(3))/3.D0)
      h2 = sqrt(2.D0*(gvec2(1)**2+gvec2(2)**2+gvec2(3)**2
     *     - gvec2(2)*gvec2(3)
     *     - gvec2(1)*gvec2(3)
     *     - gvec2(1)*gvec2(3))/3.D0)
      h1m1 = h1
      h2m1 = h2
      h1m2 = h1
      h2m2 = h2
      h1m1m2 = h1
      h2m1m2 = h2
      h1m2m1 = h1
      h2m2m1 = h2
      h1m2m1m2 = h1
      h2m2m1m2 = h2


 1100 continue
      if (h1.lt.FOLDMDIST) then
      if (h1+h2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(GVEC1,GVEC2,ip,FOLDMDIST,1,1)
      if (h1+h2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(GVEC1,m1g2,ip,FOLDMDIST,1,2)
      if (h1+h2m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(GVEC1,m2g2,ip,FOLDMDIST,1,3)
      if (h1+h2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(GVEC1,m1m2g2,ip,FOLDMDIST,1,4)
      if (h1+h2m2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(GVEC1,m2m1g2,ip,FOLDMDIST,1,5)
      if (h1+h2m2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(GVEC1,m2m1m2g2,ip,FOLDMDIST,1,6)
      endif

      if (h1m1.lt.FOLDMDIST) then
      if (h1m1+h2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1g1,GVEC2,ip,FOLDMDIST,2,1)
      if (h1m1+h2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1g1,m1g2,ip,FOLDMDIST,2,2)
      if (h1m1+h2m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1g1,m2g2,ip,FOLDMDIST,2,3)
      if (h1m1+h2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1g1,m1m2g2,ip,FOLDMDIST,2,4)
      if (h1m1+h2m2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1g1,m2m1g2,ip,FOLDMDIST,2,5)
      if (h1m1+h2m2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1g1,m2m1m2g2,ip,FOLDMDIST,2,6)
      endif

      if (h1m2.lt.FOLDMDIST) then
      if (h1m2+h2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2g1,GVEC2,ip,FOLDMDIST,3,1)
      if (h1m2+h2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2g1,m1g2,ip,FOLDMDIST,3,2)
      if (h1m2+h2m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2g1,m2g2,ip,FOLDMDIST,3,3)
      if (h1m2+h2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2g1,m1m2g2,ip,FOLDMDIST,3,4)
      if (h1m2+h2m2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2g1,m2m1g2,ip,FOLDMDIST,3,5)
      if (h1m2+h2m2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2g1,m2m1m2g2,ip,FOLDMDIST,3,6)
      endif

      if (h1m1m2.lt.FOLDMDIST) then
      if (h1m1m2+h2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1m2g1,GVEC2,ip,FOLDMDIST,4,1)
      if (h1m1m2+h2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1m2g1,m1g2,ip,FOLDMDIST,4,2)
      if (h1m1m2+h2m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1m2g1,m2g2,ip,FOLDMDIST,4,3)
      if (h1m1m2+h2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1m2g1,m1m2g2,ip,FOLDMDIST,4,4)
      if (h1m1m2+h2m2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1m2g1,m2m1g2,ip,FOLDMDIST,4,5)
      if (h1m1m2+h2m2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m1m2g1,m2m1m2g2,ip,FOLDMDIST,4,6)
      endif

      if (h1m2m1.lt.FOLDMDIST) then
      if (h1m2m1+h2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1g1,GVEC2,ip,FOLDMDIST,5,1)
      if (h1m2m1+h2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1g1,m1g2,ip,FOLDMDIST,5,2)
      if (h1m2m1+h2m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1g1,m2g2,ip,FOLDMDIST,5,3)
      if (h1m2m1+h2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1g1,m1m2g2,ip,FOLDMDIST,5,4)
      if (h1m2m1+h2m2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1g1,m2m1g2,ip,FOLDMDIST,5,5)
      if (h1m2m1+h2m2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1g1,m2m1m2g2,ip,FOLDMDIST,5,6)
      endif

      if (h1m2m1m2.lt.FOLDMDIST) then
      if (h1m2m1m2+h2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1m2g1,GVEC2,ip,FOLDMDIST,6,1)
      if (h1m2m1m2+h2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1m2g1,m1g2,ip,FOLDMDIST,6,2)
      if (h1m2m1m2+h2m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1m2g1,m2g2,ip,FOLDMDIST,6,3)
      if (h1m2m1m2+h2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1m2g1,m1m2g2,ip,FOLDMDIST,6,4)
      if (h1m2m1m2+h2m2m1.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1m2g1,m2m1g2,ip,FOLDMDIST,6,5)
      if (h1m2m1m2+h2m2m1m2.lt.FOLDMDIST)
     * FOLDMDIST = FOLDxDIST(m2m1m2g1,m2m1m2g2,ip,FOLDMDIST,6,6)
      endif

      end


C     FOLDNDIST -- return the distance between two
C     G6 vectors going through the boundary N
C
      REAL*8 FUNCTION FOLDNDIST(GVEC1,GVEC2,N,CFOLDDIST)

      implicit none

      REAL*8 GVEC1(6), GVEC2(6), CFOLDDIST
      real*8 g1(6,4), g2(6,4)
      REAL*8 pg1(6),pg2(6),mpg1(6),mpg2(6)
      REAL*8 ppg1(6),ppg2(6)
      real*8 prj(36,22), prjperp(36,22)
      integer ms(36,18)
      real*8 prjhat(36,15)
      common /projectors/ prj,prjperp,prjhat,ms
      real*8 ht1, ht2, g456dist, dist12, anorm
      integer N, i, j

      FOLDNDIST = CFOLDDIST
      
      call cpyvn(6,gvec1,g1(1,1))
      call cpyvn(6,gvec2,g2(1,1))
      if (N.eq.22) then
        call rmv6(g1(1,1),prj(1,N),pg1)
        call rmv6(g1(1,1),prjperp(1,N),ppg1)
        call rmv6(g2(1,1),prj(1,N),pg2)
        call rmv6(g2(1,1),prjperp(1,N),ppg2)
        ht1 = anorm(6,ppg1)
        ht2 = anorm(6,ppg2)
        if (ht1+ht2.lt.FOLDNDIST)
     *    FOLDNDIST= min(FOLDNDIST,
     *      sqrt(g456dist(pg1,pg2)**2
     *      +(ht1+ht2)**2))
        return
      endif
      do i = 3,5
        call imv6(gvec1,ms(1,i),g1(1,i-1))
        call imv6(gvec2,ms(1,i),g2(1,i-1))
      enddo
      do i = 1,4
        call rmv6(g1(1,i),prj(1,N),pg1)
        call rmv6(g1(1,i),prjperp(1,N),ppg1)
        call imv6(pg1,ms(1,N),mpg1)
        ht1 = anorm(6,ppg1)
        if (ht1.lt.FOLDNDIST) then
        do j = 1,4
          call rmv6(g2(1,j),prj(1,N),pg2)
          call rmv6(g2(1,i),prjperp(1,N),ppg2)
          call imv6(pg2,ms(1,N),mpg2)
          ht2 = anorm(6,ppg2)
          if (ht1+ht2.lt.FOLDNDIST) then
          dist12 = min(g456dist(pg1,pg2),
     *      g456dist(mpg1,pg2),
     *      g456dist(pg1,mpg2),
     *      g456dist(mpg1,mpg2))
          FOLDNDIST= min(FOLDNDIST,
     *      sqrt(dist12**2+(ht1+ht2)**2))
          endif
        enddo
        endif
      enddo
      end

C     FOLDXDIST -- return the distance between two
C     G6 vectors going through the given boundary
C     (7, A, C or F)
C
C     is1 and is2 are the bdy 1,2 symmetry op used
C
      REAL*8 FUNCTION FOLDXDIST(GVEC1,GVEC2,IP,CFOLDDIST,
     *  is1, is2)

      implicit none

      REAL*8 GVEC1(6), GVEC2(6), CFOLDDIST
      REAL*8 FOLDNDIST, ofolddist
      integer ip, is1, is2
      REAL*8 pg1(6),pg2(6),mpg1(6),mpg2(6)
      real*8 prj(36,22), prjperp(36,22)
      integer ms(36,18)
      real*8 prjhat(36,15)
      common /projectors/ prj,prjperp,prjhat,ms
      integer i
      real*8 ht1, ht2, g456dist, dist12
      common/xdebug/xdebug
      logical xdebug
      integer ihist(22,36)
      common /hist/ihist


      if (xdebug) then
        write(7,*) "FOLDXDIST ip, CFOLDDIST: ",ip,CFOLDDIST
        write(7,'(a,1x,6f9.2)') "FOLDXDIST gvec1: ",gvec1
        write(7,'(a,1x,6f9.2)') "FOLDXDIST gvec2: ",gvec2
      endif

      ofolddist = cfolddist

      if (ip.eq.15.or.ip.eq.22) then
        foldxdist = 
     *    foldndist(gvec1,gvec2,15,cfolddist)
        if (foldxdist.lt.ofolddist)
     *    ihist(ip,is1+(is2-1)*6)
     *    =ihist(ip,is1+(is2-1)*6)+1
        return
      endif

      do i = 1,6
        pg1(i) = ABS(GVEC1(i))
        pg2(i) = ABS(GVEC2(i))
        mpg1(i) = ABS(GVEC1(i))
        mpg2(i) = ABS(GVEC2(i))
      enddo

      if (ip.eq.7.or.ip.eq.6.or.ip.eq.8) then
        pg1(6) = GVEC1(6)
        pg2(6) = GVEC2(6)
        pg1(2) = (pg1(2)+pg1(4))/2.D0
        pg1(4) = pg1(2)
        pg2(2) = (pg2(2)+pg2(4))/2.D0
        pg2(4) = pg2(2)
        mpg1(2) = pg1(2)
        mpg1(4) = pg1(2)
        mpg2(2) = pg2(2)
        mpg2(4) = pg2(2)
        mpg1(6) = pg1(6)
        mpg2(6) = pg2(6)
        mpg1(5) = pg1(6)-pg1(5)
        mpg2(5) = pg2(6)-pg2(5)
        ht1 = abs(GVEC1(2)-ABS(GVEC1(4)))/SQRT(2.D0)
        ht2 = abs(GVEC2(2)-ABS(GVEC2(4)))/SQRT(2.D0)
      else if (ip.eq.10.or.ip.eq.9.or.ip.eq.11) then
        pg1(6) = GVEC1(6)
        pg2(6) = GVEC2(6)
        pg1(1) = (pg1(1)+pg1(5))/2.D0
        pg1(5) = pg1(1)
        pg2(1) = (pg2(1)+pg2(5))/2.D0
        pg2(5) = pg2(1)
        mpg1(1) = pg1(1)
        mpg1(5) = pg1(1)
        mpg2(1) = pg2(1)
        mpg2(5) = pg2(1)
        mpg1(6) = pg1(6)
        mpg2(6) = pg2(6)
        mpg1(4) = pg1(6)-pg1(4)
        mpg2(4) = pg2(6)-pg2(4)
        ht1 = abs(GVEC1(1)-ABS(GVEC1(5)))/SQRT(2.D0)
        ht2 = abs(GVEC2(1)-ABS(GVEC2(5)))/SQRT(2.D0)
      else if (ip.eq.13.or.ip.eq.12.or.ip.eq.14) then
        pg1(5) = GVEC1(5)
        pg2(5) = GVEC2(5)
        pg1(1) = (pg1(1)+pg1(6))/2.D0
        pg1(6) = pg1(1)
        pg2(1) = (pg2(1)+pg2(6))/2.D0
        pg2(6) = pg2(1)
        mpg1(1) = pg1(1)
        mpg1(6) = pg1(1)
        mpg2(1) = pg2(1)
        mpg2(6) = pg2(1)
        mpg1(5) = pg1(5)
        mpg2(5) = pg2(5)
        mpg1(4) = pg1(5)-pg1(4)
        mpg2(4) = pg2(5)-pg2(4)
        ht1 = abs(GVEC1(1)-ABS(GVEC1(6)))/SQRT(2.D0)
        ht2 = abs(GVEC2(1)-ABS(GVEC2(6)))/SQRT(2.D0)
      endif

      if (ht1+ht2 .ge. CFOLDDIST) then
        FOLDXDIST = CFOLDDIST
      else
        dist12 = min(g456dist(pg1,pg2),
     *    g456dist(pg1,mpg2),
     *    g456dist(mpg1,pg2),
     *    g456dist(mpg1,mpg2))
        FOLDXDIST = MIN(CFOLDDIST,SQRT((ht1+ht2)**2 +
     *    dist12**2))
        if (foldxdist.lt.ofolddist)
     *    ihist(ip,is1+(is2-1)*6)
     *    =ihist(ip,is1+(is2-1)*6)+1
      endif
      if (xdebug) then
          write(7,'(a,1x,6f9.2)') "FOLDXDIST pg1: ",pg1
          write(7,'(a,1x,6f9.2)') "FOLDXDIST mpg1: ",mpg1
          write(7,'(a,1x,6f9.2)') "FOLDXDIST pg2: ",pg2
          write(7,'(a,1x,6f9.2)') "FOLDXDIST mpg2: ",mpg2
          write(7,'(a,1x,6f9.2)') "FOLDXDIST: ",FOLDXDIST
      endif
        RETURN
      end

