#
#  Makefile for ITERATE
#
#  Herbert J. Bernstein, Bernstein + Sons
#  Lawrence C. Andrews, Thuridion, Inc.
#
#  29 Sep 1996
#
#
#  Modify the following definitions for your system
#
#  HTTPDSERVER is the name of the server on which the
#  installation is being made
#
#  *************************************************
#  *** YOU MUST CHANGE THIS DEFINITION TO PERMIT ***
#  ***              REMOTE ACCESS                ***
#  *************************************************
#
HTTPDSERVER	=	localhost
#  The following are normal defaults for a system manager
#  installation assuming an NCSA httpd default installation
#
#  BINDEST is the installation directory for the executable
#  of ITERATE
BINDEST		=	/usr/local/bin
#
#  CGIBIN is the installation directory for the cgi-bin script
#  iterate.csh
CGIBIN		=	/usr/local/etc/httpd/cgi-bin
#
#  CGIBINEXT is the external name of the directory for the
#  cgi-bin script iterate.csh
CGIBINEXT	=	/cgi-bin
#
#  HTDOCS is the installation directory for the HTML document
#  iterate.html
HTDOCS		=	/usr/local/etc/httpd/htdocs
#
#  For a user installation you need the system manager to have
#  permitted cgi-bin execution from the directory given
#  The following lines, with the "??????" replaced by a valid
#  user name are a possible start on user installation definitions
#
#USERNAME	=	??????
#BINDEST	=	/home/$(USERNAME)/bin
#CGIBIN		=	/home/$(USERNAME)/public_html/cgi-bin
#CGIBINEXT	=	/~$(USERNAME)/cgi-bin
#HTDOCS		=	/home/$(USERNAME)/public_html
#
#  Default compile flag definition to select debug mode under unix
FFLAGS	=	-g
#
#  For IBM AIX xlf compilation with full optimization try this
#FFLAGS	=	-O3 -qstrict
#FC	=	xlf
#
HTFLAGS 	=	-DFULLHTDOCS=$(CGIPATH)
#
#  For use of wwwcount2.3
#HTFLAGS	=	-DFULLHTDOCS=$(CGIPATH) -DWWWCOUNT=TRUE
#
#
#  You should not have to edit below this line
#********************************************************************
#
#
CGIPATH	=	http://$(HTTPDSERVER)$(CGIBINEXT)/iterate.csh
BINPATH	=	$(BINDEST)/ITERATE
#
all:		edit 
#
edit:	
		@/bin/echo "**************************************"
		@/bin/echo "* You must edit Makefile before      *"
		@/bin/echo "* installing ITERATE                 *"
		@/bin/echo "* Then:                              *"
		@/bin/echo "*     make edit_done                 *"
		@/bin/echo "*     make all                       *"
		@/bin/echo "**************************************"
#
edit_done:	ITERATE iterate.html iterate.csh
		touch edit
#
clean:
		-rm edit
		-rm iterate.html
		-rm ITERATE
		-rm iterate.csh
		-rm *.bak

#
iterate.html:	iterate.html.m4 Makefile
		m4 $(HTFLAGS) < iterate.html.m4 > iterate.html
#
iterate.csh:	iterate.csh.m4 Makefile
		m4 -DBINPATH=$(BINPATH) < iterate.csh.m4 > iterate.csh
#
install:	ITERATE iterate.csh iterate.html iterate.cshar iterate.shar
		-mkdir -p $(BINDEST)
		-mkdir -p $(CGIBIN)
		-mkdir -p $(HTDOCS)
		chmod 755 ITERATE
		chmod 755 iterate.csh
		cp ITERATE $(BINDEST)
		cp iterate.csh $(CGIBIN)
		cp iterate.html $(HTDOCS)
		cp iterate.cshar $(HTDOCS)
		cp iterate.shar $(HTDOCS)
#		
ciftbx.o:	ciftbx.f ciftbx.sys ciftbx.cmn ciftbx.cmf ciftbx.cmv clearfp.f
hash_funcs.o:	hash_funcs.f
iterate.o:	iterate.f ITERATE.cmn ciftbx.cmn ciftbx.cmf ciftbx.cmv
ITERATE:	ITERATE.cmn iterate.o ciftbx.o hash_funcs.o 
	$(FC) $(FFLAGS) -o ITERATE iterate.o ciftbx.o hash_funcs.o

iterate.shar:	MANIFEST README Makefile Makefile_save iterate.f ITERATE.cmn \
	ciftbx.f ciftbx.cmn ciftbx.sys ciftbx.cmf ciftbx.cmv hash_funcs.f \
	iterate.html.m4 iterate.csh.m4 cryst1-2-cif.awk clearfp.f \
	clearfp_sun.f
	-rm iterate.shar
	makekit -s5000k -m
	mv Part01 iterate.shar

iterate.cshar:	MANIFEST README Makefile Makefile_save iterate.f ITERATE.cmn \
	ciftbx.f ciftbx.cmn ciftbx.sys ciftbx.cmf ciftbx.cmv hash_funcs.f \
	iterate.html.m4 iterate.csh.m4 cryst1-2-cif.awk clearfp.f \
	clearfp_sun.f
	-rm iterate.cshar
	makekit -c -s5000k -m
	mv Part01 iterate.cshar
