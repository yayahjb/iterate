C
C
C    \ | /            /##|    @@@@  @   @@@@@   |      |              @@@
C     \|/ STAR       /###|   @      @   @     __|__    |             @   @
C  ----*----        /####|  @       @   @@@@    |      |___  __  __     @
C     /|\          /#####|   @      @   @       |      |   \   \/      @
C    / | \         |#####|    @@@@  @   @       \___/  \___/ __/\__  @@@@@
C                  |#####|________________________________________________
C                 ||#####|                 ___________________            |
C        __/|_____||#####|________________|&&&&&&&&&&&&&&&&&&&||          |
C<\\\\\\\\_ |_____________________________|&&&& Sep 26 96 &&&&||          |
C          \|     ||#####|________________|&&&&&&&&&&&&&&&&&&&||__________|
C                  |#####|
C                  |#####|                Version 2.5.1 Release
C                  |#####|
C                 /#######\ 
C                |#########|
C                    ====
C                     ||
C           An extended tool box of fortran routines for manipulating CIF data.
C                     ||
C                     ||  CIFtbx Version 2
C                     ||        by
C                     ||
C                     ||  Sydney R. Hall (syd@crystal.uwa.edu.au)
C                     ||  Crystallography Centre
C                     ||  University of Western Australia
C                     ||  Nedlands 6009, AUSTRALIA
C                     ||
C                     ||       and
C                     ||
C                     ||  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
C                     ||  Bernstein + Sons
C                     ||  5 Brewster Lane
C                     ||  Bellport, NY 11713, U.S.A.
C                     ||
C The latest program source and information is available from:
C                     ||
C Em: syd@crystal.uwa.edu.au       ,-_|\      Sydney R. Hall
C sendcif@crystal.uwa.edu.au      /     \     Crystallography Centre
C Fx: +61 9 380 1118  ||      --> *_,-._/     University of Western Australia
C Ph: +61 9 380 2725  ||               v      Nedlands 6009, AUSTRALIA
C                     ||
C                     ||
C_____________________||_____________________________________________________
C
C This is a version of CIFtbx which has been extended to work with DDL 2
C and mmCIF as well as with DDL 1.4 and core CIF dictionaries.  CIFtbx
C version 1 was written by Sydney R. Hall (see Hall, S. R., "CIF Applications
C IV.  CIFtbx: a Tool Box for Manipulating CIFs,"  J. Appl. Cryst (1993). 26,
C 482-494.  The revisions for version 2 were done by Herbert J. Bernstein
C and Sydney R. Hall (see Hall, S. R. and Bernstein, H. J., "CIFtbx 2:
C Extended Tool Box for Manipulating CIFs," J. Appl. Cryst., to appear.) 
C
C___________________________________________________________________________
C
C
C    GENERAL TOOLS
C
C
C    init_      Sets the device numbers of files.   (optional)
C               [logical function always returned .true.]
C
C               <input CIF dev number> Set input CIF device     (def=1)
C
C               <output CIF dev number>Set output CIF device    (def=2)
C
C               <diracc dev number>    Set direct access formatted
C                                      scratch device number    (def=3)
C
C               <error  dev number>    Set error message device (def=6)
C
C
C
C    dict_      Requests a CIF dictionary be used for various data checks.
C               [logical function returned as .true. if the name dictionary
C               was opened; and if the check codes are recognisable.  The
C               data item names used in the first dictionary loaded are
C               considered to be preferred by the user to aliases found
C               in dictionaries loaded in later calls]
C
C               <dictionary filename>  A CIF dictionary in DDL format
C                                      or blank if just setting flags
C                                      or resetting the dictionary
C
C               <check code string>    The codes specifying the types of 
C                                      checks to be applied to the CIF.
C
C                                      'valid'  data name validation check.
C                                      'dtype'  data item data type check.
C                                      'reset'  switch off checking flags
C                                      'close'  close existing dictionaries
C
C___________________________________________________________________________
C
C
C   CIF ACCESS TOOLS  ("the get_ing commands")
C
C
C
C    ocif_      Opens the CIF containing the required data.
C               [logical function returned .true. if CIF opened]
C
C               <CIF filename>        A blank name signals that the
C                                     currently open input CIF file
C                                     will be read.
C
C
C
C    data_      Identifies the data block containing the data to be requested. 
C               [logical function returned .true. if block found]
C
C               <data block name>     A blank name signals that the next
C                                     encountered block is used (the block
C                                     name is stored in the variable bloc_).
C
C
C    bkmrk_     Saves or restores the current position so that data from 
C               elsewhere in the cif can be examined.
C               [logical function returned as .true. on save if there was
C               room in internal storage to hold the current position, .true.
C               on restore if the bookmark number used was valid.  If the
C               argument is zero, the call is to save the position and return
C               the bookmark number in the argument.  If the argument is
C               non-zero, the call is to restore the position saved for the
C               bookmark number given.  The bookmark and the argument are
C               cleared.  The position set on return allow reprocessing of
C               the data item or loop row last processed when the bookmark
C               was placed.
C
C               NOTE:  All bookmarks are cleared by a call to data_]
C
C               <integer variable>    Bookmark number
C
C
C    find_      Find the location of the requested item in the CIF.
C               [The argument "name" may be a data item name, blank
C               for the next such item.  The argument "type" may be
C               blank for unrestricted acceptance of any non-comment
C               string (use cmnt_ to see comments), including loop headers,
C               "name" to accept only the name itself and "valu"
C               to accept only the value, or "head" to position to the
C               head of the CIF.  Except when the "head" is requested,
C               the position is left after the data item provided.  If the
C               item found is of type "name", posnam_ is set, otherwise, 
C               posval_]
C
C               <data item name>      A blank name signals that the next
C                                     item of the type specified is needed
C
C               <data item type>      blank, 'head', 'name' or 'valu'
C
C               <character variable>  Returned string is of length long_.
C
C
C
C    test_      Identify the data attributes of the named data item.
C               [logical function returned as .true. if the item is present or
C               .false. if it is not. The data attributes are stored in the
C               common variables list_, type_, dictype_, diccat_ and dicname_. 
C               The values in dictype_, diccat_ and dicname_ are valid
C               whether or not the data item is found in the input CIF, as
C               long as the named data item is found in the dictionaries
C               declared by calls to dict_.  The data item name found
C               in the input CIF is stored in tagname_.  The appropriate
C               column numbers are stored in posnam_, posval_, posend_ and (for
C               numbers) in posdec_.  The quoation mark, if any, used is
C               stored in quote_.
C
C               list_ is an integer variable containing the sequential number
C               of the loop block in the data block. If the item is not within
C               a loop structure this value will be zero.
C
C               type_ is a character*4 variable with the possible values:
C                      'numb'  for number data
C                      'char'  for character data
C                      'text'  for text data
C                      'null'  if data missing or '?' or '.'
C
C               dictype_ is a character*(NUMCHAR) variable with the type code
C               given in the dictionary entry for the named data item.  If
C               no dictionary was used, or no type code was specified, this
C               field will simply agree with type_.  If a dictionary was used,
C               this type may be more specific than the one given by type_.
C
C               diccat_ is a character*(NUMCHAR) variable with the category
C               of the named data item, or '(none)'
C
C               dicname_ is a character*(NUMCHAR) variable with the name of
C               the data item which is found in the dictionary for the
C               named data item.  If alias_ is .true., this name may
C               differ from the name given in the call to test_.  If alias_
C               is .false. or no preferred alias is found, dicname_ agrees with
C               the data item name.
C
C               tagname_ is a character*(NUMCHAR) variable with the name
C               of the data item as found in the input CIF.  It will be
C               blank if the data item name requested is not found in the
C               input CIF and may differ from the data item name provided
C               by the user if the name used in the input CIF is an
C               alias of the data item name and alias_ is .true.
C
C               posnam_, posval_, posend_  and posdec_ are integer variables
C               which may be examined if information about the horizontal
C               position of the name and data read are needed.  posnam_ is the
C               starting column of the data name found (most often 1).
C               posval_ is the starting column of the data value.  If the
C               field is numeric, then posdec_ will contain the effective
C               column number of the decimal point.  For whole numbers, the
C               effective position of the decimal point is one column to the
C               right of the field.  posend_ contains the ending column of the
C               data value.
C
C               quote_ is a character*1 varibale which may be examined to
C               determine if a quotation character was used on character data.]
C
C               <data name>           Name of the data item to be tested.
C
C
C
C    name_      Get the NEXT data name in the current data block.
C               [logical function returned as .true. if a new data name exists
C               in the current data block, and .false. when the end of the data
C               block is reached.]
C
C               <data name>           Returned name of next data item in block.
C
C
C
C    numb_      Extracts the number and its standard deviation (if appended).
C               [logical function returned as .true. if number present. If
C               .false. arguments 2 and 3 are unaltered. If the esd is not
C               attached to the number argument 3 is unaltered.]
C
C               <data name>           Name of the number sought.
C
C               <real variable>       Returned number.
C
C               <real variable>       Returned standard deviation.
C
C
C
C    numd_      Extracts the number and its standard deviation (if appended)
C               as double precision variables.
C               [logical function returned as .true. if number present. If
C               .false. arguments 2 and 3 are unaltered. If the esd is not
C               attached to the number argument 3 is unaltered.]
C
C               <data name>           Name of the number sought.
C
C               <double precision variable>
C                                     Returned number.
C
C               <double precision variable>
C                                     Returned standard deviation.
C
C
C
C    char_      Extracts character and text strings.
C               [logical function returned as .true. if the string is present.
C               Note that if the character string is text this function is 
C               called repeatedly until the logical variable text_ is .false.]
C
C               <data name>           Name of the string sought.
C
C               <character variable>  Returned string is of length long_.
C
C
C    cmnt_      Extracts the next comment from the data block.
C               [logical function returned as .true. if a comment is present.
C               The initial comment character "#" is _not_ included in the
C               returned string.  A completely blank line is treated as
C               a comment.]
C
C               <character variable>  Returned string is of length long_.
C
C
C
C    purge_     Closes existing data files and clears tables and pointers.
C               [subroutine call]        
C
C____________________________________________________________________________
C
C
C
C   CIF CREATION TOOLS ("the put_ing commands")
C
C
C
C    pfile_     Create a file with the specified file name.
C               [logical function returned as .true. if the file is opened.
C               The value will be .false. if the file already exists.]
C
C               <file name>           Blank for use of currently open file
C
C
C
C    pdata_     Put a data block command into the created CIF. 
C               [logical function returned as .true. if the block is created.
C               The value will be .false. if the block name already exists.
C               Produces a save frame instead of a data block if the
C               variable saveo_ is true during the call.  No block duplicate
C               check is made for a save frame.]
C
C               <block name>
C
C
C
C    ploop_     Put a loop_ data name into the created CIF.             
C               [logical function returned as .true. if the invocation 
C               conforms with the CIF logical structure.  If pposval_ 
C               is non-zero, the "loop_" header is positioned to 
C               that column.  If pposnam_ is non-zero, the data name is 
C               positioned to that column.]
C
C               <data name>         If the name is blank on the first call
C                                   of a loop, only the "loop_" is placed.
C
C
C
C    pchar_     Put a character string into the created CIF.             
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary, 
C               AND, if the invocation conforms to the CIF logical structure.]
C
C               <data name>         If the name is blank, do not output name.
C
C               <character string>  A character string of MAXBUF chars or less.
C
C
C
C    pcmnt_     Puts a comment into the created CIF.
C               [logical function returned as .true.  The comment character
C               "#" should not be included in the string.  A blank comment
C               is presented as a blank line without the leading "#"].
C
C               <character string>  A character string of MAXBUF chars or less.
C
C
C    pnumb_     Put a single precision number and its esd into the created CIF.
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary, 
C               AND, if the invocation conforms to the CIF logical structure.
C               The number of esd digits is controlled by the variable
C               esdlim_]
C
C               <data name>         If the name is blank, do not output name.
C
C               <real variable>     Number to be inserted.
C
C               <real variable>     Esd number to be appended in parentheses.
C
C
C    pnumd_     Put a double precision number and its esd into the created CIF.
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary, 
C               AND, if the invocation conforms to the CIF logical structure.
C               The number of esd digits is controlled by the variable
C               esdlim_]
C
C               <data name>         If the name is blank, do not output name.
C
C               <double precision variable>  
C                                   Number to be inserted.
C
C               <double precision variable>  
C                                   Esd number to be appended in parentheses.
C
C
C
C    ptext_     Put a character string into the created CIF.             
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary, 
C               AND, if the invocation conforms to the CIF logical structure.]
C               ptext_ is invoked repeatedly until the text is finished. Only
C               the first invocation will insert a data name.
C
C               <data name>         If the name is blank, do not output name.
C
C               <character string>  A character string of MAXBUF chars or less.
C
C
C    prefx_     Puts a prefix onto subsequent lines of the created CIF.
C               [logical function returned as .true.  The second argument
C               may be zero to suppress a previously used prefix, or
C               greater than the non-blank length of the string to force
C               a left margin.  Any change in the length of the prefix 
C               string flushes pending partial output lines, but does _not_
C               force completion of pending text blocks or loops.
C               This function allows the CIF output functions to be used 
C               within what appear to be text fields to support annotation 
C               of a CIF. ]
C
C               <character string>  A character string of MAXBUF chars or less.
C
C               <integer variable>  The length of the prefix string to use.
C
C
C
C
C    close_     Close the creation CIF. MUST be used if pfile_ is used.
C               [subroutine call]
C
C
C____________________________________________________________________________
C
C
C
C....The CIF tool box also provides variables for data access control:
C 
C
C    alias_      Logical variable: if left .true. then all calls to
C                CIFtbx functions may use aliases of data item names.
C                The preferred synonym from the dictionary will be
C                subsituted internally, provided aliased data names were
C                supplied by an input dictionary (via dict_).  The
C                default is .true., but alias_ may be set to .false.
C                in an application.
C
C    aliaso_     Logical variable: if set .true. then cif output 
C                routines will convert aliases to the names to preferred
C                synonyms from the dictionary.  The default is .false., but
C                aliaso_ may be set to .true. in an application.  The
C                setting of aliaso_ is independent of the setting of
C                alias_.
C
C    align_      Logical variable signals alignment of loop_ lists during
C                the creation of a CIF. The default is .true.
C
C    bloc_       Character*(NUMCHAR) variable: the current block name.
C
C    dictype_    Character*(NUMCHAR) variable: the precise data type code
C                (see test_)
C
C    diccat_     Character*(NUMCHAR) variable: the category (see test_)
C
C    dicname_    Character*(NUMCHAR) variable: the root alias (see test_)
C
C    esdlim_     Integer variable:  Specifies the upper limit of esd's
C                produced by pnumb_, and, implicitly, the lower limit.
C                The default value is 19, which limits esd's to the range
C                2-19.  Typical values of esdlim_ might be 9 (limiting
C                esd's to the range 1-9), 19, or 29 (limiting esd's
C                to the range 3-29)
C
C    file_       Character*(MAXBUF) variable: the filename of the current file.
C
C    line_       Integer variable: Specifies the input/output line limit
C                for processing a CIF. The default value is 80 characters.
C                This may be set by the program. The max value is MAXBUF
C                which has a default value of 200.
C
C    list_       Integer variable: the loop block number (see test_).
C
C    long_       Integer variable: the length of the data string in strg_.
C
C    longf_      Integer variable: the length of the filename in file_.
C
C    loop_       Logical variable signals if another loop packet is present.
C
C    pposdec_    Integer variable giving the position of the decimal point
C                for the next number to be written.
C
C    pposend_    Integer variable giving the ending column of the next
C                number or quoted character value to be written.  Used to
C                pad with zeros or blanks.
C
C    pposnam_    Integer variable giving the starting column of the next
C                name or comment or data block to be written.
C
C    pposval_    Integer variable giving the starting column of the next
C                data value to be written by pchar_, pnumb_ or pnumd_.
C                Also used to set the position of the initial "loop_"
C                in a ploop_ call or to set the position of a terminal "save_"
C                for a save frame in a pdata_ call for which saveo_ is .true.
C
C    posdec_     Integer variable giving the position of the decimal point
C                for the last number read.
C
C    posend_     Integer variable giving the ending column of the last
C                data value read, not including a terminal quote.
C
C    posnam_     Integer variable giving the starting column of the last
C                name or comment or data block read.
C
C    posval_     Integer variable giving the starting column of the last
C                data value read.  Also reports the column of the
C                terminal "save_" of a save frame.
C
C    pquote_     Character variable giving the quotation symbol to be
C                used for the next string written.
C
C    precn_      Integer variable:  Reports the record number of the last
C                line written to the output cif.  Set to zero by init_.  Also
C                set to zero by pfile_ and close_ if the output cif file name
C                was not blank.
C
C    ptabx_      Logical variable signals tab character expansion to blanks 
C                during the creation of a CIF. The default is .true.
C
C    quote_      Character variable giving the quotation symbol found
C                delimiting the last string read.
C
C    recn_       Integer variable:  Reports the record number of the last
C                line read from the direct access copy of the input cif.
C
C    save_       Logical variable signals that the current data block
C                is actually a save-frame (.true. for a save-frame).
C
C    saveo_      Logical variable signals that the output data block from
C                pdata_ is actually a save-frame (.true. for a save-frame).
C
C    strg_       Character*(MAXBUF) variable: the current data item.
C
C    tabl_       Logical variable signals tab-stop alignment of output 
C                during the creation of a CIF. The default is .true.
C
C    tabx_       Logical variable signals tab character expansion to blanks 
C                during the reading of a CIF. The default is .true.
C
C    text_       Logical variable signals if another text line is present.
C
C    type_       Character*4 variable: the data type code (see test_).
C
C
C
C_____________________________________________________________________________
C
C
C >>>>>> Set the device numbers.
C
         function init_(devcif,devout,devdir,deverr)
C
         logical   init_
         include   'ciftbx.sys'
         integer   devcif,devout,devdir,deverr
         integer   ii,kdig
         real      ytest
         double precision ztest         
C
         init_=.true.
         cifdev=devcif
         outdev=devout
         dirdev=devdir
         errdev=deverr
         recn_=0
         precn_=0
C
C        recompute decimal single precision precision
C        This is found by computing the smallest power of
C        10 which, when added to 1, produces a change
C        and then backing off by 1
C
         decprc = .1
         do ii = 1,6
         ytest = 1.+decprc/10.
         if (ytest.eq.1.) go to 100
         decprc = decprc/10.
         enddo
100      continue
         decprc=decprc*10.
C
C        recompute decimal double precision precision
C
         kdig = 1
         dpprc = .1D0
         do ii = 1,15
         ztest = 1.D0+dpprc/10.
         if (ztest.eq.1.D0) go to 200
         dpprc = dpprc/10.D0
         kdig = kdig+1
         enddo
200      continue
         dpprc=dpprc*10.D0
         write(ndpfmt,'(5h(d30.,i2,1h))') kdig-1
C
C        recompute decimal single precision minimum power of ten
C
         decmin = .1
         do ii = 1,37
         ytest = decmin/10.
         if (ytest.eq.0.) go to 300
         decmin = decmin/10.
         enddo
300      continue
C
C        recompute decimal double precision minimum power of 10
C        and its log base 10 (minexp)
C
         dpmin = .1D0
         minexp = -1
         do ii = 1,307
         ztest = dpmin/10.
         if (ztest.eq.0.D0) go to 400
         dpmin = dpmin/10.D0
         minexp = minexp-1
         enddo
400      continue
         call clearfp
         return
         end
C
C
C
C
C
C >>>>>> Read a CIF dictionary and prepare for checks
C
         function dict_(fname,checks)
C
         logical   dict_
         logical   ocif_
         logical   data_
         logical   char_
         integer   lastnb
         include  'ciftbx.sys'
         character locase*(MAXBUF)
         character fname*(*),checks*(*)
         character temp*24,codes(4)*5,name*(MAXBUF),bxname*(NUMCHAR)
         character bcname*(NUMCHAR),biname*(NUMCHAR),bname*(NUMCHAR)
         character baname*(NUMCHAR),ganame*(NUMCHAR),btname*(NUMCHAR)
         character batag*(NUMCHAR)
         integer   lbcname,lbaname,lbtname,lbname
         integer   kdict,kadict,ifind,jfind,iafind
         integer   i,j,nmatch,mycat,ksmatch,ii
C
C        Control flags for matching categories, names and types
C
C        icloop is the loop number of the block for the
C        current category
C        ictype is the type of the current category
C          0 - none found yet
C          1 - _item.category.id
C          2 - _category
C          3 - _category.id
C        inloop is the loop number of the block for the
C        current name
C        intype is the type of the current name
C          0 - none found yet
C          1 - _item.name
C          2 - _name
C        ialoop is the loop number of the block for the
C        current alias
C        iatype is the type for the current alias
C          0 - none found yet
C          1 - _item_aliases.alias_name
C        itloop is the loop number of the block for the
C        current type
C        ittype is the type of the current type
C          0 - none found yet
C          1 - _item_type.code
C          2 - _type
C
         integer icloop,ictype,inloop,intype,ialoop,iatype,
     * itloop,ittype
C
         character*4 map_type(12),map_to(12),mapped
         character*(NUMCHAR) dt(2),ct(3),nt(2),at(1),tt(2)
         data map_type
     *   /'floa','int ','yyyy','symo','ucha','ucod','name','idna',
     *    'any ','code','line','ulin'/
         data map_to
     *   /'numb','numb','char','char','char','char','char','char',
     *    'char','char','char','char'/
         data dt
     *      /'_dictionary.title               ',
     *       '_dictionary_name                '/
         data ct
     *      /'_item.category_id               ',
     *       '_category                       ',
     *       '_category.id                    '/
         data nt
     *      /'_item.name                      ',
     *       '_name                           '/
         data at
     *      /'_item_aliases.alias_name        '/
         data tt
     *      /'_item_type.code                 ',
     *       '_type                           '/
C
         data codes /'valid','dtype','reset','close'/
C
C....... Are the codes OK
C
         temp=checks
         i=0         
120      i=i+1
         if(i.ge.24)                 goto 190
         if(temp(i:i).eq.' ')        goto 120
         do 150 j=1,4
         if(temp(i:i+4).eq.codes(j)) goto 170
150      continue
         dict_=.false.
         goto 500
170      i=i+4
         if(j.eq.1) vcheck='yes'
         if(j.eq.2) tcheck='yes'
         if(j.eq.3) then
           vcheck = 'no '
           tcheck = 'no '
           goto 170
         endif
         if(j.eq.4) then
           vcheck = 'no '
           tcheck = 'no '
           ndcname = 0
           ndict = 0
           if(nname.gt.0) then
           do 180 i = 1,nname
             dtype(i)=' '
             dxtyp(i)=' '
             cindex(i)=0
             ddict(i)=0
180        continue
           endif
           dict_=.true.
           goto 500
         endif
         goto 120
C
C        if no category names have been loaded, clean up
C        the hash table for dictionary category names
C
190      if(ndcname.eq.0) then
           call hash_init(dcname,dcchain,NUMDICT,ndcname,dchash,
     *     NUMHASH)
         endif
C
C        if no dictionary names have been loaded, clean up
C        the hash table for dictionary names
C
         if(ndict.eq.0) then
           call hash_init(dicnam,dicchain,NUMDICT,ndict,dichash,
     *     NUMHASH)
         endif
C
C....... Open and store the dictionary
C
         dict_=.true.
         if(fname.eq.' ')            goto 500
         if(nname.gt.0) call err(' Dict_ must precede ocif_')
         dict_=ocif_(fname)
         if(.not.dict_)              goto 500
         dictfl='yes'
C
C....... Loop over data blocks; extract _name's, _type etc.
C
200      if(.not.data_(' '))         goto 400
         if(bloc_(1:1).eq.'_') then
           bname=locase(bloc_)
         else
           bname='_'//locase(bloc_)
         endif
         lbname=lastnb(bname)
C
C        see if this is a dictionary defining block
C
         do i = 1,2
           if(char_(dt(i),name)) goto 200
         enddo
C
Cdbg     WRITE(6,*) ndict,bloc_
C
C        Analyze loop structure for categories, names and types
C
C
C        initalize loop info
C
         icloop = -1
         inloop = -1
         ialoop = -1
         itloop = -1
         ictype = 0
         intype = 0
         iatype = 0
         ittype = 0
         bcname = ' '
         lbcname = 1
         baname = ' '
         batag = ' '
         lbaname = 1
         btname = ' '
         lbtname = 1
         biname=bloc_
         mycat=0
         loop_=.false.
         loopnl=0
         nmatch=0
         ksmatch=0
C
C        Process categories
C
         do i = 1,3
           if(char_(ct(i),name)) then
             if(ictype.ne.0)
     *         call warn(' Multiple DDL 1 and 2 category definitions ')
             ictype = i
             if(loop_) icloop = loopnl
             bcname=locase(name(1:long_))
             lbcname=long_
             call hash_store(bcname,
     *         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
             if(mycat.eq.0) then
               call err(' Dictionary category names > NUMDICT ')
             endif
C
C            if this is not a loop of categories, we expect a match
C            against the block name
C
             if(.not.loop_) then
               if(ictype.eq.1) then
                 if(bname(1:lbcname+2).ne.
     *            '_'//bcname(1:lbcname)//'.') then
                 call warn(' Category id does not match block name')
                 endif
               else
                 if(ictype.eq.2) then
                   if(bcname.ne.'dictionary_definition') then
                   if(bname(1:lbcname+2).ne.
     *               '_'//bcname(1:lbcname)//'_') then
                   if(bname(1:lbcname+2).ne.
     *               '_'//bcname(1:lbcname)//' ') then
                   call warn(' Category id does not match block name')
                   endif
                   endif
                   endif
                 endif
               endif
             endif
           endif
           loop_ = .false.
           loopnl = 0
         enddo
C
C        Process names
         do i = 1,2
         if(char_(nt(i),name)) then
           if(intype.ne.0)
     *       call warn(' Multiple DDL 1 and 2 name definitions ')
           intype = i
           bxname=locase(name(1:long_))
           if(loop_) inloop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         enddo
         if(intype.eq.0.and.ictype.ne.3)
     *     call warn (' No name defined in block')
         loop_ = .false.
         if(char_(at(1),name)) then
           iatype=1
           baname = locase(name(1:long_))
           batag = name(1:long_)
           lbaname = long_
           if(loop_) ialoop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         if(ictype.ne.3) then
           do i=1,2
             if(char_(tt(i),name)) then
               if(ittype.ne.0)
     *           call warn(' Multiple DDL 1 and 2 type definitions ')
               ittype = i
               btname = locase(name(1:long_))
               if(loop_) itloop = loopnl
             endif
             loop_ = .false.
             loopnl=0
           enddo
         endif
C
C        Now test for consistent combinations
C
         if(inloop.ne.-1) then
           if(icloop.ne.-1.and.icloop.ne.inloop)
     *       call warn(
     *       ' Categories and names in different loops')
           if(iatype.ne.0.and.ialoop.ne.inloop) then
             if(ialoop.eq.-1) then
               if(bxname.ne.bname)
     *          call warn(
     *         ' One alias, looped names, linking to first')
             else
               call warn(
     *         ' Aliases and names in different loops '
     *         //' only using first alias ')
             endif
           endif
           if(itloop.ne.-1.and.itloop.ne.inloop)
     *       call warn(
     *       ' Types and names in different loops')
         else
           if(icloop.ne.-1)
     *       call warn(
     *         ' Multiple categories for one name')
           if(itloop.ne.-1)
     *       call warn(
     *         ' Multiple types for one name')
         endif
C
C        This is the main loop
C
         if(intype.eq.0) go to 200
250      if(.not.char_(nt(intype),name)) goto 200
         kdict=ndict+1
         call hash_store(locase(name(1:long_)),dicnam,dicchain,
     *     NUMDICT,ndict,dichash,NUMHASH,ifind)
         if(ifind.eq.0) call err(' Cifdic names > NUMDICT')
         if(ifind.eq.kdict)dictag(ifind)=name(1:long_)
         if(dicnam(ifind).eq.bname) nmatch=ifind
         if(dicnam(ifind)(1:lbname).eq.bname) ksmatch=ifind
Cdbg     if(dicnam(ifind).ne.bname)
Cdbg *   call warn (' Name mismatch: '//dicnam(ifind)//bname)
         if(inloop.ge.0)then
C
C          We are in a loop of names.  If it is the same loop as
C          for categories, we need to extract the matching category
C
           if(inloop.eq.icloop) then
             mycat=0
             if(char_(ct(ictype),name)) then
               bcname=locase(name(1:long_))
               lbcname=long_
               call hash_store(bcname,
     *         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call err(' Dictionary category names > NUMDICT ')
               endif
             endif
           endif
C
C          If it is the same loop as for types, we need to extract
C          the matching type
C
           if(inloop.eq.itloop) then
             btname=' '
             if(char_(ct(ittype),name)) then
               btname=locase(name(1:long_))
               lbtname=long_
             endif
           endif
C
C          If it is the same loop as for aliases, we need to extract
C          the matching alias
C
           if(inloop.eq.ialoop) then
             baname=' '
             batag=' '
             if(char_(at(1),name)) then
               baname = locase(name(1:long_))
               batag = name(1:long_)
               lbaname = long_
             endif
           endif
         endif
C
C        now we have a name stored in dicnam at location ifind
C        the index of the category in mycat, the type in btname,
C        the alias in baname
C
C        First verify match between the name and category, if
C        we have one, or extract from the block name
C
         if (mycat.eq.0) then
         if (dcindex(ifind).eq.0) then
           if (dicnam(ifind).eq.bloc_) then
             call excat(dicnam(ifind),bcname,lbcname)
Cdbg         call warn(' Extracting category name from block name '
Cdbg *       //bloc_(1:max(1,lastnb(bloc_))))
             if(bcname(1:1).ne.' ') then
               ictype = 1
               call hash_store(bcname,
     *         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call err(' Dictionary category names > NUMDICT ')
               endif
             else
               call warn(' No category defined in block ' 
     *       //bloc_(1:max(1,lastnb(bloc_)))//' and name '
     *       //dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))
     *       //' does not match')
             endif
           endif
         endif
         else
         if (bcname(1:lbcname).ne.'dictionary_definition') then
           if (dicnam(ifind)(1:lbcname+1).ne.'_'//bcname(1:lbcname)
     *        .or.( dicnam(ifind)(lbcname+2:lbcname+2).ne.'_' .and.
     *          dicnam(ifind)(lbcname+2:lbcname+2).ne.'.' .and.
     *          dicnam(ifind)(lbcname+2:lbcname+2).ne.' ' )) then
                call warn(' Item name '//
     *          dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))//' '//
     *       ' does not match category name '//bcname(1:lbcname))
           endif
         endif
         endif
C
C        We will need the type in what follows.  cifdic.m96 defines
C        some higher level types.  We map them to primitive types
C
         mapped = btname(1:4)
         do i = 1,12
           if (btname(1:4).eq.map_type(i)) mapped = map_to(i)
         enddo
         if (mapped.ne.'char' .and.
     *       mapped.ne.'text' .and.
     *       mapped.ne.'    ' .and.
     *       mapped.ne.'null' .and.
     *       mapped.ne.'numb' ) then
             if (tcheck .eq. 'yes') call warn (' Item type '//
     *       btname(1:max(1,lastnb(btname)))//' not recognized')
         endif
C
C        There are two cases to consider, one if the name is new to
C        the dictionary, the other, if it is not
C
         if(ifind.eq.kdict) then
           aroot(ifind)=0
           alias(ifind)=0
           dcindex(ifind)=mycat
           dictyp(ifind)=mapped
           dicxtyp(ifind)=btname
         else
           if(dcindex(ifind).ne.mycat) then
             if(dcindex(ifind).eq.0) then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=ifind
255            continue
               dcindex(jfind)=mycat
               jfind=alias(jfind)
               if(jfind.ne.0) goto 255
             else
               if(mycat.ne.0.and.
     *           (vcheck.eq.'yes'.or.tcheck.eq.'yes'))
     *           call warn(' Attempt to redefine category for item')
             endif
           endif
           if(dictyp(ifind).ne.mapped .or.
     *       dicxtyp(ifind).ne.btname) then
             if(dictyp(ifind).eq.' ') then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=ifind
256            continue
               dictyp(jfind)=mapped
               dicxtyp(jfind)=btname
               jfind=alias(jfind)
               if(jfind.ne.0) go to 256
             else
               if(mapped.ne.' '.and.tcheck.eq.'yes')
     *           call warn(' Attempt to redefine type for item')
             endif
           endif
         endif
C
C        now deal with alias, if any.
C
         if(baname.ne.' ') then
           kadict=ndict+1
           call hash_store(baname(1:lbaname),dicnam,dicchain,
     *     NUMDICT,ndict,dichash,NUMHASH,iafind)
           if(iafind.eq.0) call err(' Cifdic names > NUMDICT')
           if(iafind.eq.kadict) then
             dictag(iafind)    =batag
             aroot(iafind)     =aroot(ifind)
             if(aroot(iafind).eq.0) aroot(iafind)=ifind
             alias(iafind)     =0
             alias(ifind)      =iafind
             dcindex(iafind)   =dcindex(ifind)
             dictyp(iafind)    =dictyp(ifind)
             dicxtyp(iafind)   =dicxtyp(ifind)
           else
             if(aroot(iafind).ne.0) then
               if(aroot(iafind).eq.ifind .or.
     *           aroot(iafind).eq.aroot(ifind)) then
                 call warn(' Duplicate definition of same alias')
               else
                 call warn(' Conflicting definition of alias')
               endif
             else
               if((dcindex(iafind).eq.0.or.
     *           dcindex(iafind).eq.dcindex(ifind)).and.
     *           (dictyp(iafind).eq.' '.or.
     *           (dictyp(iafind).eq.dictyp(ifind) .and.
     *            dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
               endif
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               alias(ifind)      =iafind
             endif
           endif
         endif
         if(inloop.ge.0) then
           baname = ' '
           batag = ' '
         endif
C
         if(inloop.ge.0.and.loop_) go to 250
         if(nmatch.eq.0) then
         if (ksmatch.eq.0.or.inloop.lt.0) then
         call warn(' No name in the block matches the block name')
         endif
         endif
C
C        check for aliases
C        we execute this loop only in the case of unlooped name
C        with looped alias
C
         if(inloop.lt.0.and.ialoop.ge.0) then
           loop_=.false.
           loopnl=0 
           ganame=baname
260        if(.not.char_(at(iatype),name)) goto 200
           baname=locase(name(1:long_))
           batag=name(1:long_)
           lbaname=long_
           if(baname.eq.ganame) then
             if(loop_) go to 260
             go to 200
           endif
           if(baname.ne.' ') then
             kadict=ndict+1
             call hash_store(baname(1:lbaname),dicnam,dicchain,
     *       NUMDICT,ndict,dichash,NUMHASH,iafind)
             if(iafind.eq.0) call err(' Cifdic names > NUMDICT')
             if(iafind.eq.kadict) then
               dictag(iafind)    =batag
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               alias(iafind)     =0
               alias(ifind)      =iafind
               dcindex(iafind)   =dcindex(ifind)
               dictyp(iafind)    =dictyp(ifind)
               dicxtyp(iafind)   =dicxtyp(ifind)
               ifind=iafind
             else
               if(aroot(iafind).ne.0) then
                 if(aroot(iafind).eq.ifind .or.
     *             aroot(iafind).eq.aroot(ifind)) then
                   call warn(' Duplicate definition of same alias')
                 else
                   call warn(' Conflicting definition of alias')
                 endif
               else
                 if((dcindex(iafind).eq.0.or.
     *           dcindex(iafind).eq.dcindex(ifind)).and.
     *           (dictyp(iafind).eq.' '.or.
     *           (dictyp(iafind).eq.dictyp(ifind) .and.
     *            dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
                 ifind=iafind
                 endif
                 aroot(iafind)     =aroot(ifind)
                 if(aroot(iafind).eq.0) aroot(iafind)=ifind
                 alias(ifind)      =iafind
               endif
             endif
           endif
           if(loop_) go to 260
         endif
         go to 200
C
400      bloc_=' '
         if (ndcname.ne.0) then
         do ii = 1,ndict
         if (aroot(ii).eq.0.and.dcindex(ii).eq.0)
     *     call warn(' No category specified for name '//
     *       dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
         enddo
         endif
         do ii = 1,ndict
         if (dicxtyp(ii).eq.' ') then
           dicxtyp(ii) = 'null'
           dictyp(ii) = 'null'
           if (tcheck.eq.'yes')
     *       call warn(' No type specified for name '//
     *         dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
         endif
         enddo
         close(dirdev)
         nrecd=0
         dictfl='no '
500      continue
         if(tcheck.eq.'yes') vcheck='yes'
Cdbg     WRITE(6,'(i5,3x,a,2x,a)') (i,dicnam(i),dictyp(i),i=1,ndict)
         return
         end
C
C
C
C
C
C >>>>>> Find position of last non_blank in a string
C
         function lastnb(str)
C
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) str
         integer lenn,ii
         lenn = len(str)
         do 100 ii=lenn,1,-1
         if(str(ii:ii).eq.' ') goto 100
         if(str(ii:ii).ne.tab) goto 120
100      continue
         ii=1
120      lastnb = ii
         return
         end
C
C
C
C
C
C >>>>>> Extract the item.category_id from a save frame name
C
         subroutine excat(sfname,bcname,lbcname)
C
         character*(*) sfname,bcname
         integer lbcname,ii,ic,lastnb,lenn
C
C        Note that this logic works only for item.category_id
C        not for category.id
C
         lenn = lastnb(sfname)
         bcname = ' '
         lbcname = 1
         if (lenn.eq.0.or.sfname(1:1).ne.'_') return
         do ii = 1,lenn-2
         ic = 1+lenn-ii
         if (sfname(ic:ic).eq.'.') then
           bcname = sfname(2:ic-1)
           lbcname = ic-2
           return
         endif
         enddo
         return
         end
C
C
C
C
C
C >>>>>> Open a CIF and copy its contents into a direct access file.
C
         function ocif_(fname)
C
         logical   ocif_
         integer   lastnb
         include  'ciftbx.sys'
         logical   test
         character fname*(*)
         integer   case,i
C
         save_=.false.
         jchar=MAXBUF
         lastch=0
         if(line_.gt.MAXBUF) call err(' Input line_ value > MAXBUF')
         if(nrecd.ne.0) close(dirdev)
         nrecd=0
         lrecd=0
         case=ichar('a')-ichar('A')
         tab=char(05)
         if(case.lt.0) goto 100
         tab=char(09)
         bloc_=' '
C
C....... Make sure the CIF is available to open
C
100      file_=fname
         do 120 i=1,MAXBUF
         if(file_(i:i).eq.' ') goto 140
120      continue
140      longf_=i-1
         if (longf_.gt.0) then
           inquire(file=file_(1:longf_),exist=test)
           ocif_=test
           if(.not.ocif_)      goto 200
         else
           file_ = ' '
           longf_ = 1
           ocif_ = .true.
         endif
C
C....... Open up the CIF and a direct access formatted file as scratch
C
         if (file_(1:1).ne.' ')
     *   open(unit=cifdev,file=fname,status='OLD',access='SEQUENTIAL',
     *                    form='FORMATTED')
         open(unit=dirdev,status='SCRATCH',access='DIRECT',
     *                    form='FORMATTED',recl=MAXBUF)
C
C....... Copy the CIF to the direct access file
C
160      read(cifdev,'(a,a)',end=180) buffer
         nrecd=nrecd+1
         irecd=nrecd
         if (lastnb(buffer(1:MAXBUF)).gt.line_)
     *      call warn(' Input line length exceeds line_')
         write(dirdev,'(a)',rec=nrecd) buffer
Cdbg     WRITE(6,'(i5,1x,a)') nrecd,buffer(1:70)
         goto 160
C
180      lrecd=0
         jrecd=0
         jrect=-1
         irecd=0
         recn_=0
         if (file_(1:1).ne.' ') close(cifdev)
200      return
         end
C
C
C
C
C
C >>>>>> Close off direct access file of the current CIF
C         and reset all data name tables and pointers       
C
         subroutine purge_
C
         include   'ciftbx.sys'
C
         if(nrecd.ne.0) close(dirdev)
         recn_=0
         save_=.false.
         jchar=MAXBUF
         lastch=0
         nrecd=0
         lrecd=0
         irecd=0
         nname=0
         nhash=0
         iname=0
         loopct=0
         loopnl=0
         loop_=.false.
         text_=.false.
         return
         end
C
C
C
C
C
C >>>>>> Store the data names and pointers for the requested data block
C
         function data_(name) 
C
         logical   data_
         integer   lastnb
         include  'ciftbx.sys'
         character name*(*),flag*4,temp*(NUMCHAR),ltype*4
         character ctemp*(NUMCHAR)
         character locase*(MAXBUF),isbuf*(MAXBUF)
         integer   ndata,idata,nitem,npakt,i,ii,j,k,kchar,krecd
         integer   fcatnum,lctemp,isrecd,isjchr,islast
         integer   pnname,itpos,ipp,ipj
C
         jchar=MAXBUF
         nname=0
         ndata=0
         nhash=0
         nitem=0
         idata=0
         iname=0
         loopct=0
         loopnl=0
         ltype=' '
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         data_=.false.
         loop_=.false.
         text_=.false.
         do ii = 1,MAXBOOK
         ibkmrk(1,ii)=-1
         enddo
         irecd=lrecd
         lrecd=nrecd
         if(name(1:1).ne.' ') irecd=0
         call hash_init(dname,dchain,NUMBLOCK,nname,dhash,
     *     NUMHASH)
         call hash_init(cname,cchain,NUMBLOCK,ncname,chash,
     *     NUMHASH)
         isrecd=irecd
         isjchr=jchar
         islast=lastch
         isbuf=' '
         if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
C
C....... Find the requested data block in the file
C
100      call getstr
         isjchr=jchar
         if(irecd.ne.isrecd) then
           isrecd=irecd
           islast=lastch
           isbuf=' '
           if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         endif
         if(type_.eq.'fini')           goto 500
         if(type_.ne.'text')           goto 120
110      call getlin(flag)       
         if(buffer(1:1).ne.';')        goto 110
         jchar=2
         goto 100
120      continue
         if(type_.eq.'save') then
           if(long_.lt.6) then
             if(.not.save_)
     *         call err(' Save frame terminator found out of context ')
             save_=.false.
             goto 100
           else
             if(save_)
     *         call err(' Prior save frame not terminated ')
             save_=.true.
             if(name.eq.' ')          goto 150
             if(strg_(6:long_).ne.name) goto 100
             goto 150
           endif
         endif
         if(type_.ne.'data')          goto 100
         if(name.eq.' ')              goto 150
         if(strg_(6:long_).ne.name)   goto 100
150      data_=.true.
         bloc_=strg_(6:long_)
         itpos=jchar-long_
         if(tabx_) then
         itpos=0
         do ipp=1,jchar-long_
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         endif
         posnam_=itpos
C
C....... Get the next token and identify
C
200      call getstr
Cdbg     if(dictfl.eq.'no ')
Cdbg *    WRITE(6,*) ltype,type_,loop_,nitem,ndata,idata,iname,nname
C
         if(ltype.ne.'name')                goto 201
         if(type_.eq.'numb')                goto 203
         if(type_.eq.'char')                goto 203
         if(type_.eq.'text')                goto 203
         if(type_.eq.'null')                goto 203
         if(type_.eq.'name'.and.loop_)      goto 204
         call err(' Illegal tag/value construction')
201      if(ltype.ne.'valu')                goto 204
         if(type_.eq.'numb')                goto 202
         if(type_.eq.'char')                goto 202
         if(type_.eq.'text')                goto 202
         if(type_.eq.'null')                goto 202
         goto 204
202      if(nitem.gt.0)                     goto 205
         call err(' Illegal tag/value construction')
203      ltype='valu'
         goto 205
204      ltype=type_
C
205      if(type_.eq.'name')           goto 206
         if(type_.eq.'loop')           goto 210
         if(type_.eq.'data')           goto 210
         if(type_.eq.'save')           goto 210
         if(type_.ne.'fini')           goto 220
206      if(loop_)                     goto 270
210      if(nitem.eq.0)                goto 215
C
C....... End of loop detected; save pointers
C
         npakt=idata/nitem
         if(npakt*nitem.ne.idata) call err(' Item miscount in loop')
         loopni(loopct)=nitem
         loopnp(loopct)=npakt
         nitem=0
         idata=0
215      if(type_.eq.'name')           goto 270
         if(type_.eq.'data')           goto 300
         if(type_.eq.'save')           goto 300
         if(type_.eq.'fini')           goto 300
C
C....... Loop_ line detected; incr loop block counter
C
         loop_=.true.
         loopct=loopct+1
         if(loopct.gt.NUMLOOP) call err(' Number of loop_s > NUMLOOP')
         loorec(loopct)=irecd
         loopos(loopct)=jchar-long_
         if(quote_.ne.' ') loopos(loopct)=jchar-long_-1
         itpos=0
         do ipp=1,loopos(loopct)
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         loopox(loopct)=itpos
         goto 200
C
C....... This is the data item; store char position and length
C
220      if(loop_ .and. nitem.eq.0)
     *   call err(' Illegal tag/value construction')
         loop_=.false.
C
         i=nname
         if(nitem.gt.0) i=i-nitem+mod(idata,nitem)+1
         if(i.lt.1) call err(' Illegal tag/value construction')
         if(dtype(i).ne.'test')       goto 223
         if(dictfl.eq.'yes')          goto 223
         if(tcheck.eq.'no ')          goto 223
C>>>>    if(long_.eq.1.and.strg_(1:1).eq.'?') goto 223
C>>>>    if(long_.eq.1.and.strg_(1:1).eq.'.') goto 223
         if(type_.eq.'null')          goto 223
         if(type_.eq.'numb')          goto 223
         call warn( ' Numb type violated  '//dname(i))
223      if(nitem.le.0)               goto 224
         idata=idata+1
         if(dtype(i).eq.'null') dtype(i)=type_
         if(dtype(i).eq.'numb' .and.
     *     (type_.eq.'char'.or.type_.eq.'text')) dtype(i)='char'
224      if(nname.eq.ndata)           goto 230
         ndata=ndata+1
         if(iloop(ndata).gt.1)        goto 225
         krecd=irecd
         kchar=jchar-long_-1
         if(quote_.ne.' ')kchar=kchar-1
225      continue
         if(dtype(ndata).eq.'    ') dtype(ndata)=type_
         drecd(ndata)=krecd
         dchar(ndata)=kchar
         if(nloop(ndata).gt.0)        goto 230
         nloop(ndata)=0
         iloop(ndata)=long_
C
C....... Skip text lines if present
C
230      if(type_.ne.'text')           goto 200
         if(nloop(ndata).eq.0) dchar(ndata)=0
         if(nloop(ndata).eq.0) iloop(ndata)=long_
250      call getlin(flag)
         if(buffer(1:1).eq.';') then
           jchar=2
           goto 200
         endif
         if(flag.eq.'fini') call err(' Unexpected end of data')
         goto 250
C
C....... This is a data name; store name and loop parameters
C
270      temp=locase(strg_(1:long_))
         k=0
         if(dictfl.ne.'yes' .and. ndict.gt.0) then
           call hash_find(temp,
     *       dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,k)
           if(k.ne.0) then
             if(alias_ .and. aroot(k).ne.0) temp=dicnam(aroot(k))
           endif
         endif
         pnname=nname
         call hash_store(temp,
     *   dname,dchain,NUMBLOCK,nname,dhash,
     *     NUMHASH,j)
         if(j.eq.pnname+1) then
           dtag(j)=strg_(1:long_)
           if(k.ne.0) dtag(j)=dictag(k)
           trecd(j)=irecd
           tchar(j)=jchar-long_
           if(quote_.ne.' ') tchar(j)=jchar-long_-1
           itpos=0
           do ipp=1,tchar(j)
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           xchar(j)=itpos
         endif
         if(j.eq.0)
     *     call err(' Number of data names > NUMBLOCK')
         if(k.ne.0)temp=dicnam(k)
         if(j.ne.pnname+1) then
           call warn(' Duplicate data item '//
     *     temp(1:max(1,lastnb(temp))))
           goto 200  
         endif
         dtype(nname)=' '
         dxtyp(nname)=' '
         cindex(nname)=0
         ddict(nname)=0
         ctemp='(none)'
         lctemp=6
C
         if(dictfl.eq.'yes' .or. vcheck.eq.'no ') goto 290
         j=k
         if(j.ne.0) then
           ddict(nname)=j
           cindex(nname)=dcindex(j)
           dxtyp(nname)=dicxtyp(j)
           dtype(nname)=dictyp(j)
           if(vcheck.eq.'no ')          goto 280
           if(dictyp(j).eq.'numb') then
             dtype(nname)='test'
           endif
           if(cindex(nname).ne.0) then 
             ctemp=dcname(cindex(nname))
             lctemp=lastnb(ctemp)
             goto 290
           endif   
           goto  280
         endif
         call warn(' Data name '//
     *               temp(1:max(1,lastnb(temp)))
     *               //' not in dictionary!')
280      call excat(temp,ctemp,lctemp)
         if (ctemp.eq.' '.or.'_'//ctemp.eq.temp) then
           ctemp = '(none)'
           lctemp= 6
           if (ndcname.ne.0.and.vcheck.eq.'yes')
     *       call warn(' No category defined for '
     *       //temp)
         else
           call hash_find(ctemp,
     *       dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,j)
           if(j.ne.0) then
             cindex(nname) = j
           else
             ipj=ncname
             call hash_store(ctemp(1:lctemp),
     *         cname,cchain,NUMBLOCK,ncname,chash,NUMHASH,j)
             if (j.eq.0)
     *         call err(' Number of categories > NUMBLOCK ')
             cindex(nname) = -j
             if (ndcname.gt.0.and.j.eq.ipj+1.and.vcheck.eq.'yes')
     *         call warn(' Category '//
     *         ctemp(1:lctemp)//' first implicitly defined in cif ')
           endif
         endif
C
290      lloop(nname)=0
         nloop(nname)=0
         iloop(nname)=0
         if (nitem.eq.0) fcatnum=cindex(nname)
         if(.not.loop_)               goto 200
         nitem=nitem+1
         if(nitem.gt.NUMITEM)
     *     call err(' Items per loop packet > NUMITEM')
         nloop(nname)=loopct
         iloop(nname)=nitem
         if (fcatnum.ne.cindex(nname)) then
           temp = '(none)'
           if (fcatnum.gt.0) temp=dcname(fcatnum)
           if (fcatnum.lt.0) temp=cname(-fcatnum)
           if (ctemp(1:lctemp).ne.temp(1:lastnb(temp)))
     *     call warn (' Heterogeneous categories in loop '//
     *     ctemp(1:lastnb(ctemp))//' vs '//
     *     temp(1:lastnb(temp)))
           fcatnum=cindex(nname)
         endif
         goto 200
300      continue
C
C....... Are names checked against dictionary?
C
         if(dictfl.eq.'yes')          goto 500
         if(vcheck.eq.'no '.or.ndict.eq.0) goto 500
         do 350 i=1,nname
         if(dtype(i).eq.'test') dtype(i)='numb'
350      continue
C
C....... End of data block; tidy up loop storage
C
500      lrecd=irecd-1
         if(type_.eq.'save'.and.long_.lt.6) then
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posval_=itpos
         endif
         irecd=isrecd
         jchar=isjchr
         lastch=islast
         recn_=irecd
         buffer=' '
         if(lastch.gt.0)buffer=isbuf(1:lastch)
         jrecd=irecd
         loop_=.false.
         loopct=0
         if(ndata.ne.nname) call err(' Syntax construction error')
C
Cdbg     WRITE(6,'(a)')
Cdbg *   ' data name                       type recd char loop leng'
Cdbg     WRITE(6,'(a,1x,a,4i5)') (dname(i),dtype(i),drecd(i),dchar(i),
Cdbg *              nloop(i),iloop(i),i=1,nname)
Cdbg     WRITE(6,'(3i5)') (i,loopni(i),loopnp(i),i=1,loopct)
C
         return
         end
C
C
C
C
C
C
C >>>>>> Get the attributes of data item associated with data name
C
         function test_(temp)
C
         logical   test_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  otestf*3
         character  locase*(MAXBUF)
C
         otestf=testfl
         testfl='yes'
         name=locase(temp)
         test_=.true.   
         if(otestf.eq.'no ')  goto 100
         if(name.eq.nametb)   goto 200
100      call getitm(name)        
200      list_ =loopnl
         if(type_.eq.'null') test_=.false.
         return
         end

C
C
C
C
C
C >>>>>> Set or Reference a bookmark
C
         function bkmrk_(mark)
C
         logical   bkmrk_
         include   'ciftbx.sys'
C
         integer   mark,ii,nitem
         character*4 flag
         bkmrk_=.true.
         if(mark.eq.0) then
           do ii=1,MAXBOOK
             if(ibkmrk(1,ii).lt.0)      goto 100
           enddo
           bkmrk_=.false.
           call warn(' More than MAXBOOK bookmarks requested')
           return
100        mark=ii
           ibkmrk(1,ii)=iname
           ibkmrk(2,ii)=irecd
           ibkmrk(3,ii)=jchar
           if(iname.gt.0) then
             ibkmrk(2,ii) = trecd(iname)
             ibkmrk(3,ii) = tchar(iname)
           endif
           ibkmrk(4,ii)=0
           if(iname.gt.0) then
             if(nloop(iname).ne.0.and.
     *         loopnl.eq.nloop(iname).and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               ibkmrk(2,ii)=looprd(1)
               ibkmrk(3,ii)=max(0,loopch(1)-1)
               ibkmrk(4,ii)=loopct
             endif
           endif
         else
           if(ibkmrk(1,mark).lt.0) then
             bkmrk_=.false.
             return
           endif
           iname=ibkmrk(1,mark)
           irecd=ibkmrk(2,mark)
           loopct=ibkmrk(4,mark)
           loop_=.false.
           text_=.false.
           loopnl=-1
           testfl='no '
           if(iname.gt.0) then
            if(nloop(iname).ne.0.and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               looprd(nitem+1)=ibkmrk(2,mark)
               loopch(nitem+1)=ibkmrk(3,mark)
               do ii = 1,nitem
                 lloop(ii+iname-iloop(iname))=loopct-1
               enddo
               loopct=loopct-1
               if(lloop(iname).gt.0) then
                 loop_=.true.
                 loopnl=nloop(iname)
               endif
             endif
           endif
           jchar=MAXBUF
           if(irecd.gt.0) then
             irecd=irecd-1
             call getlin(flag)
             jchar=ibkmrk(3,mark)
           endif
           ibkmrk(1,mark)=-1
           mark=0
         endif
         return
         end
C
C
C
C
C
C
C >>>>>> Find the location of the requested item in the CIF
C        The argument "name" may be a data item name, blank
C        for the next such item.  The argument "type" may be
C        blank for unrestricted acceptance of any non-comment
C        string (use cmnt_ to see comments), including loop headers,
C        "name" to accept only the name itself and "valu"
C        to accept only the value, or "head" to position to the
C        head of the CIF.  Except when the "head" is requested,
C        the position is  left after the data item provided.
C
         function find_(name,type,strg)
C
         logical   find_
         include   'ciftbx.sys'
         character  name*(*),type*(*),strg*(*),flag*4
         character  jjbuf*(MAXBUF)
         integer    jjchar,jjrecd,jjlast,jjlrec,jjjrec
C
         find_  = .false.
         strg   = ' '
         long_  = 0
         jjchar = jchar
         jjrecd = lrecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf  = ' '
         if(lastch.gt.0) jjbuf(1:lastch)=buffer(1:lastch)
         if(type.eq.'head') then
           lrecd = nrecd
           irecd=0
           jchar=MAXBUF+1
           call getlin(flag)
           if(flag.eq.'fini')       goto 300
           find_=.true.
           lrecd=jjlrec
           return
         endif
         if(name.ne.' ') then
           testfl='no '
           call getitm(name)
           if(iname.eq.0) goto 300
           if(type.eq.'valu') then
             list_=loopnl
             strg=strg_(1:long_)
             find_=.true.
             return
           endif
           if(type.eq.'name'.or.loopnl.eq.0) then
             irecd=trecd(iname)-1
             call getlin(flag)
             jchar=tchar(iname)
             posnam_=jchar+1
             call getstr
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           if(type.eq.' ') then
             irecd=loorec(loopnl)-1
             call getlin(flag)
             jchar=loopos(loopnl)
             call getstr
             posval_=loopos(loopnl)
             if(tabx_) posval_=loopox(loopnl)
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           call err(' Call to find_ with invalid arguments')
         endif
         if(name.eq.' ') then
200        call getstr
           if(type_.eq.'fini')      goto 300
           if(type.ne.' '.and.
     *      (type_.eq.'data'.or.type_.eq.'save'))   goto 300
           if(type.eq.'name'.and.type_.ne.'name')  goto 200
           if(type.eq.'valu'.and.
     *       type_.ne.'numb'.and.type_.ne.'text'
     *      .and.type_.ne.'char'.and.type_.ne.'null') goto 200
           find_=.true.
           strg=strg_(1:long_)
           if(type_.eq.'name') then
             posnam_=jchar-long_
           else
             posval_=jchar-long_
             if(quote_.ne.' ') posval_=posval_-1
           endif
           recn_=irecd
           return
         endif
    
C
C        Search failed, restore pointers
C
300      irecd  = jjrecd
         lastch = jjlast
         lrecd  = jjlrec
         jchar  = jjchar
         buffer = ' '
         if(lastch.gt.0)buffer(1:lastch)=jjbuf(1:lastch)
         jrecd  = jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_  = irecd
C           
         return
         end
C
C
C
C
C
C
C >>>>>> Get the next data name in the data block
C
         function name_(temp)
C
         logical    name_
         include   'ciftbx.sys'
         character  temp*(*)
C
         name_=.false.
         temp=' '
         iname=iname+1
         if(iname.gt.nname)  goto 100
         name_=.true.
         temp=dtag(iname)
         if(ddict(iname).ne.0) temp=dictag(ddict(iname))
100      return
         end
C
C
C
C
C
C
C >>>>>> Extract a number data item and its standard deviation
C        This version return single precision numbers
C
         function numb_(temp,numb,sdev)
C
         logical    numb_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  locase*(MAXBUF)
         real       numb,sdev
C
         name=locase(temp)
         if(testfl.eq.'no ')  goto 100
         if(name.eq.nametb)   goto 150
C
100      call getitm(name)
C
150      numb_=.false.
         if(type_.ne.'numb') goto 200
         numb_=.true.
         numb =numbtb
         if(sdevtb.ge.0.0) sdev=sdevtb
C
200      testfl='no '
         return
         end
C
C
C
C
C
C
C >>>>>> Extract a number data item and its standard deviation
C        This version returns double precision numbers
C
         function numd_(temp,numb,sdev)
C
         logical    numd_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  locase*(MAXBUF)
         double precision numb,sdev
C
         name=locase(temp)
         if(testfl.eq.'no ')  goto 100
         if(name.eq.nametb)   goto 150
C
100      call getitm(name)
C
150      numd_=.false.
         if(type_.ne.'numb') goto 200
         numd_=.true.
         numb =numbtb
         if(sdevtb.ge.0.0) sdev=sdevtb
C
200      testfl='no '
         return
         end
C
C
C
C
C
C
C >>>>>> Extract a character data item.
C
         function char_(temp,strg)
C
         logical    char_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  strg*(*),flag*4
         character  locase*(MAXBUF)
         integer    icpos,itpos,ixpos,ixtpos,ipp,iepos,ispos
C
         name=locase(temp)
         if(testfl.eq.'yes')    goto 100
         if(.not.text_)         goto 120
         if(name.ne.nametb)     goto 120
         char_=.false.
         text_=.false.
         strg=' '
         long_=0
         call getlin(flag)
         if(flag.eq.'fini')    goto 200
         if(buffer(1:1).eq.';') then
           jchar=2
           goto 200
         endif
         quote_=' '
         jchar=lastch+1
         long_=lastch
         strg_(1:long_)=buffer(1:long_)
         goto 150
C
100      if(name.eq.nametb)     goto 150
C
120      call getitm(name)
         if(type_.eq.'null') then
           char_=.false.
           text_=.false.
           strg_=' '
           long_=0
           goto 200
         endif
C
150      char_=.true.
         text_=.false.
         if(tabx_) then
           call detab
           icpos=jchar-long_
           if(quote_.ne.' ') icpos=icpos-1
           iepos=icpos+long_-1
           itpos=0
           do ipp=1,icpos
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           ispos=itpos
160        ixpos=index(buffer(icpos:iepos),tab)
           ixtpos=itpos+ixpos-1
           if(ixpos.gt.0.and.ixtpos.le.MAXBUF) then
             ixtpos=((ixtpos+7)/8)*8
             icpos=icpos+ixpos
             itpos=ixtpos+1
             if(icpos.le.iepos) goto 160
           else
           strg =
     *       bufntb(ispos:min(MAXBUF,itpos+iepos-icpos))
           long_=min(MAXBUF,itpos+iepos-icpos)-ispos+1
           if(ispos.eq.1.and.strg(1:1).eq.';')
     *       strg(1:1) = ' '
           endif
         else
           strg=strg_(1:long_)
         endif
         if(type_.eq.'char')   goto 200
         char_=.false.
         if(type_.ne.'text')   goto 200
         char_=.true.
         call getlin(flag)
         jchar=MAXBUF+1
         if(flag.eq.'fini')    goto 200
         if(buffer(1:1).eq.';')then
           jchar=2
           goto 200
         endif
         irecd=irecd-1
         text_=.true. 
C
200      testfl='no '
         return
         end
C
C
C
C
C
C
C >>>>>> Extract a comment field.               
C
         function cmnt_(strg)          
C
         logical   cmnt_
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,
     *     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   ipp,itpos,ixpos
C
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         strg=' '
         long_=0
         cmnt_=.false.
         goto 105
100      jchar=jchar+1
105      if(jchar.le.lastch)     goto 140
C
C....... Read a new line
C
110      call getlin(flag)
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=4
           cmnt_=.false.
           return
         endif
         jchar=1
         strg=char(0)
         long_=1
         posnam_=0
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
C
C....... Process this character in the line
C
150      c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         goto 300
C
C        For a tab, when not expanding to blanks, accept
C        that single character as a comment
C
190      long_=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
C
C....... Accept the remainder of the line as a comment 
C
200      long_=lastch-jchar
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
210      posnam_=itpos
         if(long_.gt.0) then
           if(tabx_) then
             call detab
             ixpos= lastnb(bufntb)
             strg = bufntb(itpos+1:ixpos)
           else
             strg = buffer(jchar+1:lastch)
           endif
         endif
         if(long_.le.0) then
           strg=' '
           long_=1
         endif
         jchar=MAXBUF+1
220      lrecd=jjlrec
         cmnt_=.true.
         return
C
C....... Found a non-comment field, restore pointers
C
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer=' '
         if(lastch.gt.0)buffer(1:lastch)=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
         end
C
C
C
C
C
C >>>>> Convert name string to lower case
C        
         function locase(name)
C
         include     'ciftbx.sys'
         character    locase*(MAXBUF)
         character    temp*(MAXBUF),name*(*)
         character    low*26,cap*26,c*1
         integer i,j
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
C
         temp=name
         do 100 i=1,MAXBUF
         c=temp(i:i)
         if(c.eq.' ') goto 200
         if(c.eq.tab) goto 200
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)        
100      continue
200      locase=temp
         return
         end
C
C
C
C
C
C >>>>>> Get the data item associated with the tag.
C
         subroutine getitm(name)
C
         include   'ciftbx.sys'
         SAVE
         character name*(*)
         character flag*4
         integer   iitem,nitem,npakt
         integer   kchar,loopi,i,j,itpos,ipp
C
C....... Find requested dataname in hash list
C
         nametb=name
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         quote_=' '
         if(name(1:1).eq.'_')       goto 100
         type_='null'
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         tagname_=' '
         strg_=' '
         long_=1
         goto 1000
100      call hash_find(nametb,
     *     dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *     iname)
         if(iname.gt.0)             goto 180
         if(dictfl.ne.'yes') then
         call hash_find(nametb,
     *     dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,j)
         if(j.ne.0) then
           dictype_=dicxtyp(j)
           if(dcindex(j).ne.0) diccat_=dcname(dcindex(j))
           dicname_=nametb
           if(aroot(j).ne.0) then
             dicname_=dictag(aroot(j))
             call hash_find(dicnam(aroot(j)),
     *         dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *         iname)
             if(iname.gt.0)      goto 180
           endif
           type_='null'
           tagname_=' '
           strg_=' '
           long_=1
           go to 1000
         endif
         endif
160      continue
         type_='null'
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         long_=1
         goto 1000
C
C
180      tagname_=dtag(iname)
         if(ddict(iname).ne.0) tagname_=dictag(ddict(iname))
         posnam_=tchar(iname)
         if(tabx_)posnam_=xchar(iname)
         if(nloop(iname).le.0)      goto 500
C
C....... Process loop packet if first item request
C
         if(nloop(iname).ne.loopnl) goto 200
         if(lloop(iname).lt.loopct) goto 300
         if(loop_)                  goto 230
200      loop_=.true.
         loopct=0
         loopnl=nloop(iname)
         nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=drecd(iname)-1
         call getlin(flag)
         jchar=max(0,dchar(iname)-1)
Cdbg     if(jchar.lt.0) write(6,'(7H dchar ,i5)') jchar
         do 220 i=1,nitem
220      lloop(i+iname-iloop(iname))=0
         goto 240
C
C....... Read a packet of loop items
C
230      nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=looprd(nitem+1)-1
         call getlin(flag)
         jchar=loopch(nitem+1)
Cdbg     if(jchar.lt.0) write(6,'(7H loopch,i5)') jchar
240      iitem=0
250      iitem=iitem+1
         if(iitem.le.nitem)     goto 255
         loopch(iitem)=jchar
         looprd(iitem)=irecd
         goto 270
255      call getstr
         loopch(iitem)=jchar-long_
         if(quote_.ne.' ')loopch(iitem)=jchar-long_-1
         loopln(iitem)=long_
         looprd(iitem)=irecd
         if(buffer(1:1).ne.';')     goto 250
260      call getlin(flag)
         if(buffer(1:1).ne.';') goto 260
         jchar=2
         goto 250
270      loopct=loopct+1
         if(loopct.lt.npakt)    goto 300
         loop_=.false.
C
C....... Point to the loop data item
C
300      lloop(iname)=lloop(iname)+1
         loopi=iloop(iname)
         irecd=looprd(loopi)-1
         call getlin(flag)
         long_=loopln(loopi)
         kchar=loopch(loopi)
         goto 550
C
C....... Point to the non-loop data item
C
500      irecd=drecd(iname)-1
         call getlin(flag)
         kchar=dchar(iname)+1
         long_=iloop(iname)
         loop_=.false.
         loopct=0
         loopnl=0
C
C....... Place data item into variable string and make number
C
550      type_=dtype(iname)
         dictype_=dxtyp(iname)
         diccat_='(none)'
         if(cindex(iname).gt.0) diccat_=dcname(cindex(iname))
         if(cindex(iname).lt.0) diccat_=cname(-cindex(iname))
         if(diccat_.eq.' ') diccat_='(none)'
         dicname_=dtag(iname)
         if(ddict(iname).ne.0) then
           if (aroot(ddict(iname)).ne.0) then
             dicname_=dictag(aroot(ddict(iname)))
           endif
         endif
         strg_(1:long_)=buffer(kchar:kchar+long_-1)
         itpos=kchar
         if(tabx_) then
         itpos=0
         do ipp=1,kchar
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         endif
         posval_=itpos
         posend_=itpos+long_-1
         jchar=kchar+long_
         if(jchar.le.MAXBUF) then
           if(buffer(jchar:jchar).ne.' ' .and.
     *       buffer(jchar:jchar).ne.tab) jchar=jchar+1
         endif
         quote_=' '
         if(kchar.gt.1) then
           if(buffer(kchar-1:kchar-1).ne.' ' .and.
     *        buffer(kchar-1:kchar-1).ne.tab) then
             quote_=buffer(kchar-1:kchar-1)
           endif
         endif
         if(type_.eq.'char' .and. kchar.eq.1 .and.
     *     buffer(1:1).eq.';') type_='text'
         if(type_.eq.'text') then
           if(buffer(1:1).eq.';') then
             strg_(1:1)=' '
           else
             type_='char'
           endif
         endif
         if(type_.eq.'numb') then
           call ctonum
           if(posdec_.gt.0) posdec_=posval_+posdec_-1
         endif
         if(quote_.ne.' ') goto 1000
         if(long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
         if(long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
C
1000     return
         end
C
C
C
C
C
C
C
C >>>>>> Read the next string from the file
C
         subroutine getstr
C
         include   'ciftbx.sys'
         integer   i,j,jj(11),im
         logical   quoted
         character c*1,num*21,flag*4
         data num/'0123456789+-.()EDQedq'/
C
         quoted=.false.
         quote_=' '
         if(irecd.gt.0.and.
     *     jchar.eq.1.and.lastch.gt.0) goto 140
100      jchar=jchar+1
         if(jchar.le.lastch)     goto 150
C
C....... Read a new line
C
110      call getlin(flag)
         type_='fini'
         dictype_=type_
         diccat_='(none)'
         dicname_=' '
Cdbg     write(6,'(/5i5,a)') 
Cdbg *              irecd,jrecd,lrecd,nrecd,lastch, buffer(1:lastch)
         if(flag.eq.'fini')  goto 500
C
C....... Test if the new line is the start of a text sequence
C
140      if(buffer(1:1).ne.';') goto 150
         type_='text'
         jchar=lastch+1
         long_=lastch
         strg_(1:long_)=buffer(1:long_)
         strg_(1:1)=' '
         goto 500
C
C....... Process this character in the line
C
150      c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 110
         if(c.eq.'''')      goto 300
         if(c.eq.'"')       goto 300
         if(c.ne.'_')       goto 200
         type_='name'
         goto 210
C
C....... Span blank delimited token; test if a number or a character 
C
200      type_='numb'
         im=0
         do 205 i=1,11
205      jj(i)=0
210      do 250 i=jchar,lastch
         if(buffer(i:i).eq.' ')       goto 400
         if(buffer(i:i).eq.tab)       goto 400
         if(type_.ne.'numb')          goto 250
         j=index(num,buffer(i:i))
         if(j.eq.0)                 type_='char'
         if(j.le.10) then
           im=im+1
           goto 250
         endif
         if(j.gt.13.and.im.eq.0) type_='char'
         jj(j-10)=jj(j-10)+1
250      continue
         i=lastch+1
         if(type_.ne.'numb') goto 400
         do 270 j=1,5
         if((jj(j).gt.1.and.j.gt.2) .or.
     *     jj(j).gt.2)             type_='char'
270      continue
         goto 400
C
C....... Span quote delimited token; assume character
C
300      type_='char'
         quoted=.true.
         jchar=jchar+1
         do 320 i=jchar,lastch
         if(buffer(i:i).ne.c)             goto 320
         if(i+1.ge.lastch)                goto 400
         if(buffer(i+1:i+1).eq.' ')       goto 400
         if(buffer(i+1:i+1).eq.tab)       goto 400
320      continue
Cdbg     write(6,'(a,4i5,a)') 
Cdbg *       '**** ',irecd,lastch,i,jchar,buffer(jchar:i)       
         call warn(' Quoted string not closed')
C
C....... Store the string for the getter
C
400      long_=i-jchar
         strg_(1:long_)=buffer(jchar:i-1)
         jchar=i
         quote_=' '
         if(quoted) then
           quote_=buffer(jchar:jchar)
           jchar =jchar+1
         endif
Cdbg     write(6,'(5x,8i5,5x,a)') 
Cdbg *   irecd,jrecd,lrecd,nrecd,lastch,i,jchar,long_,strg_(1:long_)
         if(type_.ne.'char'.or.quoted) goto 500
         if(strg_(1:5).eq.'data_') type_='data'
         if(strg_(1:5).eq.'loop_') type_='loop'
         if(long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
         if(long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
         if(strg_(1:5).eq.'save_') type_='save'
C
500      return
         end
C
C
C
C
C
C
C >>>>>> Convert a character string into a number and its esd
C
C                                          Q
C                                          D+
C                                          E-
C                                +         +
C           number string        -xxxx.xxxx-xxx(x)
C           component count CCNT 11111222223333444
C           (with at least 1 digit in the mantissa)
C
         subroutine ctonum
C
         integer   lastnb
         include  'ciftbx.sys'
         character test*22,c*1
         integer*4 m,nchar
         integer*4 ccnt,expn,msin,esin,ndec,ids,nmd
         double precision numb,sdev,ntemp,mant
         data test /'0123456789+.-()EDQedq '/
C
         numbtb=0.D0
         sdevtb=-1.D0
         numb=1.D0
         sdev=0.D0
         ccnt=0
         mant=0.D0
         expn=0.
         msin=+1
         esin=+1
         ndec=0
         ids=0
         nmd=0
         type_='char'
         posdec_=0
         if(long_.eq.1.and.
     *     index('0123456789',strg_(1:1)).eq.0) goto 500
C
C....... Loop over the string and identify components
C
C        The scan works in phases
C          ccnt = 0   processing looking for first digit
C          ccnt = 1   processing before decimal point
C          ccnt = 2   processing after decimal point
C          ccnt = 3   processing exponent
C          ccnt = 4   processing standard deviation
C
         do 400 nchar=1,long_
C
         c=strg_(nchar:nchar)
         m=index(test,c)
         if(m.eq.0)     goto 500
         if(m.gt.10)    goto 300
C
C....... Process the digits
C
         if(ccnt.eq.0)  ccnt=1
         if(ccnt.eq.2)  ndec=ndec+1
         if(ccnt.gt.2)  goto 220
         ntemp=m-1
         mant=mant*10.D0+ntemp
         nmd=nmd+1
         if(ccnt.eq.1.and.mant.ne.0.D0) ids=ids+1
         goto 400
220      if(ccnt.gt.3)  goto 240
         expn=expn*10+m-1
         goto 400
240      ntemp=m-1
         sdev=sdev*10.D0+ntemp
         sdevtb=1.D0
         goto 400
C
C....... Process the characters    . + - ( ) E D Q
C
300      if(c.ne.'.')  goto 320
         if(ccnt.gt.1) goto 500
         posdec_=nchar
         ccnt=2
         goto 400
C
320      if(nmd.eq.0.and.m.gt.13) goto 500
         if(c.ne.'(')  goto 340
         if(posdec_.eq.0) posdec_=nchar
         ccnt=4
         goto 400
C
340      if(posdec_.eq.0.and.ccnt.gt.0) posdec_=nchar
         if(c.eq.' ')  goto 400
         if(m.gt.13) m = 11
         if(ccnt.eq.3) goto 500
         if(ccnt.gt.0) goto 360
         ccnt=1
         msin=12-m
         goto 400
360      ccnt=3
         esin=12-m
C
400      continue
         if(posdec_.eq.0) posdec_=lastnb(strg_(1:long_))+1
C
C....... String parsed; construct the numbers
C
         expn=expn*esin-ndec
         if(expn+ids.gt.-minexp) then
           call warn(' Exponent overflow in numeric input')
           expn=-minexp-ids
         endif
         if(expn.lt.minexp) then
           call warn(' Exponent underflow in numeric input')
           expn=minexp
         endif
         if(expn.lt.0) numb=1./10.D0**abs(expn)
         if(expn.gt.0) numb=10.D0**expn
         if(sdevtb.gt.0.0) sdevtb=numb*sdev
         ntemp=msin
         numbtb=numb*mant*ntemp
         type_='numb'
C
500      return
         end
C
C
C
C
C
C
C >>>>>> Read a new line from the direct access file
C
         subroutine getlin(flag)
C
         integer   lastnb
         include  'ciftbx.sys'
         character flag*4
C
         irecd=irecd+1
         jchar=1
         if(irecd.eq.jrecd)  goto 200
         if(irecd.le.lrecd)  goto 100
         buffer=' '
         lastch=0
         jchar=MAXBUF+1
         jrecd=-1
         flag='fini'
         goto 200
100      read(dirdev,'(a)',rec=irecd) buffer
         recn_=irecd
         lastch=max(1,lastnb(buffer))
         jrecd=irecd
         flag=' '
200      return
         end
C
C
C
C
C
C
C >>>>>> Detab buffer into bufntb
C
         subroutine detab
C
         include   'ciftbx.sys'
         integer   icpos,itpos,ixpos,ixtpos
         if(jrecd.eq.jrect) return
         icpos=1
         itpos=1
         bufntb=' '
         if(lastch.gt.0) then
100      ixpos=index(buffer(icpos:lastch),tab)
         ixtpos=ixpos+itpos-1
         if(ixpos.gt.0.and.ixtpos.le.MAXBUF) then
           ixtpos=((ixtpos+7)/8)*8
           if(ixpos.gt.1) then
           bufntb(itpos:ixtpos)=
     *       buffer(icpos:ixpos+icpos-2)
           else
           bufntb(itpos:ixtpos)=' '
           endif
           itpos=ixtpos+1
           icpos=ixpos+icpos
           goto 100
         else
           bufntb(itpos:max(MAXBUF,itpos+lastch-icpos))=
     *       buffer(icpos:lastch)
         endif
         endif
         jrect=jrecd
         return
         end
C
C
C
C
C
C
C >>>>>> Write error message and exit.
C
         subroutine err(mess)
         character*(*) mess
         call cifmsg('error',mess)
         stop
         end
C
C
C
C
C
C
C >>>>>> Write warning message and continue.
C
         subroutine warn(mess)
         character*(*) mess
         call cifmsg('warning',mess)
         return
         end
C
C
C
C
C
C
C >>>>>> Write a message to the error device
C
         subroutine cifmsg(flag,mess)
C
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) flag
         character*(*) mess
         character*(MAXBUF)  tline
         character*5   btype
         integer       ll,ls,ltry,ii,i
C
         btype = 'data_'
         if(save_) btype = 'save_'
         tline= ' ciftbx '//flag//': '
     *   //file_(1:longf_)//' '//btype
     *   //bloc_(1:max(1,lastnb(bloc_)))//' line:'
         ll = max(1,lastnb(tline))
         write(errdev,'(a,i7)')tline(1:ll),irecd
         ll=len(mess)
         ls=1
100      if(ll-ls.le.79) then
           write(errdev,'(1X,a)') mess(ls:ll)
           return
         else
           ltry = min(ll,ls+79)
           do ii = ls+1,ltry
           i = ltry-ii+ls+1
           if(mess(i:i).eq.' ') then
             write(errdev,'(1X,a)') mess(ls:i-1)
             ls=i+1
             if(ls.le.ll) go to 100
             return
           endif
           enddo
           write(errdev,'(1X,a)') mess(ls:ltry)
           ls=ltry+1
           if(ls.le.ll) go to 100
           return
         endif  
         end
C
C
C
C
C >>>>>> Create a named file.
C
         function pfile_(fname)
C
         logical   pfile_
         include   'ciftbx.sys'
         logical   test
         integer   i
         character fname*(*)
C
C....... Test if a file by this name is already open.
C
         if(pfilef.eq.'yes') call close_
         pfilef='no '
         file_=fname
         do 120 i=1,MAXBUF
         if(file_(i:i).eq.' ') goto 140
120      continue
140      if (i.gt.1) then
           inquire(file=file_(1:i-1),exist=test)
           pfile_=.false.
           longf_ = i-1
           if(test)            goto 200
         else
           file_ = ' '
           pfile_ = .true.
           longf_ = 1
         endif
C
C....... Open up a new CIF
C
         if (file_(1:1) .ne. ' ')  then 
         open(unit=outdev,file=fname,status='NEW',access='SEQUENTIAL',
     *                    form='FORMATTED')
         precn_=0
         endif
         pfile_=.true.  
         pfilef='yes'
         nbloc=0
         pchar=1+lprefx
         pcharl=0
         obuf=prefx
         obuf(pchar:MAXBUF)=' '
200      return
         end
C
C
C
C
C
C >>>>>> Store a data block command in the CIF
C        Call with blank name to close current block only
C
         function pdata_(name) 
C
         logical   pdata_
         include  'ciftbx.sys'
         character name*(*),temp*(MAXBUF)
         character dbloc(100)*(NUMCHAR)
         integer   i
C
         pdata_=.true.
         if(ploopn.ne.0)     call eoloop
         if(ptextf.eq.'yes') call eotext
         if(psaveo) then
           pchar=-1
           if(pposval_.ne.0) then
             pchar=lprefx+1
             call putstr(' ')
             pchar=lprefx+pposval_
             pposval_=0
           endif
           call putstr('save_')
           psaveo=.false.
         endif
C
C....... Check for duplicate data name
C
         temp=name
         if(temp.eq.' ')        goto 200
         if(saveo_)             goto 130
         pdata_=.false.
         do 120 i=1,nbloc
         if(temp.eq.dbloc(i))   goto 200
120      continue
C
C....... Save block name and put data_ statement
C
         nbloc=nbloc+1
         if(nbloc.le.100) dbloc(nbloc)=temp
130      pchar=-1
         temp='data_'//name
         if(saveo_) temp='save_'//name
         psaveo=saveo_
         if(pposnam_.gt.0) then
           pchar=lprefx+1
           call putstr(' ')
           pchar=lprefx+pposnam_
           pposnam_=0
         endif
         call putstr(temp)         
         pchar=lprefx
         pdata_=.true.
C
200      return
         end
C
C
C
C
C
C
C >>>>>> Put a number into the CIF, perhaps with an esd appended
C
         function pnumb_(name,numb,sdev)
C
         logical    pnumb_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         real       numb,sdev
         double precision dnumb,dsdev,dprec
C
         pnumb_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call eotext
C
         if(name(1:1).eq.' ')   goto 120
         if(vcheck.eq.'no ')    goto 100
         call dcheck(temp,'numb',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0)
     *       temp=dictag(aroot(xdchk))
         endif
         pnumb_=flag
100      if(ploopn.ne.0)        call eoloop
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         call putstr(temp)
C
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         dprec=decprc
         dnumb=numb
         dsdev=sdev
         call putnum(dnumb,dsdev,dprec)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
C
150      pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         return
         end
C
C
C
C
C
C
C >>>>>> Put a double precision number into the CIF, perhaps 
C        with an esd appended
C
         function pnumd_(name,numb,sdev)
C
         logical    pnumd_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         double precision numb,sdev
C
         pnumd_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call eotext
C
         if(name(1:1).eq.' ')   goto 120
         if(vcheck.eq.'no ')    goto 100
         call dcheck(temp,'numb',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0)
     *       temp=dictag(aroot(xdchk))
         endif
         pnumd_=flag
100      if(ploopn.ne.0)        call eoloop
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         call putstr(temp)
C
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         call putnum(numb,sdev,dpprc)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
C
150      pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         return
         end
C
C
C
C
C
C
C >>>>>> Put a character string into the CIF.
C
         function pchar_(name,string)      
C
         logical    pchar_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR),string*(*)
         character  line*(MAXBUF),strg*(MAXBUF)
         integer    i,j
C
         pchar_=.true.
         flag  =.true.
         tflag =.true.
         temp  =name
         if(ptextf.eq.'yes') call eotext
C
         if(name(1:1).eq.' ')   goto 110
         if(vcheck.eq.'no ')    goto 100
         call dcheck(temp,'char',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0)
     *       temp=dictag(aroot(xdchk))
         endif
         pchar_=flag
100      if(ploopn.ne.0)        call eoloop
         pchar=-1
         if(pposnam_.gt.0) pchar=posnam_+lprefx
         call putstr(temp)
C
110      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         line=string
         do 120 i=MAXBUF,2,-1
         if(line(i:i).ne.' ') goto 130
120      continue
130      if(pposval_.ne.0.and.pposend_.ge.pposval_)
     *      i=max(i,pposend_-pposval_+1)
         if(pquote_.ne.' ')   goto 150
         do 140 j=i,1,-1
         if(line(j:j).eq.' ') goto 150
140      continue
         if((line(1:1).eq.'_' 
     *     .or. line(i:i).eq.'_'
     *     .or. line(1:1).eq.''''
     *     .or. line(1:1).eq.'"'
     *     .or. line(1:1).eq.';')
     *     .and.line(1:i).ne.'''.'''
     *     .and.line(1:i).ne.'''?'''
     *     .and.line(1:i).ne.'"."'
     *     .and.line(1:i).ne.'"?"') goto 150
         strg=line(1:i)
         goto 200
150      if(pquote_.eq.';')       goto 190
         if(pquote_.eq.'''')      goto 165
         if(pquote_.eq.'"')       goto 185
         do 160 j=1,i
         if(line(j:j).eq.'''')    goto 170
160      continue
165      strg=''''//line(1:i)//''''
         i=i+2
         pquote_=''''
         goto 200
170      do 180 j=1,i
         if(line(j:j).eq.'"')     goto 190
180      continue
185      strg='"'//line(1:i)//'"'
         i=i+2
         pquote_='"'
         goto 200
190      pchar=-1
         strg='; '//line(1:i)
         i=i+2
         ptextf='yes'
         call putstr(strg(1:i))
         pchar=-1
         ptextf='no '
         call putstr(';')
         pchar=lprefx
         call putstr(' ')
         call warn(' Converted pchar_ output to text for: '
     *     //strg(3:i))
         goto 210
C
200      if(pposval_.ne.0) then
           pchar=pposval_+lprefx
           if(pquote_.ne.' ') pchar=pchar-1
         endif
         call putstr(strg(1:i))   
210      if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if((.not.tflag).and.line(1:i).ne.'.'.and.
     *     line(1:i).ne.'?'.and.pquote_.eq.' ') then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
250      pposval_=0
         pposdec_=0
         pposnam_=0
         pposend_=0
         pquote_=' '
         return
         end
C
C
C
C
C
C >>>>>> Put a comment in the output CIF 
C
         function pcmnt_(string)     
C
         logical    pcmnt_
         include   'ciftbx.sys'
         character  string*(*), temp*(MAXBUF)
C
         if(ptextf.eq.'yes') call eotext
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if(string.eq.' '.or.
     *     (string.eq.char(0)) .or.
     *     (string.eq.tab.and.(.not.ptabx_))) then
           if(string.eq.' ') pchar=-1
           call putstr(string)
           if(string.eq.' ') call putstr(char(0))
         else
           temp='#'//string
           call putstr(temp)
           call putstr(char(0))
         endif
         pcmnt_=.true.
         pposnam_=0
         if(string.ne.tab)pchar=lprefx+1
         return
         end
C
C
C
C
C
C
C
C >>>>>> Put a text sequence into the CIF.
C
         function ptext_(name,string)      
C
         logical    ptext_
         integer    lastnb
         include   'ciftbx.sys'
         logical    flag,tflag
         integer    ll
         character  name*(*),temp*(NUMCHAR),string*(*),store*(NUMCHAR)
         character  temp2*(MAXBUF)
         data store/'                                '/
C
         ptext_=.true.
         flag  =.true.
         tflag =.true.
         ll=lastnb(string)
         if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         temp=name
         if(ptextf.eq.'no ')    goto 100
         if(temp.eq.store)      goto 150
         call eotext
C
100      if(name(1:1).ne.' ')   goto 110
         if(ptextf.eq.'yes')    goto 150
         goto 130
C
110      if(ploopn.ne.0)        call eoloop
         if(vcheck.eq.'no ')    goto 120
         call dcheck(name,'char',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0)
     *       temp=dictag(aroot(xdchk))
         endif
         ptext_=flag
120      pchar=-1
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         call putstr(temp)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
130      ptextf='yes'
         store=temp
         if(string(1:1).eq.' '.and.ll.gt.1) then
           pchar=-1
           temp2=';'//string(2:ll)
           call putstr(temp2)
           pchar=-1
           return
         endif
         pchar=-1
         call putstr(';')
         pchar=-1
         if(string.eq.' ') return
150      pchar=-1
         call putstr(string(1:max(1,ll)))
         pchar=-1
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         return
         end
C
C
C
C
C
C
C >>>>>> Put a loop_ data name into the CIF.
C
         function ploop_(name)      
C
         logical    ploop_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
C
         ploop_=.true.
         flag  =.true.
         if(ptextf.eq.'yes')    call eotext
         if(ploopf.eq.'no ')    call eoloop
         temp=' '
         if(name(1:1).eq.' ')   goto 100
C
         if(tabl_.and.pposnam_.eq.0) then
           temp='    '//name
         else
           temp=name
         endif
         if(vcheck.eq.'no ')    goto 100
         call dcheck(name,'    ',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0) then
             if(tabl_.and.pposnam_.eq.0) then
               temp='    '//dictag(aroot(xdchk))
             else
               temp=dictag(aroot(xdchk))
             endif
           endif
         endif
         ploop_=flag
100      if(ploopn.ne.0)        goto 120
         ploopf='yes'
         pchar=-1
         if(pposval_.ne.0) then
           pchar=lprefx+1
           call putstr(' ')
           pchar=pposval_+lprefx
         else
           if(pposnam_.ne.0) then
             pchar=lprefx+1
             call putstr(' ')
             pchar=pposnam_+lprefx+1
           endif
         endif
         call putstr('loop_')
         pchar=-1
         if(name(1:1).eq.' ') then
           ploopn=-1
           return
         endif
120      if(pposnam_.ne.0) pchar=pposnam_+lprefx
         call putstr(temp)
         if(flag)               goto 130
         if(.not.tabl_) pchar=lprefx+57
         call putstr('#< not in dictionary')
130      pchar=lprefx+1
         ploopn=max(ploopn,0)+1
C
150      return
         end
C
C
C
C
C
C >>>>>> Create or clear a prefix string
C        Any change in the length of the prefix string flushes
C        pending text, if any,  loops and partial output lines
C
         function prefx_(strg,lstrg)
C
         logical    prefx_
         include   'ciftbx.sys'
         character  strg*(*)
         integer    lstrg,mxline
C
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         if(lstrg.ne.lprefx.and.pcharl.gt.0) then
           pchar=-1
           call putstr(' ')
         endif
         if (lstrg.le.0) then
           prefx=' '
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx
           lprefx=0
         else
           if(lstrg.gt.mxline) then
             call warn(' Prefix string truncated')
           endif
           prefx=strg
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx+lstrg
           obuf(1:min(mxline,lstrg))=prefx
           lprefx=lstrg
           if(mxline-lprefx.lt.NUMCHAR) then
             call warn(' Output prefix may force line overflow')
           endif
         endif
         prefx_=.true.
         return
         end
C
C
C
C
C
C
C >>>>>> Close the CIF
C
         subroutine close_
C
         include   'ciftbx.sys'
C
         if(ptextf.eq.'yes') call eotext
         if(ploopn.ne.0)     call eoloop
         if(pcharl.ge.lprefx+1) then
           pchar=-1
           call putstr(' ')
         endif
         if (file_(1:1) .ne. ' ') then
           close(outdev)
           precn_=0
         endif
         return
         end
C
C
C
C
C
C >>>>>> Put the string into the output CIF buffer 
C
         subroutine putstr(string)     
C
         integer    lastnb
         include   'ciftbx.sys'
         SAVE
         character  string*(*),temp*(MAXBUF),bfill*(MAXBUF)
         character  temp2*(MAXBUF)
         integer    i,ii,mxline,ioffst,ifree,icpos,itpos
         integer    ixpos,ixtpos,it,im,kbin,kpass
         logical    pflush,waslop
         data       waslop /.false./

C
         bfill = ' '
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         temp=string
         temp2=temp
         pflush=.false.
         if(pchar.lt.0) pflush=.true.
C
         do 100 i=MAXBUF,1,-1
         if(temp(i:i).eq.' ')              goto 100
         if(ptabx_.and.temp(i:i).eq.tab)    goto 100
         goto 110
100      continue
         i=0
         it=i
C
C....... Organise the output of loop_ items
C
110      if(i.eq.0)             goto 130
         if(i.eq.1.and.string.eq.tab) goto 130
         if(i.eq.1.and.string.eq.char(0)) then
           pcharl=MAXBUF
           goto 200
         endif           
         if(temp(1:1).eq.'#')   goto 130
         if(ploopf.eq.'yes')    goto 130
         if(ptextf.eq.'yes')    goto 130
         if(ploopn.le.0)        goto 130
         ploopc=ploopc+1
         if(align_.or.tabl_) then
           if(ploopc.gt.ploopn) then
             if(pcharl.gt.lprefx) pflush=.true.
             ploopc=1
             if(pchar.gt.0) pchar=lprefx+1
           endif
           if(pchar.lt.0)    goto 130
           if(tabl_) then
           kbin=(mxline-lprefx)/8
           if(ploopn.lt.kbin) then
             if(kbin/(ploopn+1).gt.1) then
             pchar=9+lprefx+
     *         (ploopc-1)*8*(kbin/(ploopn+1))
             else
             pchar=1+lprefx+
     *         (ploopc-1)*8*(kbin/ploopn)
             endif
           else
             if(ploopc.le.kbin) then
               pchar=1+lprefx+(ploopc-1)*8
             else
               kpass=(ploopc-kbin-1)/(kbin-1)+1
               pchar=2*kpass+1+lprefx+
     *           mod(ploopc-kbin-1,kbin-1)*8
             endif
           endif
           else
             if(ptabx_) then
             icpos=1
             itpos=1
120          ixpos=index(temp(icpos:i),tab)
             ixtpos=(pchar+itpos-1+ixpos)
             ixtpos=((ixtpos+7)/8)*8
             if(ixpos.gt.0) then
               if(ixpos.gt.1) then
                 temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
               else
                 temp2(itpos:ixtpos-pchar+1)=' '
               endif
               icpos=ixpos+1
               itpos=ixtpos+2-pchar
               if(icpos.le.i) goto 120
               it=itpos-1
             else
               temp2(itpos:itpos+i-icpos)=temp(icpos:i)
               it=itpos+i-icpos
             endif
             endif
             if((pchar+i).gt.mxline+1.or.
     *          (ptabx_.and.pchar+it.gt.mxline+1)) then
               if(pcharl.gt.lprefx)pflush=.true.
               pchar=lprefx+1
             endif
           endif
         else
           if(ploopc.le.ploopn)   goto 130
           ploopc=1
         endif
C
C....... Is the buffer full and needs flushing?
C
130      if(i.eq.1.and.string.eq.tab) then
           if(pcharl.gt.lprefx) then
             if(obuf(pcharl:pcharl).eq.' ') pcharl=pcharl-1
           endif
         endif
         if(pchar.le.pcharl.and.pcharl.gt.lprefx) pflush=.true.
         pchar=max(lprefx+1,pchar)
         if((ploopf.eq.'yes'.or.ploopn.le.0).and.tabl_)
     *     pchar=((pchar-lprefx+6)/8)*8+1+lprefx
         if(ptabx_) then
         icpos=1
         itpos=1
135      ixpos=index(temp(icpos:i),tab)
         ixtpos=(pchar+itpos-1+ixpos)
         ixtpos=((ixtpos+7)/8)*8
         if(ixpos.gt.0) then
           if(ixpos.gt.1) then
             temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
           else
             temp2(itpos:ixtpos-pchar+1)=' '
           endif
           icpos=ixpos+1
           itpos=ixtpos+2-pchar
           if(icpos.le.i) goto 135
           it=itpos-1
         else
           temp2(itpos:itpos+i-icpos)=temp(icpos:i)
           it=itpos+i-icpos
         endif
         endif
         if((pchar+i).gt.mxline+1.or.
     *     (ptabx_.and.pchar+it.gt.mxline+1)) then
            pflush=.true.
            pchar=mxline+1-i
            pchar=max(lprefx+1,pchar)
         endif
         if(.not.pflush)  goto 150
140      if(pcharl.gt.lprefx) then
           if(waslop.or.(.not.tabl_)) goto 145
           ioffst=0
           pcharl=max(lastnb(obuf(1:pcharl)),lprefx+1)
           ifree=mxline-pcharl
           if(ifree.gt.0) then
           im=numtab+2
           if(numtab.gt.0.and.numtab.le.MAXTAB) then
             if(obuf(itabp(numtab):itabp(numtab)).eq.'#')
     *         im=im-1
           endif
           if(ifree.ge.16.and.im.lt.4.and.
     *       (obuf(1+lprefx:1+lprefx).ne.'#'
     *        .and.obuf(1+lprefx:1+lprefx).ne.';'
     *        .and.obuf(1+lprefx:1+lprefx).ne.'_'
     *        .and.obuf(1+lprefx:1+lprefx).ne.' '
     *        .and.obuf(1+lprefx:5+lprefx).ne.'data_'
     *        .and.obuf(1+lprefx:5+lprefx).ne.'save_'
     *        .and.obuf(1+lprefx:5).ne.'loop_')) then
             temp(1+lprefx:pcharl)=obuf(1+lprefx:pcharl)
             obuf(1+lprefx:pcharl+8)=
     *         bfill(1:8)//temp(1+lprefx:pcharl)
             ioffst = 8
             ifree=ifree-8
             pcharl=pcharl+8
           endif
           do ii=1,min(MAXTAB,numtab)
             icpos=itabp(ii)+ioffst
             if(icpos.gt.pcharl)   goto 145
             if(im.lt.4) then
             itpos=(max(icpos-lprefx,
     *         ii*(mxline-lprefx)/im)+6)/8
             itpos=itpos*8+1+lprefx
             else
             itpos=(max(icpos-lprefx,
     *         ii*(mxline-lprefx)/im)+4)/6
             itpos=itpos*6+1+lprefx
             endif
             if((obuf(icpos:icpos).eq.''''.or.
     *          obuf(icpos:icpos).eq.'"').and.
     *          itpos.gt.icpos) itpos=itpos-1
             if(itpos-icpos.gt.ifree) itpos=icpos+ifree
             if(itpos.gt.icpos) then
               temp(1:pcharl-icpos+1)=
     *           obuf(icpos:pcharl)
               if(i.lt.numtab) then
                 ixpos=itabp(ii+1)+ioffst
                 if(ixpos.gt.icpos+itpos-icpos+1) then
                   if(obuf(ixpos-(itpos-icpos+1):ixpos-1).eq.
     *               bfill(1:itpos-icpos+1)) then
                     temp(ixpos-itpos+1:pcharl-itpos+1)=
     *               obuf(ixpos:pcharl)
                     pcharl=pcharl-(itpos-icpos)
                   endif
                 endif
               endif
               obuf(icpos:pcharl+itpos-icpos)=
     *           bfill(1:itpos-icpos)//temp(1:pcharl-icpos+1)
               ifree=ifree-(itpos-icpos)
               ioffst=ioffst+itpos-icpos
               pcharl=pcharl+itpos-icpos
             endif
             if(ifree.le.0)      goto 145
           enddo
           endif
145        pcharl=max(1,lastnb(obuf))
           write(outdev,'(a)') obuf(1:pcharl)
         else
           if(precn_.gt.0) then
           if(lprefx.gt.0) then
           write(outdev,'(a)') obuf(1:lprefx)
           else
           write(outdev,'(a)')
           endif
           else
           precn_=precn_-1
           endif
         endif
         waslop=.false.
         precn_=precn_+1
         do ii = 1,MAXTAB
           itabp(ii)=0
         enddo
         numtab=0
         if(lprefx.gt.0) then
           obuf=prefx(1:lprefx)
         else
           obuf=' '
         endif
C
C....... Load the next item into the buffer
C
150      pcharl=pchar+i
         if(ptabx_) pcharl=pchar+it
         waslop= ploopf.eq.'no '.and.ploopn.gt.0.and.align_
         if(i.eq.0) then
           if(pcharl.eq.lprefx+1.and.
     *       obuf(lprefx+1:lprefx+1).eq.' ') pcharl=pcharl-1
             pchar=pcharl+1
           goto 200
         endif
         if(ptabx_) then
           obuf(pchar:pcharl)=temp2(1:it)
         else
           if(string.eq.tab) pcharl=pcharl-1
           obuf(pchar:pcharl)=string(1:i)
         endif
         if(pchar.gt.1+lprefx) then
           numtab=numtab+1
           if(numtab.le.MAXTAB) itabp(numtab)=pchar
         endif
         pchar=pcharl+1
         if(pchar.gt.mxline+2) then
           call warn(' Output CIF line longer than line_')
         endif
C
200      return
         end
C
C
C
C
C
C >>>>>> Convert the number and esd to string nnnn(m), limited
C        by relative precision prec
C
         subroutine putnum(numb,sdev,prec)  
C
         include   'ciftbx.sys'
         character  string*30,temp*30,c*1,sfmt*8
         double precision numb,sdev,prec,xxnumb,xsdev,slog
         integer    i,iexp,ifp,ii,jj,j,jlnz,jn,kexp,m,ixsdev,islog
         integer    kdecp,ibexp
C
         kdecp=0
         if (sdev.gt.abs(numb)*prec) then
           if (esdlim_.ne.esdcac) then
C
C            determine the number of digits set by esdlim_
C
             if (esdlim_.lt.9 .or.esdlim_.gt.99999) then
               call warn(' Invalid value of esdlim_ reset to 19')
               esdlim_ = 19
             endif
C
C            determine the number of esd digits
C
             esddig = 1.+alog10(float(esdlim_))
             esdcac = esdlim_
           endif
C
C          determine kexp, the power of 10 necessary
C          to present sdev as an integer in the range
C          (esdlim_/10,esdlim_]
C
           slog = dlog10(sdev)
           islog = slog+1000.
           islog = islog-1000
           kexp = -islog+esddig
C
C          Adjust exponent kexp, so that sdev*10**kexp
C          is in the interval (esdlim_/10,esdlim_]
C
 20        if (kexp.lt.minexp) then
             call warn(' Underflow of esd')
             ixsdev = 0
             go to 30
           endif
           if (kexp.gt.-minexp) then
             call warn(' Overflow of esd')
             ixsdev = 99999
             go to 30
           endif
           xsdev = sdev*10.D0**kexp
           ixsdev = xsdev+.5
           if (ixsdev.gt.esdlim_) then
             kexp = kexp -1
             go to 20
           endif
           if (ixsdev.lt.(esdlim_+5)/10) then
             kexp = kexp+1
             go to 20
           endif
C
C          We need to present the number to the same scaling
C          at first, but will adjust to avoid Ennn notation
C          if possible
C
 30        xxnumb = dabs(numb)*10.d0**kexp+.5
           if(xxnumb*prec .gt.1.D0) then
             call warn(' ESD less than precision of machine')
             ixsdev=0
           endif
           if(numb.lt.0.d0) xxnumb = -xxnumb
           write(string,ndpfmt)xxnumb
C
C          Extract the power of 10
C
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = string(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 40
               endif
             endif
           enddo
           call err(' Internal error in putnum')
C
C          Scan the rest of the string shifting the
C          decimal point to get an integer
C
40         ifp = 0
           j=1
           do ii = 1,i-1
           c = string(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               temp(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.temp(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 50
               endif
             else
               if(c.eq.'.') then
                 ifp=1
                 if(iexp.le.0) goto 50
               endif
             endif
           endif
           enddo
C
C          The string from 1 to j-1 has an integer
C          If iexp < 0, we present a 0.  If iexp > 0
C          we pad with zeros
C
50         if(j.eq.1.or.iexp.lt.0) then
             temp(1:1)='0'
             j=2
             iexp = 0
           endif
           if (iexp.gt.0) then
             do ii = 1,iexp
             temp(j:j)='0'
             j=j+1
             enddo
             iexp=0
           endif
           string=temp(1:j-1)
C
C          We have the number for which the presentation
C          would be nnnnnE-kexp.  If kexp is gt 0, we can
C          decrease it and introduce a decimal point
C
           jj=0
           if(index('0123456789',temp(1:1)).eq.0) jj=1
           if(kexp.gt.0.and.kexp.lt.j-jj+8) then
             if(kexp.lt.j-1) then
               string=temp(1:j-1-kexp)//'.'//
     *         temp(j-kexp:j-1)
               kexp = 0
               j=j+1
             else
               if(jj.ne.0)string(1:1)=temp(1:1)
               string(1+jj:1+jj)='.'
               do ii=1,kexp-(j-1-jj)
                 string(1+jj+ii:1+jj+ii)='0'
               enddo
               string(2+jj+(kexp-(j-1-jj)):30)=
     *           temp(1+jj:j-1)
               j=j+1+kexp-(j-1-jj)
               kexp=0
             endif
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.eq.0) kdecp=j
           if(kexp.ne.0) then
             write(temp(1:5),'(i5)') -kexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
C
C          if there is a standard deviation
C          append it in parentheses
C
           if(ixsdev.ne.0) then
             write(temp(1:5),'(i5)') ixsdev
             string(j:j)='('
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
             string(j:j)=')'
             j=j+1
           endif
         else
C
C          There is no standard deviation, just write numb
C          But limit to the digits implied by prec
C
           slog = dlog10(min(.1D0,max(prec,dpprc)))
           islog = slog+1000.5
           islog = islog-1000
           kexp = -islog
           write(sfmt,'(5h(D30.,i2,1h))') kexp
           write(temp,sfmt)numb
C
C          Now have the number in the form 
C          [sign][0].nnnnnnnnDeee
C          which, while sufficient, is not neat
C          we reformat for the case 0<=eee<=kexp
C
C
C          Extract the power of 10
C
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = temp(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 140
               endif
             endif
           enddo
           call err(' Internal error in putnum')
C
C          Scan the rest of the string shifting the
C          decimal point to get a number with exponent 0,
C          if possible
C
140        ifp = 0
           j=1
           do ii = 1,i-1
           jn=ii
           c = temp(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               string(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.string(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 150
               endif
             else
               if(c.eq.'.') then
                 ifp = -1
                 if(iexp.le.0) goto 150
               endif
             endif
           endif
           enddo
150        string(j:j)='.'
           ifp = j
           j = j+1
           jlnz = j-1
155        do ii = jn+1,i-1
             c = temp(ii:ii)
             if (c.ne.' ')then
               m=index('0123456789',c)
               if(m.ne.0) then
                 string(j:j)=c
                 j=j+1
                 if(m.ne.1)jlnz=j
                 if(m.eq.1.and.ifp.ge.1.and.
     *             pposdec_.ne.0.and.pposend_.ne.0) then
                   if(j-1-ifp-min(iexp,0).le.pposend_-pposdec_)
     *               jlnz=j
                 endif
               else
                 goto 160
               endif
             endif
           enddo
160        j=jlnz
           if(j.eq.1) then
            string(1:1)='0'
            j=2
           endif
           if(iexp.lt.0.and.iexp.gt.-7.and.ifp.lt.j-1.and.
     *       ifp.ne.0.and.j-ifp-iexp.le.kexp) then
             temp(1:ifp)=string(1:ifp)
             do ii = 1,-iexp
               temp(ifp+ii:ifp+ii) = '0'
             enddo
             temp(ifp-iexp+1:j-iexp-1) = string(ifp+1:j-1)
             j = j-iexp
             iexp=0
             string(1:j-1) = temp(1:j-1)
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.eq.0) kdecp=j
           if(iexp.ne.0) then
             write(temp(1:5),'(i5)')iexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
         endif
C
         if(j.lt.1) then
           string(1:1)='0'
           j=2
         endif
         if(kdecp.lt.1)kdecp=j
         if(pposdec_.ne.0) then
           pchar=lprefx+pposdec_-kdecp+1
         else
           if(pposval_.ne.0)pchar=lprefx+pposval_
         endif
         call putstr(string(1:j-1))
         return
         end
C
C
C
C
C
C >>>>>> Check dictionary for data name validation    
C
         subroutine dcheck(name,type,flag,tflag)
C
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR),
     *              locase*(MAXBUF),type*4
C
         flag=.true.
         tflag=.true.
         temp=locase(name)
         call hash_find(temp,
     *     dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,xdchk)
         if(xdchk.eq.0) goto 150
         if(tcheck.eq.'no ')          goto 200
         if(type.eq.dictyp(xdchk))    goto 200
         if(type.eq.'    ')           goto 200
         if(dictyp(xdchk).eq.'text' .and. type.eq.'char') goto 200
         tflag=.false.
         goto 200
150      flag=.false.
200      continue
         return
         end
C
C
C
C
C
C >>>>>> End of text string
C
         subroutine eotext
C
         include   'ciftbx.sys'
C
         if(ptextf.ne.'yes') then
           call warn(' Out-of-sequence call to end text block')
           return
         endif
         ptextf='no '
         pchar=-1
         call putstr(';')
         call putstr(char(0))
         return
         end
C
C
C
C
C
C >>>>>> End of loop detected; check integrity and tidy up pointers
C
         subroutine eoloop
C
         include   'ciftbx.sys'
         integer   i
C
         if(ploopn.eq.0)          goto 200
         if(ploopn.eq.-1) then
           call putstr('_DUMMY')
           ploopn=1
           ploopc=0
           call warn(
     *       ' Missing: missing loop_ name set as _DUMMY')
         endif
         if(ploopn.eq.ploopc)     goto 200
         do 150 i=ploopc+1,ploopn
150      call putstr('DUMMY')
         call warn(    
     *         ' Missing: missing loop_ items set as DUMMY')
C
200      ploopc=0
         ploopn=0
         return
         end
C
C
C
C
C
C
C >>>>>> Set common default values
C
         block data
C
         include   'ciftbx.sys'
         data cifdev     /1/
         data outdev     /2/
         data dirdev     /3/
         data errdev     /6/
         data loopct     /0/
         data nhash      /0/
         data ndict      /0/
         data nname      /0/
         data nbloc      /0/
         data ploopn     /0/
         data ploopc     /0/
         data ploopf     /'no '/
         data ptextf     /'no '/
         data pfilef     /'no '/
         data testfl     /'no '/
         data vcheck     /'no '/
         data tcheck     /'no '/
         data align_     /.true./
         data tabl_      /.true./
         data tabx_      /.true./
         data ptabx_     /.true./ 
         data text_      /.false./
         data loop_      /.false./
         data ndcname    /0/
         data ncname     /0/
         data save_      /.false./
         data saveo_     /.false./
         data alias_     /.true./
         data aliaso_    /.false./
         data dchash     /NUMHASH*0/
         data dichash    /NUMHASH*0/
         data dhash      /NUMHASH*0/
         data dcchain    /NUMDICT*0/
         data aroot      /NUMDICT*0/
         data cindex     /NUMBLOCK*0/
         data line_      /80/
         data lastch     /0/
         data dictype_   /' '/
         data dicname_   /' '/
         data diccat_    /' '/
         data tagname_   /' '/
         data prefx      /' '/
         data lprefx     /0/
         data esdlim_    /19/
         data esdcac     /19/
         data esddig     /2/
         data esdfmt     /'(e12.2)'/
         data edpfmt     /'(d12.2)'/
         data ndpfmt     /'(d30.14)'/
         data decprc     /1.e-6/
         data dpprc      /1.d-14/
         data decmin     /1.e-37/
         data dpmin      /1.d-307/
         data minexp     /-307/
         data itabp      /MAXTAB*0/
         data jrect      /-1/
         data numtab     /0/
         data recn_      /0/
         data precn_     /0/
         data posnam_    /0/
         data posval_    /0/
         data posdec_    /0/
         data posend_    /0/
         data pposnam_   /0/
         data pposval_   /0/
         data pposdec_   /0/
         data pposend_   /0/
         data quote_     /' '/
         data pquote_    /' '/
         data ibkmrk     /MAXBOOK*-1,MAXBOOK*-1,
     *                    MAXBOOK*-1,MAXBOOK*-1/

         end
C
C
C       change the following include to include 'clearfp_sun.f'
C       for use on a SUN
C
        include 'clearfp.f'
