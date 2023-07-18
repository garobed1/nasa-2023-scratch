#!/bin/csh -f

# $Id: quest2sboom.csh,v 1.1 2020/07/23 17:55:17 mnemec Exp $

# prepare sBOOM runs

# MN, July 2020 

set files = ( sig.dat c3d_close_dpp.py )

set source = 'templates'

if (! -e $source/mk_presb.csh ) then
  exit 0
endif

if (! -d cases ) then
  exit 0
endif

cd cases

foreach x (`\ls -1d case.??.?????`)
  cd $x

  \rm -rf sboom >& /dev/null

  mkdir sboom
  cd sboom

  # now in cases/case.00.00000/sboom
  
  foreach f ( $files )
    if ( -e ../../../$source/$f ) then
      \cp -f ../../../$source/$f .
    endif
  end

  # make sboom input file
  ../../../$source/mk_presb.csh
  
  cd ../..
end

# done with cases
cd ..

exit 0
