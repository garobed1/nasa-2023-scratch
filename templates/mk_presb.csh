#!/bin/csh -f

# $Id: mk_presb.csh,v 1.1 2020/07/23 17:57:34 mnemec Exp $

# generate sBOOM input file by getting parameters from various sources

# MN, July 2020

set file = 'presb.input'
set src  = '.'
set phi  = 0

set height = '54000.0'
set rf     = '1.9'
set ground = '360.9285'

set mach = '1.4'

set RoLs = 1
set Lref = 1
set dist = '270.0'

# set ftmfac = 0.3048006
set lineSensor = 'lineSensor_LINE1.dat'

while ( $#argv )
  set phi = $argv[1]
  shift
end # while

if ( -e ../boom.inputs.json ) then
  set RoLs = ( `grep '^\"RoLs\"' ../boom.inputs.json | sed 's/^[^[]*\[//' | sed 's/\]//' | sed 's/,/ /g'` )
  set Lref = `grep '^\"Lref\"'   ../boom.inputs.json | sed 's/^[^:]*://'  | sed 's/,//'`
  set lineSensor = 'lineSensor_RoL'$RoLs[1]'.a.phi'$phi'.dat'
  set dist = `echo "$Lref*$RoLs[1]" | bc -l`
endif

if ( -e ../input.cntl ) then
  set mach = `grep '^Mach' ../input.cntl | awk '{print $2}'`
endif

if ( -e ../CFD/input.cntl ) then
  set mach = `grep '^Mach' ../CFD/input.cntl | awk '{print $2}'`
endif

if ( -e ../QUEST.dat ) then
  grep '^altitude' ../QUEST.dat >& /dev/null
  if ( 0 == $status ) then
    set height = `grep '^altitude' ../QUEST.dat | awk '{print $2}'`
  endif

  grep '^reflection_factor' ../QUEST.dat >& /dev/null
  if ( 0 == $status ) then
    set rf = `grep '^reflection_factor' ../QUEST.dat | awk '{print $2}'`
  endif

  grep '^ground_height' ../QUEST.dat >& /dev/null
  if ( 0 == $status ) then
    set ground = `grep '^ground_height' ../QUEST.dat | awk '{print $2}'`
  endif
endif

# begin file

cat << PARAMS >! $file
sig.dat         1. Name of wave file, from $lineSensor
$mach           2. Mach number
$height         3. Cruise Altitude (ft)
$dist           4. Starting distance of the wave from the aircraft axis (ft) Lref = $Lref
20000           5. Number of points to be used in the waveform during propagation
0               6. Number of zero-padding points - these number of points are padded in front
1               7. Non-linear flag  - yes (1) or no (0)
1               8. Thermoviscous flag  - yes (1) or no (0)
1               9. Relaxation flag  - yes (1) or no (0)
1.0             10.Initial step size - this is dynamically controlled within sboom
$rf             11.Reflection factor at ground level
500             12.Number of points needed in the resampled signature (!!no longer supported!!)
8e-6            13.slope tolerance for outputting the ground signature without zero paddings
-1              14.input tolerance
0               15.Propagate Pressure waveform(0), F-function(1) or Equivalent Area(2)
$ground         16.Height of the ground in feet
1               17.Number of Azimuths
$phi               18.azimuthal angle for off-track analysis
PARAMS

#if ( 0 == 1 ) then
if ( -e "$src/../TEMP_profile.txt") then # NOTE: We're checking just for existence of a specifically named user profile
  echo '1               19. User Input Temperature Flag: 0 => Std. Atmosphere, 1 => user input' >> $file
  # NOTE: ATT Data already given in Fahrenheit, altitude given in 1000s of feet. Will convert accordingly
  cat $src/../TEMP_profile.txt | awk '{if ( $1*1000*0.3048006 < 90001 ) {printf("%8.5f %8.13f\n", $1*1000*0.3048006, $2)}}' >! tmp.txt 
  @ lines = `wc -l tmp.txt | awk '{print $1}'`
  printf "%-5d %9s If temperature flag is 1, then input the number of user altitudes (m) temperatures (F)\n" "$lines" >> $file
  cat $file tmp.txt > ${file}.tmp
  \mv -f ${file}.tmp $file
  \rm -f tmp.txt
else
  echo '0               19. User Input Temperature Flag: 0 => Std. Atmosphere, 1 => user input' >> $file
endif

if ( -e "$src/../WINDX_profile.txt" &&  -e "$src/../WINDY_profile.txt") then
  echo '1               20.User Input Wind Flag: 0 => No wind, 1 => user input wind profile' >> $file
  # NOTE: I don't know the units for wind velocity, I think it's feet/s? Need to change
  cat $src/../WINDX_profile.txt | awk '{if ( $1*1000*ftmfac < 90001 ) {printf("%8.5f %8.13f\n", $1*1000*0.3048006, $7)}}' >! tmp.txt
  @ lines = `wc -l tmp.txt | awk '{print $1}'`
  printf "%-5d %9s number of altitude (m) X-wind (m/s) pairs\n" "$lines" >> $file
  cat $file tmp.txt > ${file}.tmp
  \mv -f ${file}.tmp $file
  \rm -f tmp.txt
  printf "%-5d %9s number of altitude (m) Y-wind (m/s) pairs\n" "$lines" >> $file
  cat $src/../WINDY_profile.txt | awk '{if ( $1*1000*0.3048006 < 90001 ) {printf("%8.5f %8.13f\n", $1*1000*0.3048006, $7)}}' >! tmp.txt
else
  echo '0               20.User Input Wind Flag: 0 => No wind, 1 => user input wind profile' >> $file
endif

if ( -e "$src/../HUMIDITY_profile.txt") then # NOTE: We're checking just for existence of a specifically named user profile
  echo '1               21.User Input relative humidity flag: 0 => Std. Atmosphere, 1 => user input RH profile' >> $file
  cat $src/../HUMIDITY_profile.txt | awk '{if ( $1*1000*0.3048006 < 90001 ) {printf("%8.5f %8.13f\n", $1*1000*0.3048006, $2)}}' >! tmp.txt 
  @ lines = `wc -l tmp.txt | awk '{print $1}'`
  printf "%-5d %9s If humidity flag is 1, then input the number of user altitudes (m) relative humidity (percent)\n" "$lines" >> $file
  cat $file tmp.txt > ${file}.tmp
  \mv -f ${file}.tmp $file
  \rm -f tmp.txt
else
  echo '0               21.User Input relative humidity flag: 0 => Std. Atmosphere, 1 => user input RH profile' >> $file
endif


cat << REST >> $file
0.0             22.Heading angle in degrees. 180 heading means away from X axis, 90 heading means away from Y-axis
0.0             23.Climb angle in degrees
0               24.outflag
2               25.Input in Inches ( 0-feet, 1-inches, 2-meters)
0               26.Adjoint mode flag, if 1 adjoint sboom will be run along with sboom
1               27.Cost function Mode (Four choices 1/2/3/4 - See user manual)
0               28.rate of change of mach number (dM/dt in  - 1/sec)
0               29.turn rate (degrees/sec)
0               30.climb rate (degrees/sec)
0               31.3D ellipsoid (0= flat earth, 1=ellipsoid)
REST

exit 0
