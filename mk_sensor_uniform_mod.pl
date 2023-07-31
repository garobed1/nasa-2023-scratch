#!/usr/bin/perl

# $Id: mk_sensor_uniform.pl,v 1.4 2023/06/17 03:22:29 mnemec Exp $

# interpolates points along a line to a uniform spacing

# M. Nemec, Feb 2020, NASA ARC

use strict;
use warnings;
use Carp qw{ croak };
use English qw{ -no_match_vars };
use FileHandle;
use Getopt::Long;
use File::Basename;

local $PROGRAM_NAME = basename($PROGRAM_NAME);
my $FAILED = "$PROGRAM_NAME FAILED";

my $help;
my $verbose;
my $lsf; # signature file
my $cmdline = $PROGRAM_NAME;
my $trimh =  9999999; # default big, not active
my $triml = -9999999; # default big, not active
my $np;

my $usage = <<"END_USAGE";

Usage: $PROGRAM_NAME
    -v      verbose
    -f      signature file name (lineSensor, sBOOM, Tecplot)
    -n      number of uniformly spaced points in output file

Option to trim the output data:
    -triml  trim low side of signature to define metric section, e.g. avoid
            upstream region, default <-9999999>
    -trimh  trim high side of signature to define metric section, e.g. avoid
            trailing wake, default <9999999>
END_USAGE

if (@ARGV) {
  foreach (@ARGV) {
    $cmdline .= /\s/ ? " \'" . $_ . "\'" : " " . $_;
  }
  GetOptions(
    'f=s'    => \$lsf,
    'n=i'    => \$np,
    'triml=f'=> \$triml,
    'trimh=f'=> \$trimh,
    'v'      => \$verbose,
    'h'      => \$help,
    'help'   => \$help,
    ''       => \$help
  );
}
die "$usage\n" if $help;

my $DIV0_EPSILON = 1.e-12;

my $file_type = 0;

print "#\n# Command: $cmdline\n#\n" if $verbose;

my ($r_x,$r_f,$ox,$oy,$oz,$dx,$dy,$dz) = read_data_file($lsf);

$np = $#$r_x if ( ! $np );

my $ds = @{ $r_x }[$#$r_x] / ( $np-1 ); # desired uniform mesh spacing

print "# Uniform spacing $ds\n" if $verbose;

my @uniform_x;
my $length = 0.0;
my $i = 0;
while ( $i != $np ) {
  $i++;
  push @uniform_x, $length;
  $length += $ds;
}

# interpolate lineSensor onto uniform_x
my @uniform_f = linear_interp(\@uniform_x, $r_x, $r_f);

# output

my $outfile;

if ( $lsf =~ /lineSensor_/ ) {

  # output Cart3D lineSensor file

  if ( $lsf =~ m/ee_lineSensor_/ ) {
    ($outfile = $lsf) =~ s/ee_lineSensor_/uniform_ee_lineSensor_/;
  }
  else {
    ($outfile = $lsf) =~ s/lineSensor_/uniform_lineSensor_/;
  }

  print "# Writing $outfile\n" if $verbose;

  my $ofh = FileHandle->new( $outfile, '>' );
  croak "ERROR: Cannot write\n" unless ( defined $ofh );

  print $ofh "# $cmdline\n";
  print $ofh "# ", scalar localtime;
  print $ofh "\n#\n";

  print $ofh "# Verbatum comments from source lineSensor file\n";

  my $ifh = FileHandle->new( $lsf, '<' );
  croak "ERROR: Cannot read\n" unless ( defined $ofh );

  while (<$ifh>) {
    next unless /\w/; # skip unless line contains alphanumeric char

    if ( /^#/ ) {
      print $ofh "$_";
    }
    else {
      last;
    }
  }
  undef $ifh;

  my $i = 0;
  foreach my $loc (@uniform_x) {
    my $x = $ox + $loc/$uniform_x[$#uniform_x]*($dx-$ox);
    my $y = $oy + $loc/$uniform_x[$#uniform_x]*($dy-$oy);
    my $z = $oz + $loc/$uniform_x[$#uniform_x]*($dz-$oz);

    $ofh->printf("%16.9e %16.9e %16.9e %20.13e",$x,$y,$z,$loc);
    foreach my $val ( @{ $uniform_f[$i] } ) {
      $ofh->printf(" %16.9e",$val);
    }
    $i++;
    print $ofh "\n";
  }
  undef $ofh;
}
else {
  # output sBOOM ground signature
  if ( $lsf =~ m/modSig/ ) {
    ($outfile = $lsf) =~ s/modSig/uniform_modSig/;
  }
  else {
    ($outfile = $lsf) =~ s/SBground/uniform_SBground/;
  }
  

  print "# Writing $outfile\n" if $verbose;

  my $ofh = FileHandle->new( $outfile, '>' );
  croak "ERROR: Cannot write\n" unless ( defined $ofh );

  print $ofh "# $cmdline\n";
  print $ofh "# ", scalar localtime;
  print $ofh "\n#\n";

  print $ofh "## Verbatum comments from source file\n";

  my $ifh = FileHandle->new( $lsf, '<' );
  croak "ERROR: Cannot read\n" unless ( defined $ofh );

  while (<$ifh>) {
    next unless /\w/; # skip unless line contains alphanumeric char

    if ( /^[[:alpha:]]/ ) {
      print $ofh "# $_";
    }
    elsif ( /^#/ ) {
      print $ofh "$_";
    }
    else {
      last;
    }
  }
  undef $ifh;
  print $ofh "##\n";

  my $i = 0;
  foreach my $loc (@uniform_x) {
    $ofh->printf("%20.13e",$loc);
    foreach my $val ( @{ $uniform_f[$i] } ) {
      $val = 0.0 if ( abs($val) < 1.e-12 );
      $ofh->printf(" %20.13e",$val);
    }
    $i++;
    print $ofh "\n";
  }

  undef $ofh;
}

print "#\n" if $verbose;

exit 0; # all done

sub linear_interp {

  # we interpolate f from locations at x to fi at desired locations xi

  # xi ... desired location
  # x  ... data location
  # f  ... data values

  my ($r_xi, $r_x, $r_f ) = @_;

  my @fi; # interpolated data, return value

  my $il = 0;
  my $iu = $#$r_f;
  my $ip = int(( $iu -$il )/2);

  my ($j, $df, $st, $sn);

  for ($j = 0; $j <= $#$r_xi; $j++) {

    # outside current bracket?
    if ( ( 1 == ( $iu - $il ) ) && ( $$r_xi[$j] > $$r_x[$iu] ) ) {
      $il = $iu;
      $iu = $#$r_x;
      $ip = int( ($iu - $il)/2) + $il;
    }

    # find bracket via binary search
    while ( ( $iu - $il ) > 1 ) {
      if ( $$r_xi[$j] < $$r_x[$ip] ) {
        $iu = $ip;
      }
      else {
        $il = $ip;
      }
      $ip = int(($iu - $il)/2) + $il;
    }

    if ( $iu < $il ) {
      print("ERROR: Target sensor interpolation: upper and lower bound swapped\n");
      exit 1;
    }

    # boundary check: if point on sensor is less than location of first
    # target value, use direct injection
    if ( $$r_xi[$j] < $$r_x[$il] ) {
      $df = 1.;
    }
    # boundary check: if point on sensor is greater than location of last
    # target value, use direct injection
    elsif ( $$r_xi[$j] > $$r_x[$iu] ) {
      $df = 0.;
    }
    # usual case when we have a valid bracket
    else {
      $st = $$r_x[$iu] - $$r_x[$il];
      if ( $st < $DIV0_EPSILON ) {
        print("Target sensor interpolation: coincident points detected\n");
        $df = 0.5;
      }
      else {
        $sn = $$r_xi[$j] - $$r_x[$il];
        $df = ( $st - $sn ) / $st; # distance fraction
      }
    }

    if ( $df < 0. || $df > 1. ) {
      printf("Target sensor interpolation: unable to interpolate target value\n");
      exit 1;
    }

    # r_f contains references to all variables of line sensor

    my $f_il = $$r_f[$il];
    my $f_iu = $$r_f[$iu];

    my $i;
    my @vals;

    # loop over 'i' interpolates all variables of line sensor

    for ($i=0; $i <= $#$f_il; $i++) {
      push @vals, $df*$f_il->[$i] + (1.-$df)*$f_iu->[$i];
    }

    $fi[$j] = \@vals; # array of references to arrays
  }

  return @fi;
}

sub read_data_file {
  my $f  = shift;

  my (@rf, @rx);

  print "# Reading $f\n" if $verbose;

  my $ifh = FileHandle->new( "$f", '<' );
  croak "ERROR: Cannot read $f\n" unless ( defined $ifh );

  while (<$ifh>) {
    next unless /\w/; # skip unless line contains alphanumeric char
    if ( /ZONE\s+T/ || /zone t/ ) {
      # Tecplot
      $file_type = 1;
      last;
    }
    elsif ( /^# Line Sensor/ ) {
      # Cart3D lineSensor
      $file_type = 2;
      last;
    }
    elsif ( /Number of ground signatures/ ) {
      # sBOOM ground signature file
      $file_type = 3;
      last;
    }
    elsif ( /Number of Waveforms/ ) {
      # sBOOM near signature file, treat as ground sig?
      $file_type = 5;
      last;
    }
  }

  # back to beginning
  seek $ifh, 0, 0;

  # if still unknown file type, see if it looks like x,y file
  if ( 0 == $file_type ) {
    while (<$ifh>) {
      next if /^#/;     # skip comments
      next unless /\w/; # skip unless line contains alphanumeric char
      if ( /[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?/ ) {
        $file_type = 4;
        last;
      }
    }
    seek $ifh, 0, 0;
  }

  my ($ox,$oy,$oz);
  my ($dx,$dy,$dz);

  if ( 2 == $file_type ) {
    # lineSensor file
    while (<$ifh>) {
      next unless /\w/; # skip unless line contains alphanumeric char
      chomp;

      if ( /^# Line Origin / ) {
        my @line = split ' ';
        $oz = pop @line;
        $oy = pop @line;
        $ox = pop @line;
      }
      elsif ( /^# Line Destination / ) {
        my @line = split ' ';
        $dz = pop @line;
        $dy = pop @line;
        $dx = pop @line;
      }

      next if /^#/;     # skip comments

      my @line = split ' ';

      shift @line; # skip x
      shift @line; # skip y
      shift @line; # skip z

      my $loc = shift @line; # save distance
      if ( $loc >= $triml && $loc <= $trimh ) {
        # trim signature
        push @rx, $loc;        # save distance
        push @rf, \@line;      # save all variables
      }
    }

    # prepend start of lineSensor
    if ( 0.0 < $rx[0] ) {
      unshift @rx, '0.0';
      unshift @rf, $rf[0];
    }

    # postpend end of lineSensor
    my $s = sqrt( ($dx - $ox)*($dx - $ox) + ($dy - $oy)*($dy - $oy) + ($dz - $oz)*($dz - $oz) );

    if ( $s > $rx[$#rx] ) {
      push @rx, $s;
      push @rf, $rf[$#rf];
    }
  }
  elsif ( 3 == $file_type or 5 == $file_type ) {
    # sBOOM ground signature file
    while (<$ifh>) {
      next unless /^\s*[-+]?[0-9]/;  # skip unless it is a digit
      chomp;

      my @line = split ' ';

      my $loc = shift @line; # save time
      if ( $loc > $triml && $loc < $trimh ) {
        # trim signature
        push @rx, $loc;
        push @rf, \@line;   # save overpressure
      }
    }

    $oz = 0.0;
    $oy = 0.0;
    $ox = 0.0;

    if ( $triml > -9999999 ) {
      unshift @rx, $triml;
      unshift @rf, $rf[0];
    }

    if ( $trimh < 9999999 ) {
      $dx = $trimh;
      push @rx, $trimh;
      push @rf, $rf[$#rf];
    }
    else {
      $dx = $rx[$#rx];
    }
    $dy = 0.0;
    $dz = 0.0;

    if ( $triml > -9999999 ) {
      my $i;
      for ($i=0; $i<=$#rx; $i++) {
        $rx[$i] -= $triml;
      }
    }
  }
  else {
    print "ERROR: unknown data file $f\n";
    exit 1;
  }
  undef $ifh;

  return (\@rx, \@rf, $ox, $oy, $oz, $dx, $dy, $dz);
}
