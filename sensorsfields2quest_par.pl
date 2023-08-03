#!/usr/bin/env perl

# $Id: sensorsfields2quest.pl,v 1.4 2022/11/28 21:55:08 mnemec Exp $

# Format Cart3D lineSensors and sBOOM ground signatures as QUEST xml files

# Also format KL-generated atmospheric fields as QUEST xml files

# M. Nemec, G. Bedonian NASA Ames, Feb 2023

use strict;
use warnings;
use English qw{ -no_match_vars };
use FileHandle;
use File::Basename;
use File::Spec;
use File::Find;
use Getopt::Long;
use Cwd;
use XML::LibXML;

my $sec;
my $help;
my $cmdline = $PROGRAM_NAME;

if (@ARGV) {
  foreach (@ARGV) {
    $cmdline .= /\s/ ? " \'" . $_ . "\'" : " " . $_;
  }
  GetOptions(
    'p=s'   => \$sec,
    'h'     => \$help,
    'help'  => \$help,
    ""      => \$help
  );
}
die "Usage: $PROGRAM_NAME [-v]
Parses loud.dat and output.out files to generate QUEST xml data files
    -p      p0000th index to run\n"
    if ($help);



my $dbdir;
if ( -d 'cases' ) {
  $dbdir = 'cases';
}
else {
  exit 0;
}

my @caseDirs;
opendir(my $dh, $dbdir) || die "Can't open directory: $!";
while (readdir $dh) {
  my $dir = File::Spec->catdir($dbdir, $_);
  next unless (-d $dir);
  next unless ( $dir =~ /case\.\d\d\.${sec}\d+/ );
  push @caseDirs, $dir;
}
closedir $dh;

exit 0 unless (@caseDirs);

my @lineSensors;
my %sensorNames;

find sub { push @lineSensors, $File::Find::name, if (-f and
/^uniform_lineSensor_/ and ( $File::Find::dir =~ /adapt/ and $File::Find::dir =~ /FLOW/
) ) }, @caseDirs;

find sub { push @lineSensors, $File::Find::name, if (-f and
/^uniform_SBground/ and ( $File::Find::dir =~ /sboom/ ) ) }, @caseDirs;

find sub { push @lineSensors, $File::Find::name, if (-f and
/^uniform_modSig/ and ( $File::Find::dir =~ /sboom/ ) ) }, @caseDirs;

find sub { push @lineSensors, $File::Find::name, if (-f and
/^TEMP_profile/)}, @caseDirs;

find sub { push @lineSensors, $File::Find::name, if (-f and
/^HUMIDITY_profile/)}, @caseDirs;

if ( ! @lineSensors ) {
  find sub { push @lineSensors, $File::Find::name, if (-f and /^uniform_lineSensor_/ ) }, @caseDirs;
}

exit 0 unless (@lineSensors);

foreach my $sensor ( @lineSensors ) {
  my $file        = basename($sensor);
  my $sensor_path = dirname($sensor);

  my $name;
  ( $name = $file ) =~ s/uniform_//;
  $name =~ s/\.dat|\.sig|\.txt//;

  if (not exists( $sensorNames{$name} ) ) {
    my $loudnessdir = "QUEST_${name}";
    mkdir $loudnessdir;
    $sensorNames{$name} = $loudnessdir;
    print "# Created directory $loudnessdir\n";
  }

  my $ee_sensor = $sensor;
  $ee_sensor =~ s/_lineSensor_/_ee_lineSensor_/;

  print "Processing $sensor_path\n";

  my @dirs = File::Spec->splitdir( $sensor_path );

  my $case;
  while (@dirs ) {
    $case = shift @dirs;
    last if ( $case =~ /\./ );
  }

  my @scratch    = split '\.', $case;

  my $questfile  = "${name}.00000." . $scratch[1] . '.' . $scratch[2] . '.xml';

  my (@abscissas, @ordinates);

  if ( $name =~ /^lineSensor/ ) {
    @abscissas = ( 'X', 'Y', 'Z', 'Distance Along Sensor' );
    @ordinates = ( 'dP/Pinf', 'RHO', 'U', 'V', 'W', 'P', 'T', 'SoundSpeed' );
  }
  elsif ( $name =~ /^SBground/ ) {
    @abscissas = ( 'Time (ms)' );
    @ordinates = ( 'Overpressure (psf)' );
  }
  elsif ( $name =~ /^modSig/ ) {
    @abscissas = ( 'Time (ms)' );
    @ordinates = ( 'Overpressure (psf)' );
  }
  elsif ( $name =~ /^TEMP/ ) {
    @abscissas = ( 'Altitude (1000 ft)' );
    @ordinates = ( 'Temperature (F)' );
  }
  elsif ( $name =~ /^HUMIDITY/ ) {
    @abscissas = ( 'Altitude (1000 ft)' );
    @ordinates = ( 'Relative Humidity (%)' );
  }
  else {
    print "ERROR: unknown file type\n";
    exit 1;
  }

  my @re;
  if (-e "$ee_sensor") {
    my $efh = FileHandle->new( "$ee_sensor", '<' );
    die "ERROR: Cannot read $ee_sensor, stopped " unless ( defined $efh );

    while (<$efh>) {
      next if /^#/;      # skip comments
      next unless /\w/;  # skip empty lines
      chomp;
      ( my @scratch ) = split ' ', $_;
      push @re, $scratch[4];
#      push @re, $scratch[5]; # no smoothing
    }

    undef $efh;
  }

  my $fh = FileHandle->new( "$sensor", '<' );
  die "ERROR: Cannot read $sensor, stopped " unless ( defined $fh );

  # output QUEST xml data file
  my $doc = XML::LibXML::Document->createDocument( "1.0", "ISO-8859-1" );

  my $element = $doc->createElement('root');
  $doc->setDocumentElement($element);

  my $child = $doc->createElement("${name}");
  $element = $element->addChild($child);

  $child = $doc->createElement("data");
  $element = $element->addChild($child);

  $child = $doc->createElement("abscissas");
  $child->appendText(scalar @abscissas);
  $element->addChild($child);

  my $i = 0;
  foreach ( @abscissas ) {
    $child = $doc->createElement("abscissa${i}_name");
    $child->appendText("$_");
    $element->addChild($child);
    $i++;
  }

  $child = $doc->createElement("ordinates");
  $child->appendText(scalar @ordinates);
  $element->addChild($child);

  $i = 0;
  foreach (@ordinates) {
    $child = $doc->createElement("ordinate${i}_name");
    $child->appendText("$_");
    $element->addChild($child);
    $child = $doc->createElement("ordinate${i}_type");
    if ( /Pinf/ and @re ) {
      $child->appendText("1");
    }
    else {
      $child->appendText("0");
    }
    $element->addChild($child);
    $i++;
  }

  $i = 0;
  while (<$fh>) {
    next if /^#/;      # skip comments
    next unless /\w/;  # skip empty lines
    chomp;
    ( my @scratch ) = split ' ', $_;

    if ( $i > 0 ) {
      $child = $doc->createElement("data");
      $element = $element->addChild($child);
    }

    $child = $doc->createElement("index");
    $child->appendText("$i");
    $element->addChild($child);

    foreach (0...$#abscissas) {
      $child = $doc->createElement("abscissa${_}_val");
      $child->appendText(shift @scratch);
      $element->addChild($child);
    }

    foreach (0...$#ordinates) {
      $child = $doc->createElement("ordinate${_}_val");
      $child->appendText(shift @scratch);
      $element->addChild($child);
      if ( ($ordinates[$_] =~ /Pinf/) and @re ) {
        $child = $doc->createElement("ordinate${_}_error");
        $child->appendText(shift @re);
        $element->addChild($child);
      }
    }

    $element = $element->parentNode;
    $i++;
  }

  undef $fh;

  my $doc_string = $doc->toString(2); # pretty format
  $doc_string =~ s/>\n\s*\n/>\n/g;    # remove excessive new-line sequences
  undef $doc;

  $fh = FileHandle->new( File::Spec->catfile($sensorNames{$name},
                                             $questfile), '>' );
  print $fh $doc_string;
  undef $fh;
}

# create jar bundle for faster reading
for my $path ( split /:/, $ENV{PATH} ) {
  if ( -e "${path}/jar"
       && -f "${path}/jar"
       && -x "${path}/jar" )
  {
    print "Creating jar bundle using '$path/jar cf', inspect with '$path/jar tf'\n";

    foreach my $name (keys %sensorNames) {
      chdir $sensorNames{$name};
      # the command apparently gets too long if done with wildcards
      system("jar cMf ../$name.bundle.xml .");# $name.*.??.*.xml");
      print "Done file $name.bundle.xml\n";
      unlink glob "$name.*.??.*.xml";
      chdir File::Spec->updir();
    }
    last;
  }
}

exit 0;
