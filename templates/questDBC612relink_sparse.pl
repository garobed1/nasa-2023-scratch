#!/usr/bin/env perl

# $Id: questDBC612relink.pl,v 1.7 2023/09/29 22:06:42 mnemec Exp $

# Parses database file generated by QUEST
# Separates CFD cases from other disciplines to save computational cost

# Modified to link to CFD cases in the C612A-P1NC-uq directory, specifically
# This allows us to change the sBOOM prepconfig/database without even
# needing to rerun the CFD

# Downside: You are restricted to whatever UQ was performed for the CFD-only
# cases. We may be able to lower nested quadrature levels, though

# To summarize
# CFD-only database: C612A-P1NC-uq,  
# 1 param
# Mach: 1.4, 0.015, DCC L4 (ACTUAL RUNS DO NOT DELETE)
# This database: weathersim-C612A-uq
# Mach: 1.4, 0.015, DCC L4 or lower (links to runs in previous database)
# KLTEMPHUMID

# 2ND ALTERNATIVE
# It seems sparse databases behave differently, so here is a second
# alternative to questDB, this goes 0000111122223333 ... 9999 etc.

# NOTE: NEITHER OF THESE WILL WORK WITH MIXED DENSE/SPARSE WITHIN SBOOM

# CFD keywords: Mach, alpha, beta, gamma, InletPressRatioBC,
# InletVelocityBC, RoL, Phi

# M. Nemec, NASA ARC, Dec 2019, updated to QUEST V4, Sept 2022
# G. Bedonian, July 2023

use strict;
use warnings;
use English qw{ -no_match_vars };
use Getopt::Long;
use FileHandle;
use File::Basename;
use File::Spec;
use File::Path;

my $help;
my $verbose;
my $force;
my $ifile = 'database';
my $ofile = 'QUEST.dat';

my $usage = <<"END_USAGE";
Usage: $PROGRAM_NAME
    -i      QUESTPrep database <database>
    -f      Force overwrite of cases directory
    -v      Verbose

Script parses database file generated by QUEST ver. 4 and generates a
cases directory. If CFD parameters are present, a CFD database
directory is generated as well. CFD random parameters are Mach, alpha,
beta, gamma, InletPressRatioBC, InletVelocityBC, RoL and Phi.
END_USAGE

if (@ARGV) {
  GetOptions(
    'i=s'  => \$ifile,
    'f'    => \$force,
    'v'    => \$verbose,
    'h'    => \$help,
    'help' => \$help,
    ''     => \$help
  );
}
die "$usage\n" if $help;

local $PROGRAM_NAME = basename($PROGRAM_NAME);
my $FAILED = "$PROGRAM_NAME FAILED";

my $fh = FileHandle->new( "$ifile", '<' );
die "$FAILED Cannot read $ifile, stopped" unless ( defined $fh );

my ($n_independent_evals, $n_evals, $n_vars, $n_seq);
my ($quest_major_rev, $quest_minor_rev);
my $have_prepconfig;

my @var_names;
my (@cfd_params, @cfd_cases);
my (%cases, %cfd, %case2cfd);
my $n_cfd = -1;
my %pconfig;

if ( -d 'cases' and not $force ) {
  print "WARNING: 'cases' directory exists, use -f to overwrite\n";
  exit 0;
}
elsif ( -d 'cases' ) {
  rmtree('cases');
}

print " o Parsing QUEST $ifile\n" if $verbose;

my $line    = 0;
my $plength = 0; # used for STDOUT
while (<$fh>) {
  next unless /\w/; # skip unless line contains alphanumeric character
  $line++;
  next if /^#/;     # skip comments
  chomp;

  if ( 3 == $line ) {
    # line 3 has 4 fields, see QUEST manual
    ($n_independent_evals, $n_evals, $quest_major_rev, $quest_minor_rev) = split ' ';
  }
  elsif ( 4 == $line ) {
    # line 4 has names of random variables
    @var_names = split ' ';
    my $ivar=0;
    foreach (@var_names) {
      if ( flowParam($_) ) {
        push @cfd_params, $ivar;
      }
      $ivar++;
    }
  }
  elsif ( $line > 4 && $line <= (4+$n_evals)) {
    # list of required realizations
    my ($i, $j, $level, $seq, $nlev, @var_vals) = split ' ';

    my $case = sprintf("%02d.%05d",$level,$j);
    push @{$cases{$case}}, @var_vals;

    if ( @cfd_params) {
      my $run; # string representing cfd run
      my %pv;  # parameter values

      foreach (@cfd_params) {
        $run .= "_" if $run;
        $run .= "$var_names[$_]" . "$var_vals[$_]";
        $pv{"$var_names[$_]"} = $var_vals[$_];
      }

      if ( not exists $cfd{$run} ) {
        $cfd{$run} = \%pv;
        $n_cfd++;
        $cfd_cases[$n_cfd] = $run; # explicitly store hash keys so we can keep order
      }
      $case2cfd{$case} = $n_cfd;
    }
  }
  elsif ( $line == 4+$n_evals+1) {
    if (/^1/) {
      $have_prepconfig = 1;
    }
  }
  elsif ( $line > 4+$n_evals+1) {

    next if /==== begin prepconfig configuration/;
    last if /==== end prepconfig configuration/;
    last if /==== begin custom PDF/;
    last unless ( $have_prepconfig );

    # these lines define the independent parameters
    my ($param, $mean, $minval, $maxval, $sigma, $minsigma, $maxsigma,
        @scratch) = split ' ';

    my %v;

    # parname : m : ml : mu : σ : σl : σu : ω0 : ω0,l : ω0,u : ω1 :
    #           ω1,l : ω1,u : nlevels : n_eval : pdf : algorithm
    # parameter name string
    # mean or s1 pdf parameter
    # mean or s1 pdf parameter lower bound
    # mean or s1 pdf parameter upper bound
    # standard deviation or s2 or 1/2-width pdf parameter
    # standard deviation or s2 or 1/2-width pdf parameter lower bound
    #          standard deviation or s2 or 1/2-width pdf parameter upper bound
    #          shape parameter 0 (log-normal,beta pdfs)
    # shape parameter 0 lower bound (log-normal,beta pdfs)
    # shape parameter 0 upper bound (log-normal,beta pdfs)
    # shape parameter 1 (beta pdfs)
    # shape parameter 1 lower bound (beta pdfs)
    # shape parameter 1 upper bound (beta pdfs)
    # number of levels or samples
    # total number of evaluations
    # pdftype
    # algorithm type

    $v{mean}   = $mean;
    $v{min}    = $minval;
    $v{max}    = $maxval;
    $v{sigma}  = $sigma;
    $v{minsig} = $minsigma;
    $v{maxsig} = $maxsigma;

    $pconfig{$param} = \%v;

    if ( $verbose ) {
      if ( length($param) > $plength ) {
        $plength = length($param);
      }
    }
  }
}

undef $fh; # done parsing database file

$n_vars = @var_names;

# check if file looks good

die "$FAILED: mismatch in number of cases (should be $n_evals), stopped"
    unless ( keys %cases == $n_evals );

foreach my $case ( keys %cases ) {
  die "$FAILED: mismatch in random variable values for realization $case, stopped"
      unless ( @{$cases{$case}} == $n_vars );
}

if ( $verbose && $plength > 0 ) {
  print " o Random variables:\n";
  foreach ( sort keys %pconfig ) {
    my $par = $pconfig{$_};
    printf "   %-${plength}s -> (%g,%g)\n", $_, $par->{mean}, $par->{sigma};
  }
}

mkdir 'cases';
chdir 'cases' or die "$FAILED cannot 'chdir cases', stopped";

# set up CFD database
my $db;
$db = cfdDB() if @cfd_cases;

# create a directory for each realization, i.e. each case to run

print " o Creating directories for each realization\n" if $verbose;

foreach my $case ( sort keys %cases ) {
  my $dirname = 'case.'.$case;
  mkdir $dirname;
  chdir $dirname or die "$FAILED cannot chdir $dirname, stopped";
  $fh = FileHandle->new( "$ofile", '>' );
  die "$FAILED Cannot write $ofile, stopped" unless ( defined $fh );
  print $fh "# QUEST Random Variables\n";
  print $fh "# Created by $PROGRAM_NAME on ";
  print $fh scalar localtime;
  print $fh "\n";
  print $fh "# Random_Variable_Name Random_Variable_Value\n";
  foreach my $par ( 0 ... $#{$cases{$case}}) {
    print $fh $var_names[$par]," ",$cases{$case}[$par],"\n";
  }
  undef $fh;

  if ($db) { # link correct CFD run
    my $cfdir = 'case.00.' . sprintf("%05d",$case2cfd{$case});
    symlink "../../../C612A-P1NC-uq/cases/$db/$cfdir", "CFD";
    if (lstat "CFD" and not stat "CFD") {
      print "ERROR: problems linking to CFD database\n";
      exit 1;
    }
  }

  chdir File::Spec->updir();
}

chdir File::Spec->updir(); # out of cases

exit 0; # all done

sub flowParam {
  my $param = shift;
  my $rv    = 0;

  if ( $param =~ /Mach/ || $param =~ /alpha/ || $param =~ /beta/ ||
       $param =~ /gamma/ ||
       $param =~ /InletPressRatioBC/ ||
       $param =~ /InletVelocityBC/ ||
       $param =~ /RoL/ || $param =~ /Phi/ ) {
    $rv = 1;
  }

  return $rv;
}


# NOTE: Don't create this?
sub cfdDB {
  # name CFD database
  my $db = 'CFD';
  foreach ( sort keys %pconfig ) {
    my $par = $pconfig{$_};
    if ( flowParam($_) ) {
      $db .= "_$_"."_m". $par->{mean} . "_s" . $par->{sigma};
    }
  }

  print " o Creating CFD database $db\n" if $verbose;

  # create CFD database directories (nothing but the manifest in it)
  mkdir $db;
  chdir $db or die "$FAILED cannot chdir $db, stopped";

  # record manifest of runs
  my $fh = FileHandle->new( 'MANIFEST.txt', '>' );
  die "$FAILED Cannot write MANIFEST.txt, stopped" unless ( defined $fh );
  print $fh "# Directory-Parameter Pairs\n";
  print $fh "# Created by $PROGRAM_NAME on ";
  print $fh scalar localtime;
  print $fh "\n";

  my $i=0;
  foreach (@cfd_cases) {
    print $fh 'case.00.' . sprintf("%05d",$i), " $_\n";
    $i++;
  }
  undef $fh;

  chdir File::Spec->updir(); # return from CFD dir

  return $db;
}
