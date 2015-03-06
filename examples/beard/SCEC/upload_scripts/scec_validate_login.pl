#!/usr/bin/perl
# This script modified from:
#    http://www.willmaster.com/library/manage-forms/using_perl_to_submit_a_form.php


# Validate the users scec login data

use strict;
use warnings;

my $SCEC_CVWS = "http://scecdata.usc.edu/cvws/cgi-bin/cvws.cgi";

my $NARGS = $#ARGV+1;
if ($NARGS != 3 && $NARGS != 4)
{
  print "Usage:   scec_validate_login.pl USER_NAME PASSWORD PROBLEM_NUMBER\n";
  print "Usage:   scec_validate_login.pl USER_NAME USER_NUMBER PASSWORD PROBLEM_NUMBER\n";
  print "Example: scec_validate_login.pl user pass 31\n";
  print "Example: scec_validate_login.pl user 2 pass 31\n";
  exit;
}

my $k = 0;
my $UNAME = $ARGV[$k];

my $URNM  = $UNAME;
if($NARGS == 4)
{
  $k = $k+1;
  $URNM  = "$UNAME.$ARGV[$k]";
}

$k = $k+1;
my $PSSWD = $ARGV[$k];

$k = $k+1;
my $TPV   = "tpv$ARGV[$k]";

my %LIST_DATA = (
  "G1090$URNM" => "Select",
  "u"          => $UNAME,
  "p"          => $PSSWD,
  "m"          => $TPV,
  "urv"        => $URNM,
  "o"          => "1005"
);

# Modules with routines for making the browser.
use LWP::UserAgent;
use HTTP::Request::Common;

# Create the browser that will post the information.
my $Browser = new LWP::UserAgent;

# Post the information to the CGI program.
my $Page = $Browser->request(POST $SCEC_CVWS,\%LIST_DATA);

if ($Page->is_success)
{

  if(index($Page->content,"Forbidden Operation") != -1)
  {
    print "Username or password is likely wrong\n";
    exit 1;
  }
  if(index($Page->content,"Benchmark Not Found") != -1)
  {
    print "Benchmark $TPV not found\n";
    exit 2;
  }
}
else
{
  print $Page->message;
  exit 3;
}

exit 0;

