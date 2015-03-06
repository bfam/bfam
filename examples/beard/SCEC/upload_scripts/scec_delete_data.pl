#!/usr/bin/perl
# This script modified from:
#    http://www.willmaster.com/library/manage-forms/using_perl_to_submit_a_form.php

# This script deletes the all the data in the given scec test problem


use strict;
use warnings;

my $SCEC_CVWS = "http://scecdata.usc.edu/cvws/cgi-bin/cvws.cgi";

my $NARGS = $#ARGV+1;
if ($NARGS != 3 && $NARGS != 4)
{
  print "Usage:   scec_delete_data.pl USER_NAME PASSWORD PROBLEM_NUMBER\n";
  print "Usage:   scec_delete_data.pl USER_NAME PASSWORD USER_NUMBER PROBLEM_NUMBER\n";
  print "Example: scec_delete_data.pl user pass 31\n";
  print "Example: scec_delete_data.pl user pass 31 2\n";
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


print "Confirm the following delete information (y/n):\n";
print "   user name: $UNAME\n";
print "   user num:  $URNM\n";
print "   problem:   $TPV\n";

chomp(my $input = <STDIN>);
if(lc(substr($input,0,1)) ne 'y')
{
  print "cancelling with input: $input\n";
  exit 0;
}

my %LIST_DATA = (
  "G1090$URNM" => "Select",
  "u"          => $UNAME,
  "p"          => $PSSWD,
  "m"          => $TPV,
  "urv"        => $URNM,
  "o"          => "1005"
);

# It's a good habit to always use the strict module.
use strict;

# Modules with routines for making the browser.
use LWP::UserAgent;
use HTTP::Request::Common;

# Create the browser that will post the information.
my $Browser = new LWP::UserAgent;

# Post the information to the CGI program.
my $Page = $Browser->request(POST $SCEC_CVWS,\%LIST_DATA);

if ($Page->is_success)
{

  # Check the info provided
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

  # Delete the contour plot
  my %DELETE_CMD = (
    "G1036"=>"Delete",
    "u"=>$UNAME,
    "p"=>$PSSWD,
    "m"=>$TPV,
    "s"=>"cplot",
    "urv"=>$URNM,
    "o"=>"1005",
  );
  print "DELETE Contour Plot\n";
  my $PageDelete = $Browser->request(POST $SCEC_CVWS,\%DELETE_CMD);

  # Delete all the station data
  foreach (split(/\n/,$Page->content))
  {
    $_ =~ /"G1027([^"]+)"/ || next;
    my ($station) = ($1);
    print "Deleting $station....";

    my %DELETE_CMD = (
      "G1025"=>"Delete",
      "u"=>$UNAME,
      "p"=>$PSSWD,
      "m"=>$TPV,
      "s"=>$station,
      "urv"=>$URNM,
      "o"=>"1005",
    );
    my $PageDelete = $Browser->request(POST $SCEC_CVWS,\%DELETE_CMD);

    if(index($PageDelete->content,"File Deleted") != -1)
    {
      print "sucessful!\n";
    }
    elsif(index($PageDelete->content,"Data File Not Found") != -1)
    {
      print "fail (Data File Not Found)\n";
    }
    else
    {
      print "fail (unknown reason)\n";
    }
  }
}
else { print $Page->message; }

# end of script

