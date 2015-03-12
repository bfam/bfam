#!/usr/bin/perl
# This script modified from:
#    http://www.willmaster.com/library/manage-forms/using_perl_to_submit_a_form.php

# This script deletes the all the data in the given scec test problem


use strict;
use warnings;

my $SCEC_CVWS = "http://scecdata.usc.edu/cvws/cgi-bin/cvws.cgi";

my $NARGS = $#ARGV+1;
if ($NARGS < 4)
{
  print "Usage:   scec_download_images.pl USER_NAME PASSWORD PROBLEM_NUMBER USERS\n";
  print "Example: scec_download_images.pl user pass 31\n";
  exit;
}

my $IMG_W = "1600";
my $IMG_H = "1600";

my $UNAME = $ARGV[0];

my $PSSWD = $ARGV[1];

my $TPV   = "tpv$ARGV[2]";

my $USERS = $ARGV[3];

for(my $k=4; $k < $NARGS; $k++)
{
  $USERS = "$USERS*$ARGV[$k]";
}

print "Confirm the following delete information (y/n):\n";
print "   user name: $UNAME\n";
print "   problem:   $TPV\n";
print "   users:     $USERS\n";

chomp(my $input = <STDIN>);
if(lc(substr($input,0,1)) ne 'y')
{
  print "cancelling with input: $input\n";
  exit 0;
}

my %LIST_DATA = (
  "G1090$UNAME" => "Select",
  "u"          => $UNAME,
  "p"          => $PSSWD,
  "m"          => $TPV,
  "urv"        => $UNAME,
  "o"          => "1005"
);

# Modules with routines for making the browser.
use LWP::UserAgent;
use LWP::Simple;
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

  foreach (split(/\n/,$Page->content))
  {
    $_ =~ /"G1027([^"]+)"/ || next;
    my ($station) = ($1);

    my @datatypes = ();
    if(index($station,"fault") != -1)
    {
      @datatypes = ( "h-shear-stress", "h-slip", "h-slip-rate", "n-stress", "v-shear-stress", "v-slip", "v-slip-rate");
    }
    elsif(index($station,"body") != -1)
    {
      @datatypes = ( "h-disp", "h-vel", "n-disp", "n-vel", "v-disp", "v-vel" );
    }

    foreach(@datatypes)
    {
      my $datatype = $_;
      my %DOWNLOAD_CMD = (
        "G1074".$datatype=>"Graph",
        "u"=>$UNAME,
        "p"=>$PSSWD,
        "m"=>$TPV,
        "o"=>"1005",
        "ppr"=>$IMG_H,
        "ppc"=>$IMG_W,
        "plw"=>"1",
        "psu"=>"1",
        "mus"=>$USERS,
        "mss"=>$station
      );

      my $filename = $TPV."_".$station."_".$datatype.".gif";
      print "Downloading station: $filename\n";
      my $PageDownload = $Browser->request(POST $SCEC_CVWS,\%DOWNLOAD_CMD);
      foreach (split(/\n/,$PageDownload->content))
      {
        $_ =~ /".*(http:\/\/scecdata.usc.edu\/cvws\/cgi-bin\/cvws.cgi\/cvws.gif[^"]*).*"/ || next;
        my ($url) = ($1);
        $url =~ s/amp;//g;
        getstore($url,$filename);
        $filename = "legend_".$filename;
      }
    }
  }

    {
      my %DOWNLOAD_CMD = (
        "G1080cplot"=>"Graph",
        "u"=>$UNAME,
        "p"=>$PSSWD,
        "m"=>$TPV,
        "o"=>"1005",
        "ppr"=>$IMG_H,
        "ppc"=>$IMG_W,
        "plw"=>"1",
        "psu"=>"1",
        "mus"=>$USERS,
      );

      my $filename = "trup.gif";
      print "Downloading trup: $filename\n";
      my $PageDownload = $Browser->request(POST $SCEC_CVWS,\%DOWNLOAD_CMD);
      foreach (split(/\n/,$PageDownload->content))
      {
        $_ =~ /".*(http:\/\/scecdata.usc.edu\/cvws\/cgi-bin\/cvws.cgi\/cvws.gif[^"]*).*"/ || next;
        my ($url) = ($1);
        $url =~ s/amp;//g;
        getstore($url,$filename);
        $filename = "legend_".$filename;
      }
    }
}
