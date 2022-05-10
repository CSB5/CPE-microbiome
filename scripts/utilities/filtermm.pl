#!/usr/bin/perl
use strict;
use warnings;


my ($nb_mismatch) = @ARGV;

while (<STDIN>) {

    chomp $_;

   if ($_ =~ /^\@/){
        print "$_\n";
    }else{

        my @words = split /\s+/, $_;
        my $mdz = $words[12];
        my @segmt = split /[0-9]+/, $mdz;

#        print STDERR "$mdz $#segmt\n";
        if($#segmt<=$nb_mismatch){
            print "$_\n";
        }
    }

}
