#!/usr/bin/env perl
#
# based on aregier's scripts/annotate_original.pl 
# from: ssh://git/srv/git/gatk-workflow.git

use strict;
use warnings;
use Data::Dumper;

my $original = $ARGV[0];
my $cadd = $ARGV[1];

open(CADD, "<$cadd");

my $scores;
my $raw_scores;
my $header1;
my $header2;
while (my $line = <CADD>) {
    chomp $line;
    if ($line =~ "^##INFO=<ID=CADD,") {
        $header1 = $line;
    }
    elsif ($line =~ "^##INFO=<ID=CADD_RAW,") {
        $header2 = $line;
    }
    elsif ($line =~ "^#") {
    }
    else {
        my @fields = split("\t", $line);
        my @infos = split(";", $fields[7]);
        my $ref = $fields[3];
        my $alt = $fields[4];
        my $chr;
        my $start;
        my $cadd;
        my $cadd_raw;
        my %info_dict = map{split("=", $_)} grep {$_ =~ "="} @infos;
        $cadd = $info_dict{"CADD"};
        $cadd_raw = $info_dict{"CADD_RAW"};
        $chr = $info_dict{"OriginalContig"};
        $start = $info_dict{"OriginalStart"};
        die unless defined $cadd and defined $cadd_raw and defined $chr and defined $start;
        $scores->{$chr}->{$start}->{$ref}->{$alt} = $cadd;
        $raw_scores->{$chr}->{$start}->{$ref}->{$alt} = $cadd_raw;
    }
}
close(CADD);

open(VCF, "<$original");

while (my $line = <VCF>) {
    if ($line =~ "^##") {
        print $line;
    }
    elsif ($line =~ "^#CHROM") {
        print "$header1\n$header2\n$line";
    }
    else {
        chomp $line;
        my @fields = split("\t", $line);
        my $cadd;
        my $cadd_raw;
        if (defined $scores->{$fields[0]}->{$fields[1]}->{$fields[3]}->{$fields[4]}) {
            $cadd = "CADD=".$scores->{$fields[0]}->{$fields[1]}->{$fields[3]}->{$fields[4]};
            $cadd_raw = "CADD_RAW=".$raw_scores->{$fields[0]}->{$fields[1]}->{$fields[3]}->{$fields[4]};
        }
        elsif (defined $scores->{$fields[0]}->{$fields[1]}->{rev_comp($fields[3])}->{rev_comp($fields[4])}) {
            $cadd = "CADD=".$scores->{$fields[0]}->{$fields[1]}->{rev_comp($fields[3])}->{rev_comp($fields[4])};
            $cadd_raw = "CADD_RAW=".$raw_scores->{$fields[0]}->{$fields[1]}->{rev_comp($fields[3])}->{rev_comp($fields[4])};
        }
        else {
            $cadd = "CADD=.";
            $cadd_raw = "CADD_RAW=.";
        }
        $fields[7] = join(";", $fields[7], $cadd, $cadd_raw);
        print join("\t", @fields)."\n";
    }
}

sub rev_comp {
    my $allele = shift;
    $allele =~ tr/atcgATCG/tagcTAGC/;
    reverse($allele);
}
