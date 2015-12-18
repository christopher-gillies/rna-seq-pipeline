#!/usr/bin/perl
use strict;
use warnings;
print `mvn clean`;
print `mvn package -DskipTests=True`;

my $f =  `ls ./target/*.jar`;

$f =~ /\.\/target\/(rna-seq-pipeline-\d+\.\d+\.\d+)/;
chomp($f);
my $cmd = "mkdir -p ./release/;cp $f ./release/$1.jar\n";
print "$cmd";
print `$cmd`;