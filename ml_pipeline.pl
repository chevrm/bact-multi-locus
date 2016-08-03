#!/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use Bio::SeqIO;

## usage:
##  perl ml_pipeline.pl /PATH/TO/GENOMES/*.fna

## Setup paths
my $script_dir = abs_path($0);
my @p = split(/\//, $script_dir);
$script_dir =~ s/$p[-1]//;
my $hmmdb = $script_dir.'genprop0799/genprop0799.hmmdb';
my $cutfi = $script_dir.'genprop0799/genprop0799.cutoffs.tsv';

## Read in pHMM cutoffs
my $trust = 1; ## Set to 0 to use noise cutoff, 1 to use trusted cutoff (default=1)
my %cutoff = ();
open my $cfh, '<', $cutfi or die $!;
while(<$cfh>){
    next if($_ =~ m/^#/);
    chomp;
    my ($fam, $t, $n) = split(/\t/, $_);
    if($trust == 1){
	$cutoff{$fam} = $t;
    }else{
	$cutoff{$fam} = $n;
    }
}
close $cfh;
my @famorder = sort keys %cutoff;

## Read in pHMM lengths
my %lenof = ();
my $lastseen = undef;
open my $dfh, '<', $hmmdb or die $!;
while(<$dfh>){
    chomp;
    if($_ =~ m/^NAME\s+(\S+)/){
	$lastseen = $1;
    }elsif($_ =~ m/^LENG\s+(\d+)/){
	$lenof{$lastseen} = $1;
    }
}
close $dfh;

## Process each genome
system("rm cat.faa") if(-e 'cat.faa');
my %nuc = ();
foreach my $genome (@ARGV){
    my @a = split(/\//, $genome);
    my $pref = $a[-1];
    $pref =~ s/\.fna//;
    print "Processing $pref\n";
    ## Predict ORFs
    print "\tPredicting ORFs...";
    system("prodigal -t $pref.ptrain -c -i $genome > /dev/null 2>&1");
    system("prodigal -c -i $genome -a $pref.faa -d $pref.orf -t $pref.ptrain > /dev/null 2>&1");
    print "DONE!\n";
    ## Run HMM search
    print "\tAssigning gene families...";
    system("hmmscan --tblout tmp.hmmtbl.out --noali $hmmdb $pref.faa > /dev/null");
    ## Parse hmm tbl
    open my $hfh, '<', 'tmp.hmmtbl.out' or die $!;
    my %famhit = ();
    while(<$hfh>){
	next if($_ =~ m/^#/);
	chomp;
	my ($target, $tacc, $query, $qacc, $fullevalue, $fullscore, @rest) = split(/\s+/, $_);
	if($cutoff{$target} <= $fullscore){
	    if(exists $famhit{$target}){
		if($famhit{$target}{'score'} < $fullscore){
		    $famhit{$target}{'score'} = $fullscore;
		    $famhit{$target}{'id'} = $query;
		}
	    }else{
		$famhit{$target} = {
		    'score' => $fullscore,
		    'id' => $query
		};
	    }
	}
    }
    close $hfh;
    print "DONE!\n";
    ## Pull sequences and make prealignment string
    print "\tGenerating concatenated sequence...";
    my %porf = ();
    my $pfa = new Bio::SeqIO(-file=>"$pref.faa", -format=>'fasta');
    while(my $seq = $pfa->next_seq){
	my $s = $seq->seq;
	$s =~ s/\*//g;
	$porf{$seq->id} = $s;
    }
    open my $cfh, '>>', 'cat.faa' or die $!;
    print $cfh '>'."$pref\n";
    my %norf = ();
    my $nfa = new Bio::SeqIO(-file=>"$pref.orf", -format=>'fasta');
    while(my $seq = $nfa->next_seq){
	my $s = $seq->seq;
	$s =~ s/\w{3}$//;
	$norf{$seq->id} = $s;
    }
    my $nseq = '';
    foreach my $fam (@famorder){
	if(exists $famhit{$fam}){
	    print $cfh $porf{$famhit{$fam}{'id'}};
	    $nseq .= $norf{$famhit{$fam}{'id'}};
	}else{
	    foreach(1..$lenof{$fam}){
		print $cfh 'X';
		$nseq .= 'NNN';
	    }
	}
    }
    print $cfh "\n";
    close $cfh;
    $nuc{$pref} = $nseq;
    print "DONE!\n"
}

## Make alignment
print "\nGenerating alignment...";
system("mafft --quiet --namelength 90 cat.faa > cat.afa");
## Make nucl alignment
my $afa = new Bio::SeqIO(-file=>'cat.afa', -format=>'fasta');
open my $aln, '>', 'cat.aln' or die $!;
while(my $seq = $afa->next_seq){
    my @pseq = split(//, $seq->seq);
    my @nseq = ( $nuc{$seq->id} =~ m/.{3}/g );
    my $nucaln = '';
    my $ni = 0;
    for (my $pi = 0;$pi<scalar(@pseq);$pi+=1){
	if($pseq[$pi] eq '-'){
	    $nucaln .= '---';
	}else{
	    $nucaln .= $nseq[$ni];
	    $ni += 1;
	}
    }
    print $aln '>'.$seq->id."\n$nucaln\n";
}
close $aln;
print "DONE!\n";

## Tree
system("rm RAxML*") if(-e 'RAxML_bestTree.cat.tree');
my $bootstraps = 1000;
print "Generating tree...";
system("raxmlHPC -s cat.aln -n cat.tree -x 3524627 -p 3524627 -f a -y -N $bootstraps -m GTRGAMMA");
print "DONE!\n";

## Cleanup
my @rm = (
    '*.orf', '*.faa', '*.ptrain',
    'tmp*'
    );
foraech(@rm){
    system("rm $_");
}
