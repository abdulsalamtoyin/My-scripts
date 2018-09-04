#!/usr/local/bin/perl -w

###################################################################################################################
# 	
#	CuffDiffFiltering.pl
#
#	Reads CuffDuff output *_exp.diff (gene_exp.diff) and organizes the data by genes (row)and comparison (column)
#	
# 	Andreas Gisel September 2015
#
###################################################################################################################

use strict;

my $file_in = shift;	# *_exp.diff file
my $file_out = shift;	# output file name
my $val_order = shift;	# order of the sample value in the output (B+D+E+A+C)
my $comp_order = shift;	# order of the expression comparison (foldchange) in the output (B-A+B-C+.....)
my $file_annot = shift; # annotation file

unless(open(IN, $file_in))
{
	print "cannot open $file_in!!\n";
	exit;
}

my %gene_locus;
my %gene_fc;
my %gene_val;

my @val_order = split/\+/, $val_order;
my @comp_order = split/\+/, $comp_order;

while(<IN>)
{
	next if(/text_id/);
	chomp;
	
	#test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
	#XLOC_000001	XLOC_000001	Manes.01G000100	Chromosome01:3-803	A	B	OK	1.84212	9.08311	2.30182	2.42948	0.0016	0.00857572	yes

	
	my @data = split/\t/;
	$gene_locus{$data[2]} = $data[3];
	# check for 'inf' foldchange
	if($data[9] eq "-inf")
	{
		$data[9] = -100;
	}
	elsif($data[9] eq "inf")
	{
		$data[9] = -100;
	}
	
	$gene_fc{$data[2]}{"$data[4]-$data[5]"} = "$data[6]\t$data[9]\t$data[11]\t$data[13]";
	$gene_val{$data[2]}{$data[4]} = $data[7];
	$gene_val{$data[2]}{$data[5]} = $data[8];
	
}

unless ( open(OUT, ">$file_out"))
{
	print "Cannot open file \"$file_out\"\n\n";
	exit;
}

# write header of the table
print OUT "Gene\tGene Locus\tChr\tStart\tEnd";

for my $gene ( keys %gene_val )
{
	foreach my $sample (@val_order)
#	for my $sample (sort keys %{ $gene_val{$gene}})
	{
		print OUT "\t$sample Expression";
	}
	last;
}
	
for my $gene (keys %gene_fc)
{
	foreach my $comp (@comp_order)
#	for my $comp ( sort keys %{ $gene_fc{$gene} } )
	{
		print OUT "\t$comp Check\t$comp foldchange\t$comp p-value\t$comp significant";
	}
	last;
}
print OUT "\n";

my %annot;
	
if($file_annot)
{
	unless(open(IN, $file_annot))
	{
		print "cannot open $file_annot!!\n";
		exit;
	}
	
	while (<IN>)
	{
		my $line = $_;
		chomp $line;
		my @data = split;
		# 32359217	Manes.01G000300	Manes.01G000300.1	Manes.01G000300.1.p	PF07851	PTHR21433	KOG4758			GO:0016021	AT4G10430.3		TMPIT-like protein
		$annot{$data[1]} = $line;
	}
}
	

while ( my ($gene,$locus) = each %gene_locus ) 
{
	$locus =~ /(Chromosome\d\d):(\d+)-(\d+)/;
	print OUT "$gene\t$locus\t$1\t$2\t$3";
	
	foreach my $sample (@val_order)
#	for my $sample ( sort keys %{ $gene_val{$gene} } )
	{
		print OUT "\t$gene_val{$gene}{$sample}";
	}
	
	foreach my $comp (@comp_order)
#	for my $comp ( sort keys %{ $gene_fc{$gene} } )
	{
		$comp =~ /(\w)-(\w)/;
		my $revcomp = "$2-$1";
		
		if($gene_fc{$gene}{$comp})
		{
			print OUT "\t$gene_fc{$gene}{$comp}";
		}
		elsif($gene_fc{$gene}{$revcomp}) # inverse order of comparison and value of FC
		{
			my @data = split/\t/,$gene_fc{$gene}{$revcomp};
			my $rev = $data[1]*-1;
			print OUT "\t$data[0]\t$rev\t$data[2]\t$data[3]";
		}
		else
		{
			print "Ops - wrong comparison!\n";
			print "$gene - $comp - $revcomp\n";
			#exit;
		}
	}
	
	if($file_annot)
	{
		if($annot{$gene})
		{
			print OUT "\t$annot{$gene}",
		}
	}
	
	print OUT "\n"
}


print "done";
