use warnings;
use strict;


my %count_overall_kmer;
my $region_size_setting = 100;

my $kmer_size = 6;

my $bedfi = $ARGV[0];
my $species = $ARGV[1];

my $outfi = $bedfi.".".$kmer_size."mer.txt";
open(OUT,">$outfi");


&count_kmers($bedfi,"clip");

for my $kmer (keys %count_overall_kmer) {
    next if ($kmer eq "all");
    print OUT "$kmer\t".$count_overall_kmer{$kmer}{"clip"}."\t".sprintf("%.5f",100*$count_overall_kmer{$kmer}{"clip"}/$count_overall_kmer{"all"}{"clip"})."\n";
}

close(OUT);

sub count_kmers {
    
    my $bedfi = shift;
    my $label = shift;

    my %data_by_chr;
    open(F,$bedfi);
    while (<F>) {
        chomp($_);
	
	my @tmp = split(/\t/,$_);
	my $chr = $tmp[0];
	my $pos = $tmp[1]."-".$tmp[2];
	my $str = $tmp[5];
	my $region = $chr.":".$pos.":".$str;
#        my ($chr,$start,$stop,$count) = split(/\t/,$_);
        push @{$data_by_chr{$chr}},$region;
    }
    close(F);

    for my $chr (keys %data_by_chr) {
        print STDERR "doing $chr\n";

        my $chr_fi = "/storage/vannostrand/genomes/".$species."/chromosomes/".$chr.".fa";
        

	if (-e $chr_fi) {
	    open(C,$chr_fi) || die "no $chr_fi\n";
	} else {
	    my $chr_fi = "/storage/vannostrand/genomes/".$species."/chromosomes/".$chr.".fa.gz";
	    if (-e $chr_fi) {
		open(C,"gunzip -c $chr_fi |") || die "no $chr_fi\n";
	    } else {
		print STDERR "dying - couldn't find $chr_fi for chromosome $chr\n";
	    }
	}

        my $chr_seq;
        while (<C>) {
            chomp($_);
            next if ($_ =~ /^\>/);
            $chr_seq .= uc($_);
        }
        close(C);

        for my $bed_region (@{$data_by_chr{$chr}}) {
            my ($bedchr,$bedpos,$bedstrand) = split(/\:/,$bed_region);
            my ($bedstart,$bedstop) = split(/\-/,$bedpos);
	    my $region;                                                                                                                 

#	    print "chr $chr bedstart $bedstart bedstop $bedstop\n";

	    if ($bedstrand eq "+") {   
		$region = substr($chr_seq,$bedstart,$bedstop-$bedstart);
	    } elsif ($bedstrand eq "-") {                                                                                               
		$region = &rc(substr($chr_seq,$bedstart,$bedstop-$bedstart));
	    } else {
		print STDERR "strand error $bedstrand\n";
	    }

	    unless ($region) {
		print STDERR "error - no sequence $chr $chr_fi $bed_region\n";
	    }
	    
	    for (my $i=0;$i<length($region)-$kmer_size+1;$i++) {
		my $kmer = substr($region,$i,$kmer_size);
		next if ($kmer =~ /N/);

		$count_overall_kmer{$kmer}{$label}++;
		$count_overall_kmer{"all"}{$label}++;
	    }
	}
    }
}




sub rc {
    my $seq = shift;
    my $rc = reverse($seq);
    $rc =~ tr/ACGT/TGCA/;
    return($rc);
}


