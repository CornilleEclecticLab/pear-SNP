use strict;
use warnings;
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use File::Path qw(make_path);  # Importing make_path function to create directories

# By Xilong CHEN

my ($te_gff);


GetOptions(
	"gff=s" => \$te_gff,
	);

# Create 'regions' directory if it doesn't exist
my $regions_dir = 'regions';
unless (-d $regions_dir) {
    make_path($regions_dir) or die "Failed to create directory $regions_dir: $!";
}

# First pass = calculating the coverage
my %te=();

my $gffin = Bio::Tools::GFF->new(-file => $te_gff, -gff_version => 3);
while( my $feature = $gffin->next_feature()) {
	next unless ($feature->source_tag =~ /REPET_TEs$/);
	if ($feature->primary_tag eq "match") {
		my ($name)=$feature->get_tag_values("Target");
		$name =~ s/\s+.*//;
		push @{$te{$name}}, $feature;
	}
}
$gffin->close();

foreach my $te (keys %te) {
	open BED, ">regions/$te.bed" or die "failed to open regions.$te.bed\n";
	foreach my $ft (@{$te{$te}}) {
		print BED join ("\t", $ft->seq_id, $ft->start(), $ft->end()), "\n";
	}
	close BED;
}
