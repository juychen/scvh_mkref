#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long qw(GetOptions);

#########################################################################################
BEGIN { unshift @INC, "/code/utilities"; }

use FAlite;
#use DataBrowser;

my $samtools = "samtools";
my $gffread = "/code/utilities/gffread/gffread";
my $minimap2 = "/code/utilities/minimap2/minimap2";
my $dir = "/";
my $human_fa = "/code/utilities/hg38.fa";
my $human_gtf = "/code/utilities/hg38.unique_gene_names.gtf";
my $source = "viruSITE.NCBIprokaryotes";
my $viruSITE = "/code/viruSITE_human_host.txt";
my $prokaryotes = "/code/prokaryotes.csv";

#######################################################################################
my $output_prefix = "human_host_viruses_microbes";
my $min_length_exon = 50;
my $threads = 20;
my $seq_divergence = 20;
my $output_dir = "$dir/human_host_viruses_reference_set";
mkdir $output_dir unless (-d $output_dir);
my $output_dir_fa = "$output_dir/fa";
mkdir $output_dir_fa unless (-d $output_dir_fa);
my $output_dir_gtf = "$output_dir/gtf";
mkdir $output_dir_gtf unless (-d $output_dir_gtf);

GetOptions(
	'o|output_dir=s' => \$output_dir,
	't|threads=i' => \$threads,
	'd|database=s' => \$source,
	'human_fa=s' => \$human_fa,
	'human_gtf=s' => \$human_gtf,
	'viruSITE=s' => \$viruSITE,
	'prokaryotes=s' => \$prokaryotes,
	'output_prefix=s' => \$output_prefix,
	'min_length_exon=i' => \$min_length_exon,
	'seq_divergence=i' => \$min_length_exon,
) or die_usage();

sub die_usage {
die "
Usage: scvh_map_reads.pl [Options] <vmh_genome_dir> <R2> <R1> or <vmh_genome_dir> <.bam file>

Options:                                                                                                                                Defaults
-o/--output-dir	<string>   	the output directory                                                                                          [<$output_dir>]
-t/--threads <int>         	number of threads to run alignment with                                                                       [<$threads>]
-d/--database <string>     	select virus or virus and prokaryotes database, can be 'viruSITE' or 'viruSITE.NCBIprokaryotes'               [<$source>]
--output_prefix <string>  	Prefix of the output file, can be 'human_host_viruses' or 'human_host_viruses_microbes'                       [<$output_prefix>]
--human_fa <string>  		Input path of the human fa file                                                          					  [<$human_fa>]
--human_gtf <string>  		Input path of the human gtf file                                                          					  [<$human_gtf>]
--viruSITE <string>  		Input path of the viruSITE accesion file                                                          	  		  [<$viruSITE>]
--prokaryotes <string>  	Input path of the prokaryotes acccesion file                                                          		  [<$prokaryotes>]
--min_length_exon <int>  	minimap2 param:  See --min_length_exon in minimap2 manual                              			  			  [<$min_length_exon>]
--seq_divergence <int>  	minimap2 param:  See --seq_divergence in minimap2 manual                              						  [<$seq_divergence>]
";
}

download($viruSITE);
download($prokaryotes);

my $set = "with_hg38";

merge($set);

remove_amb_viral_exon_annot_using_minimap2($set);

#die;


sub remove_amb_viral_exon_annot_using_minimap2 {
	my $set = shift;
	#extract viral exon seqs
	my $viral_fa = "$output_prefix.$source.fa";
	my %fas;
	open IN, "<$output_dir/$viral_fa" or die "can't open $output_dir/$viral_fa\n";
	my $fasta = new FAlite(\*IN);

	while(my $entry = $fasta->nextEntry) {
		my $header = $entry->def;
		my ($ref) = $header =~ /^>(\S+)/;
	
		my $seq = $entry->seq;
		$fas{$ref} = $seq;
	}
	close IN;
	
	my $viral_gtf = "$output_prefix.$source.gtf";
	open my $IN, "<$output_dir/$viral_gtf" or die "can't open $output_dir/$viral_gtf\n";
	
	my $feature = "exon";
	
	my %remove_due_to_exceed_ref_coord;
	
	my $output_file = $output_dir . "/$output_prefix.$source.$feature.fa";
	open my $OUT, ">$output_file" or die "can't open $output_file\n";
	while (<$IN>) {
		my @line = split("\t", $_);
		#next unless ($line[2] eq $feature);
		if ($line[2] ne $feature) {
			
			my ($ref,$start,$end) = ($line[0],$line[3],$line[4]);
			if ( ($start > length($fas{$ref})) ||  ($end > length($fas{$ref})) ) {
				print "$start-$end in $ref not defined since $ref is " . length($fas{$ref}) . " bp long\n";
				$remove_due_to_exceed_ref_coord{"$ref:$start-$end"}++;
			}
			next;
		}
	
		my ($ref,$start,$end) = ($line[0],$line[3],$line[4]);
		my $length = $end - $start + 1;
		my ($gene) = $_ =~ /gene_id\s+\"(.*?)\"/;
	
		my $name = "$ref:$start-$end" . "_" . $gene;
		
		if ($length < $min_length_exon) {
			print "$name with length $length is < $min_length_exon, will not be evaluated by minimap2\n";
			next;
		}
		
		if ( ($start > length($fas{$ref})) ||  ($end > length($fas{$ref})) ) {
			print "$gene with $start-$end in $ref not defined since $ref is " . length($fas{$ref}) . " bp long\n";
			$remove_due_to_exceed_ref_coord{"$ref:$start-$end"}++;
			next;
		}
		
		my $sequence = substr($fas{$ref}, $start-1, $length);
		print $OUT ">$name\n$sequence\n";
	
	}
	close $IN;
	close $OUT;
	
	my $minimap2_output_sam = $output_dir . "/hg38.$source.$feature.sam";
	
	my $preset;
	if ($seq_divergence == 5) {
		$preset = "asm5";
	} elsif ($seq_divergence == 20) {
		$preset = "asm20";
	}
	
	################
	#map to $human_fa

	#map to $output_dir_set
	my $output_dir_set = "$output_dir/$set";
	mkdir $output_dir_set unless (-d $output_dir_set);

	my $minimap2_output_bam = $output_dir . "/hg38.$source.$feature.bam";
	my $minimap2_output_sorted_bam = $output_dir . "/hg38.$source.$feature.sorted.bam";
	#######################

	system("$minimap2 -ax $preset -t $threads $human_fa $output_file > $minimap2_output_sam");
	system("$samtools view -Sb -F 4 $minimap2_output_sam > $minimap2_output_bam");
	system("$samtools sort -o $minimap2_output_sorted_bam $minimap2_output_bam");
	system("$samtools index $minimap2_output_sorted_bam");
	
	my %to_remove;
	open $IN, "$samtools view $minimap2_output_sorted_bam |" or die "can't open $minimap2_output_sorted_bam\n";
	while (<$IN>) {
		chomp;
		my @line = split("\t", $_);
		my $exon_name = $line[0];
		my ($acc, $start, $end, $gene_id) = $exon_name =~ /(.*?):(\d+)-(\d+)_(.*)/;
		$to_remove{"acc"} = $acc;
		$to_remove{"start"} = $start;
		$to_remove{"end"} = $end;
		$to_remove{"gene_id"} = $gene_id;
	}
	close $IN;
	
	my $human_viral_gtf = "$output_dir_set/$output_prefix.$source.$set.gtf";
	$output_file = "$output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.temp";
	open $OUT, ">$output_file" or die "can't open $output_file\n";
	open $IN, "<$human_viral_gtf" or die "can't open $human_viral_gtf\n";
	while (<$IN>) {
		chomp;
		my @line = split("\t", $_);
		my ($ref,$start,$end) = ($line[0],$line[3],$line[4]);
		my $locus = "$ref:$start-$end";
		my ($gene) = $_ =~ /gene_id\s+\"(.*?)\"/;
		my ($transcript) = $_ =~ /transcript_id\s+\"(.*?)\"/;
		if ( ($ref eq $to_remove{"acc"}) && ($start == $to_remove{"start"}) && ($end == $to_remove{"end"}) && ($gene eq $to_remove{"gene_id"}) ) {
			print "Removed $_ for mapping to host" . "\n";
			next;
		} elsif (exists $remove_due_to_exceed_ref_coord{$locus}) {
			print "Removed $_ for start or end coords exceeding the reference\n";
			next;
		} elsif ( ($gene =~ /\s/) || ($transcript =~ /\s/) ) { 
			print "Removed $_ for gene or transcript id containing whitespace\n";
		} else {
			print $OUT $_ . "\n";
		}
	}
	close $IN;
	close $OUT;
	
	my $chr_prefix = "chr";
	
	
	my %transcript_id_to_gene_id;
	my %transcript_id_to_gene_name;
	
	my %transcript_id_feature_type;

	my %transcript_id_to_strand;
	my %transcript_id_to_coords;
	my %transcript_id_to_ref;
	my %transcript_id_to_source;
	#add "transcript" feature back to GTF for transcript_ids without the transcript feature (removed due to coords exceeding)
	open $IN, "<$output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.temp" or die "can't open $output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.temp\n";
	while (<$IN>) {
		my @line = split("\t", $_);
		my ($ref, $source, $feature_type, $start, $end, $strand) = ($line[0], $line[1], $line[2], $line[3], $line[4], $line[6]);
		
		if ($ref =~ /^$chr_prefix/) {
			next;
		}
		
		my ($gene_name_name_portion) = $_ =~ /gene_name \"(.*?)\";/;
		my ($gene_id_name_portion) = $_ =~ /gene_id \"(.*?)\";/;
		my ($transcript_id_name_portion) = $_ =~ /transcript_id \"(.*?)\";/;

		if ($gene_name_name_portion) {
			$transcript_id_to_gene_name{$transcript_id_name_portion} = $gene_name_name_portion;
		}
		
		$transcript_id_to_gene_id{$transcript_id_name_portion} = $gene_id_name_portion;
		
		$transcript_id_feature_type{$transcript_id_name_portion}{$feature_type}++;
		
		$transcript_id_to_strand{$transcript_id_name_portion} = $strand;
		
		$transcript_id_to_coords{$transcript_id_name_portion}{$start}++;
		$transcript_id_to_coords{$transcript_id_name_portion}{$end}++;
		
		$transcript_id_to_ref{$transcript_id_name_portion} = $ref;
		
		$transcript_id_to_source{$transcript_id_name_portion} = $source;
	}
	close $IN;
	
	#browse(\%transcript_id_feature_type);
	#die;
	
	#my $appended_transcript_entries;
	
	open $OUT, ">$output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.appended_transcripts" or die "can't open $output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.appended_transcripts\n";
	
	for my $transcript_id (keys %transcript_id_feature_type) {
		if (exists $transcript_id_feature_type{$transcript_id}{"transcript"}) {
			next;
		}
		
		if (!exists $transcript_id_feature_type{$transcript_id}{"exon"}) {
			next;
		}
		
		#lowest start and highest end
		my @coords = sort { $a <=> $b } keys %{$transcript_id_to_coords{$transcript_id}};
		#print "min: $coords[0]\n";
		#print "max: $coords[-1]\n";
		
		my $transcript_gtf_line = $transcript_id_to_ref{$transcript_id} . "\t" . $transcript_id_to_source{$transcript_id} . "\t" . "transcript" . "\t" . $coords[0] . "\t" . $coords[-1] . "\t"
		. ".\t" . $transcript_id_to_strand{$transcript_id} . "\t.\t" . "transcript_id \"$transcript_id\"; gene_id \"" . $transcript_id_to_gene_id{$transcript_id} . "\";";
		
		if ($transcript_id_to_gene_name{$transcript_id}) {
			$transcript_gtf_line .= " gene_name \"" . $transcript_id_to_gene_name{$transcript_id} . "\";";
		}
		#$appended_transcript_entries .= "$transcript_gtf_line\n";
		print $OUT $transcript_gtf_line . "\n";
		
	}
	
	close $OUT;
	
	system("cat $output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.temp $output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf.appended_transcripts > $output_dir_set/$output_prefix.$source.$set.removed_amb_viral_exon.gtf");
	
}

sub merge {
	my $set = shift;
	my $output_dir_set = "$output_dir/$set";
	mkdir $output_dir_set unless (-d $output_dir_set);

	#copy 10x barcodes whitelist
	# my $whitelist = "/home/asdfken/data/Cellranger_barcodes_whitelist/737K-august-2016.txt";
	# system("cp $whitelist $output_dir_set");
	
	my $viral_gtf_temp = "$output_prefix.$source.gtf.temp";

	my @gtf_files = `ls $output_dir_gtf/*force_gene_id.gtf`;
	chomp(@gtf_files);

	
	my $merge_gtf_command = "cat " . join(" ", @gtf_files) . " > $output_dir/$viral_gtf_temp";
	system("$merge_gtf_command");
	
	#check gene_ids and gene_names are actually unique, and make them all unique if not already
	my $viral_gtf = check_viral_gtf_ids_names_unique($viral_gtf_temp);
	
	my $viral_fa = "$output_prefix.$source.fa";
	my @fa_files = `ls $output_dir_fa/*.fa`;
	chomp(@fa_files);
	my $merge_fa_command = "cat " . join(" ", @fa_files) . " > $output_dir/$viral_fa";
	system("$merge_fa_command");

	#merge human with virus ref
	my $human_viral_gtf = "$output_prefix.$source.$set.gtf";
	my $human_viral_fa = "$output_prefix.$source.$set.fa";
	system("cat $human_fa $output_dir/$viral_fa > $output_dir_set/$human_viral_fa");
	system("cat $human_gtf $output_dir/$viral_gtf > $output_dir_set/$human_viral_gtf");
}

sub check_viral_gtf_ids_names_unique {
	my $viral_gtf_temp = shift;
	
	my %uniqueness;
	
	open my $IN, "<$output_dir/$viral_gtf_temp" or die "can't open $output_dir/$viral_gtf_temp\n";
	while (<$IN>) {
		my @line = split("\t", $_);
		my ($acc, $feature) = ($line[0], $line[2]);
		
		next unless ($feature eq "transcript");
		
		my ($gene_name_name_portion) = $_ =~ /gene_name \"(.*?)\";/;
		my ($gene_id_name_portion) = $_ =~ /gene_id \"(.*?)\";/;
		
		if ($gene_name_name_portion) {
			$uniqueness{"gene_name"}{$gene_name_name_portion}{$acc}{$_}++;
		}
		
		if ($gene_id_name_portion) {
			$uniqueness{"gene_id"}{$gene_id_name_portion}{$acc}{$_}++;
		}
	}
	close $IN;
	
	my %need_changing;
	

	for my $name_type (keys %uniqueness) {
		print $name_type . ":\n";
		for my $name_portion (keys %{$uniqueness{$name_type}}) {
			my $times_appeared = scalar(keys  %{$uniqueness{$name_type}{$name_portion}});
			if ($times_appeared > 1) {
				#print join("", (keys  %{$uniqueness{$name_type}{$name_portion}})) . "\n";
				print "$name_portion appears in " . join(" ", (keys  %{$uniqueness{$name_type}{$name_portion}})) . "\n";
				
				$need_changing{$name_type}{$name_portion}++;
			}
		}
	}
	
	my $viral_gtf = "$output_prefix.$source.gtf";
	open my $OUT, ">$output_dir/$viral_gtf" or die "can't open $output_dir/$viral_gtf\n";
	open $IN, "<$output_dir/$viral_gtf_temp" or die "can't open $output_dir/$viral_gtf_temp\n";
	while (<$IN>) {
		chomp;
		my @line = split("\t", $_);
		my ($acc, $feature) = ($line[0], $line[2]);
		
		my $final_out_line .= $_;
		
		my ($transcript_id_name_portion) = $_ =~ /transcript_id \"(.*?)\";/;
		my ($gene_name_name_portion) = $_ =~ /gene_name \"(.*?)\";/;
		my ($gene_id_name_portion) = $_ =~ /gene_id \"(.*?)\";/;
		my ($transcript_id) = $_ =~ /transcript_id (\".*?\";)/;
		my ($gene_name) = $_ =~ /gene_name (\".*?\";)/;
		my ($gene_id) = $_ =~ /gene_id (\".*?\";)/;
		
		if (exists $need_changing{"gene_id"}{$transcript_id_name_portion}) {
			my $replacement = "\"$transcript_id_name_portion.$acc\";";
			print "in $acc: replaced transcript_id $transcript_id with $replacement\n";
			
			$final_out_line =~ s/transcript_id $transcript_id/transcript_id $replacement/;
		}
		
		if (exists $need_changing{"gene_id"}{$gene_id_name_portion}) {
			my $replacement = "\"$gene_id_name_portion.$acc\";";
			print "in $acc: replaced gene_id $gene_id with $replacement\n";
			
			$final_out_line =~ s/gene_id $gene_id/gene_id $replacement/;
		}
		
		if ($gene_name) {
			if (exists $need_changing{"gene_name"}{$gene_name_name_portion}) {
				my $replacement = "\"$gene_name_name_portion.$acc\";";
			
				print "in $acc: replaced gene_name $gene_name with $replacement\n";
			
				$final_out_line =~ s/gene_name $gene_name/gene_name $replacement/;
			}
		}
		
		print $OUT $final_out_line . "\n";
	}
	close $IN;
	close $OUT;
	
	system("rm $output_dir/$viral_gtf_temp");
	return $viral_gtf;
	
}

sub download {
	my $input_accs = shift;
	print "Open $input_accs\n";
	open my $IN, "<$input_accs" or die "can't open input $input_accs \n";
	
	my $type;
	if ($input_accs =~ /viruSITE/) {
		$type = "virus";
	} else {
		$type = "microbe";
	}
	
	#viruSITE_human_host.txt
	#"Virus name"; "Genome length"; "Genome type"; "Genome segments"; "Taxonomy ID"; "RefSeq ID"; "Updated"
	
	#prokaryotes.csv
	##Organism Name,Organism Groups,Strain,BioSample,BioProject,Assembly,Level,Size(Mb),GC%,Replicons,WGS,Scaffolds,CDS,Release Date,GenBank FTP,RefSeq FTP
	
	my $header = <$IN>;
	while (<$IN>) {
		chomp;
		my @line;
		my $acc;
		if ($type eq "virus") {
			@line = split(";", $_);
			$acc = $line[5];
		
			if ($acc =~ /[A-Z]+_\d+/) {

			} else {
				print "Skipped $acc\n";
				next;
			}
			
		} elsif ($type eq "microbe") {
			@line = split(",", $_);
			my $acc_info = $line[9];
			($acc) = $acc_info =~ /.*?:(.*?)\//;
			
			if ($acc =~ /[A-Z]+_.*/) {

			} else {
				print "Skipped $acc\n";
				next;
			}
		}

		system("curl \"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$acc\" --output $output_dir_gtf/$acc.gff");
			
		#die;
		
		#remove UTR
		my $output_gff = "$output_dir_gtf/$acc.noUTR.gff";
		open my $OUT_GFF, ">$output_gff" or die "can't open $output_gff\n";
		open my $IN_GFF, "<$output_dir_gtf/$acc.gff" or die "can't open $output_dir_gtf/$acc.gff";
		while (<$IN_GFF>) {
			if ($_ =~ /^#/) {
				print $OUT_GFF $_;
				next;
			}
			my @line = split("\t", $_);
			my $feature = $line[2];
			if ($feature) {
				if ($feature =~ /UTR/) {
					print "Skipped UTR $_";
				} else {
					print $OUT_GFF $_;
				}
			} else {
				print $OUT_GFF $_;
			}
		}
		close $IN_GFF;
		close $OUT_GFF;
		
		system("$gffread $output_gff -T --force-exons --keep-genes -o $output_dir_gtf/$acc.gtf");
	
		open my $OUT_GTF, ">$output_dir_gtf/$acc.force_gene_id.gtf";
		open my $IN_GTF, "<$output_dir_gtf/$acc.gtf";
		while (<$IN_GTF>) {
			chomp;
			my ($transcript_id) = $_ =~ /transcript_id (\".*?\";)/;
			if ($transcript_id) {
				
				#print "$transcript_id\n";
				#die;
				
				my $final_out_line;
				
				my ($gene_name_name_portion) = $_ =~ /gene_name \"(.*?)\";/;
				
				my ($gene_id_name_portion) = $_ =~ /gene_id \"(.*?)\";/;
				
				my ($gene_name) = $_ =~ /gene_name (\".*?\";)/;
				
				my ($gene_id) = $_ =~ /gene_id (\".*?\";)/;
				
				if($gene_id){
					$final_out_line .= $_;
				} else {
					print "appended $transcript_id to $_ in $acc\n";
					$final_out_line .= $_ . " gene_id " . $transcript_id;
					$gene_id = $transcript_id;
				}
				
				if ($gene_name) {
					if ($gene_name eq $gene_id) {
						
					} else {
						my $final_gene_name = "\"$gene_name_name_portion.$gene_id_name_portion\";";
						
						print "replaced gene_name $gene_name with $final_gene_name\n";
						
						$final_out_line =~ s/gene_name $gene_name/gene_name $final_gene_name/;
					}
				}
				
				print $OUT_GTF $final_out_line . "\n";
				
			} else {
				print "No transcript ID for $_ in $acc, skipped\n";
			}
		}

		close $OUT_GTF;
		close $IN_GTF;
		
		system("curl \"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=$acc\" --output $output_dir_fa/$acc.fa");

	}
	close $IN;
	
}
