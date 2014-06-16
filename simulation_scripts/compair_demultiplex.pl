#!/bin/perl

$reference=$ARGV[0];
$trimmed=$ARGV[1];
$readlength=$ARGV[2];
$barcode=$ARGV[3];

open REFERENCE, $reference or die $!;
open FASTQ, $trimmed or die $!;

my $count = 0;
my $correct = 0;
my $not_trimmed = 0;
my $not_needed = 0;
my $short_trimmed = 0;
my $questions = 0;
my $short_trim_needed = 0;
my $questions_trim_needed = 0;
my $wrong = 0;

my $reference_line = <REFERENCE>;
chomp $reference_line;
my $reference_sequence = <REFERENCE>;
chomp $reference_sequence;
my $fastq_line;
my $fastq_sequence;

while ($fastq_line = <FASTQ>){
chomp $fastq_line;
    $fastq_sequence = <FASTQ>;
chomp $fastq_sequence;
    my $rubish = <FASTQ>;
    $rubish = <FASTQ>;
    my $fastqloc = $fastq_line;
    $fastqloc =~ s/@//;
    $demBarcode = substr($fastqloc, rindex($fastqloc, "|") + 1);
    if (index($fastqloc, "_") > 0){
        $demBarcode = substr($demBarcode, 0, rindex($demBarcode, "_"));
    }
    $fastqloc = substr($fastqloc, 0, rindex($fastqloc, "-"));

    while ($fastqloc ne $reference_line){
        $reference_line = <REFERENCE>;
        chomp $reference_line;
        $reference_sequence = <REFERENCE>;
        chomp $reference_sequence;
    }

    $count = $count + 1;
    if ($barcode ne $demBarcode){
        $wrong = $wrong + 1;
    }

    if (length($reference_sequence) < length($fastq_sequence)){
        $not_trimmed = $not_trimmed + 1;

    }
    if (length($reference_sequence) > length($fastq_sequence) && length($fastq_sequence) < 81){
        $short_trimmed = $short_trimmed + 1;
if (length($reference_sequence)<$readlength){$short_trim_needed = $short_trim_needed + 1;}
    }
    if (length($reference_sequence) > length($fastq_sequence) && length($fastq_sequence) < $readlength && length($fastq_sequence) > 80){
        $questions = $questions + 1;
if(length($reference_sequence) < $readlength){$questions_trim_needed = $questions_trim_needed + 1;}
    }
    if (length($reference_sequence) == length($fastq_sequence) && length($fastq_sequence) < $readlength){
        $correct = $correct + 1;
    }
    if (length($reference_sequence) >= length($fastq_sequence) && length($fastq_sequence) >= $readlength){
        $not_needed = $not_needed + 1;
    }
}

#print "Total Count:\t" . $count . "\n";
#print "Not needed:\t" . $not_needed . "\n";
#print "Correct trim:\t" . $correct . "\n";
#print "Not trimmed:\t" . $not_trimmed . "\n";
#print "To Short:\t" . $short_trimmed . "\n";
#print "Questionable:\t" . $questions . "\n";

#print "\n\nShort trim needed:\t" . $short_trim_needed . "\n";
#print "questions trim needed:\t" . $questions_trim_needed . "\n";

print "" . $trimmed . "\t" . ($not_needed + $correct) . "\t" . $count . "\t" . $wrong . "\t";

close REFERENCE;
close FASTQ;
