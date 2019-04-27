#!/usr/bin/perl
use Getopt::Std;
getopts("dhj:kn:p:tw");

# Print help information if requested
if ($opt_h) {
    print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
    print "\nOptions:\n";
    print "-d        (Don't duplicate frames that already exist\n";
    print "-h        (Print this information)\n";
    print "-j <n>    (Render every n frames)\n";
    print "-k        (Keep output directory)\n";
    print "-p <n>    (Render output in parallel using n procs)\n";
    print "-n <n>    (Render up to this number of files)\n";
    print "-t        (Add tracers)\n";
    print "-w        (Render PNGs only without making a movie)\n";
    exit 0;
}

die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

# Determine the number of processors
if(!defined $opt_p) {
    $uname=`uname`;
    if($uname=~/Linux/) {
        $nodes=`lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l`;
        chomp $nodes;
    } elsif($uname=~/Darwin/) {
        $nodes=`sysctl -n hw.physicalcpu_max`;
        chomp $nodes;
    } else {
        $nodes=4;
    }
} else {
    $nodes=$opt_p;
}

# Make output directory
$e=@ARGV[0];
$ebase=$e;
$ebase=~s/\.out$// or die "Filename should end in '.out'\n";
$odir=$ebase.".frames";
mkdir $odir unless -e $odir;

# Set color range
if($#ARGV==3) {
    $cb="[@ARGV[2]:@ARGV[3]]";
} else {
    $cb="[*:*]";
}

$P=1;$queue=1 if $opt_p==1;$a=0;$k=0;$da=$opt_j?$opt_j:1;

# Open header file if available
if(-e "$e/header") {
    open A,"$e/header" or die "Error reading header\n";
    $_=<A>;
    ($t_s,$t_e,$frn)=split;
    close A;
    $header=1;
}

# Set gnuplot header and read information about image cropping
open A,"gp_headers/t".($opt_g?$opt_g:1).".gnuplot" or die "Can't find gnuplot header";

# Read the rest of the file, skipping and altering lines as necessary
$gpn=0;
while(<A>) {

    # Header substitution
    if($header) {s/^HEA://;} else {next if /^HEA:/;}

    # Tracer substitution
    if($opt_t) {s/^TRA://g;} else {next if m/^TRA:/;}
    $gp[$gpn++]=$_;
}
close A;
$gpn--;

# Loop over the available frames
while(-e "$e/$ARGV[1].$a") {

    # Terminate if the specified frame limit has been reached
    last if defined $opt_n && $a>$opt_n;

    # Prepare input and output filenames
    $za=sprintf "_%04d",$a;
    $of="$odir\/fr$za";
    $infile="$e/$ARGV[1].$a";

    # Skip existing file if -d option is in use
    if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
        print "$a (skipped)\n";
        $a+=$da;
        next;
    }

    # Prepare tracer filenames
    $tfile="$e/trace.$a";

    # Create temporary Gnuplot file
    print "Frame $a (thread $P)\n";
    open B,">$odir/temp$P.gnuplot";
    foreach $i (0..$gpn) {
        $_=$gp[$i];

        # File substitutions for tracers
        s/TRACERS/$tfile/g if $opt_t;

        # File substitutions for timing information
        if($header) {
            $ti=$t_s+($t_e-$t_s)*$a/$frn;
            $tif=sprintf "%.2f",$ti;
            s/TIME/$tif/;
        }

        # General file substitutions
        s/CBRANGE/$cb/g;
        s/INFILE/$infile/g;
        s/OUTFILE/$of/g;
        print B;
    }
    close B;

    # Make a fork to create the graph
    exec "gnuplot $odir/temp$P.gnuplot 2>/dev/null" if (($pid[$P]=fork)==0);

    # Wait for one of the forked jobs to finish
    if ($queue) {
        $piddone=wait;$P=1;
        $P++ while $piddone!=$pid[$P] && $P<=$nodes;
        die "PID return error!\n" if $P>$nodes;
    } else {
        $P++;$queue=1 if $P>=$nodes;
    };
    $a+=$da;
}

# Wait for all the remaining forked jobs to finish
wait foreach 1..($queue?$nodes:$h-1);

# Additional code to automatically make a movie
unless ($opt_w) {
    $mf=$ebase."_".$ARGV[1];
    unlink "$mf.mov";
    system "ffmpeg -r 24 -y -i $odir/fr_%4d.png -preset slow -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart $mf.mov";
}

# Delete the temporary output directory
system "rm -rf $odir" unless $opt_k || $opt_w;
