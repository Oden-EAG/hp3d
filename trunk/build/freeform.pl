#!/usr/bin/perl -w

use strict;
use Pod::Usage;
use Getopt::Long;

my $opt_b;
my ($opt_help, $opt_man);
GetOptions("help", "man", "b")
    or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

my $comment_char = '!';
my $cont_char = '&';
my $old = "";
my $comment = "";
my $in_string = 1;
#
# main loop
#
while (<>) {
    my $new = $_;
#
# Delete trailing blanks and tabs
#
    $new =~ s/[ \t]*$//;
#
# Save comments, converting to '!' comments.  Note that "comments"
# include C preprocessor lines and blank lines.
#
    if ($new =~ /^[*c#!]|^$/i) {
if ($new =~ /^[*c]/i) {
    substr($new,0,1) = $comment_char;
}
$comment .= $new;
next;
    }
#
# Replace tabs with spaces
#
    $new =~ s/\t/        /g;
#
# Look for continues, make sure continuation is '&' if backward
# compatible or in a string
#
    if (substr($new,5,1) ne " " and substr($new,5,1)) {
        if ($opt_b || $in_string) {
            substr($new,5,1) = $cont_char;
        } else {
            substr($new,5,1) = ' ';
        }

#
# Check for ! comments in previous line, put & in column 73 and before !
#
        my($pos, $pad, $len);
        if ( ($pos = index($old,"!")) >= $[ ) {
            if ($opt_b && $pos < 72) {
                $pad = 72 - $pos;
            } else {
                $pad = 1;
            }
            substr($old,$pos,0) = ' ' x $pad . $cont_char;
        } else {
            $len = length($old);     # includes a '\n'
            if ($opt_b && $len < 73) {
                $pad = 73 - length($old);
            } else {
                $pad = 1;
            }
            substr($old,$len-1,0) = ' ' x $pad . $cont_char;
        }
    }
#
# Print $old and any "comments"
#
    print $old;
    print $comment;
    $comment = "";
    $old = $new;
}
#
# Print the last $old and "comments"
#
print $old;
print $comment;
__END__

freeform - Fortran 90 reformatter

=head1 SYNOPSIS

freeform [--help] [--man] [-b] [file ...] > outfile

=head1 DESCRIPTION

Convert to F90 free format - start comments with '!' and put '&'
at end of line for continuations.

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<-help>

Print more details about the arguments.

=item I<-man>

Print a full man page.

=item I<-b>

Maintain compatibility with fixed format.

=back

=head1 BUGS

The $in_string flag is always true because I'm not sure
how to check for it.

=head1 AUTHOR

Kate Hedstrom
kate@arsc.edu

=cut
