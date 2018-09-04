$in=$ARGV[0];
open(IN,"<$in");
$first=1;
while (<IN>) {
chomp $_;
@W=split(/[\s]+/,$_);
if (! $first) {
if ($W[0] =~ /$prevqid/ && $prevqid =~ /$W[0]/) {
if (($W[1] > $prevqs && $W[2] < $prevqe) || ($W[1] > $prevqs && $W[2]
== $prevqe) || ($W[1] == $prevqs && $W[2] < $prevqe)) {
print "\n$_\tNested";
}
elsif (($prevqs > $W[1] && $prevqe < $W[2]) || ($prevqs > $W[1] &&
$W[2] == $prevqe) || ($W[1] == $prevqs && $prevqe < $W[2])) {
 print "\tNested\n$_";
 }
else {
print "\n$_";
}
}
else {
print "\n$_";
}
}
else {
print "\n$_";
}
$prevqid=$W[0];
$prevqs=$W[1];
$prevqe=$W[2];
$first=0;
}

