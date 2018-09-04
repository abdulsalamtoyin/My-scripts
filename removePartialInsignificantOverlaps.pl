$in=$ARGV[0];
open(IN,"<$in");
$first=1;
while (<IN>) {
chomp $_;
@W= split(/[\s]+/,$_);
if (! $first) {
if ($W[0] =~ /$prevqid/ && $prevqid =~ /$W[0]/) {
if ($W[1] > $prevqs && $W[2] > $prevqe) { # Partial
if (($W[1] - $prevqs) < 10 && ($W[2]-$prevqe) > 10) {
print "\tPartial\n$_";
}
elsif (($W[1] - $prevqs) > 10 && ($W[2] - $prevqe) < 10) {
print "\n$_\tPartial";
}
elsif ((($W[1] - $prevqs) <= 10 && ($W[2] - $prevqe) <= 10)
||(($W[1] - $prevqs) >= 10 && ($W[2] - $prevqe) >= 10)) {
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
}
else {
print "\n$_";
}
$prevqid=$W[0];
$prevqs=$W[1];
$prevqe=$W[2];
$first=0;
}
print "\n";
