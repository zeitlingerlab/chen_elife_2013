MEME=/data/apps/meme/bin/meme

$MEME meme/prembt_active.fasta -mod zoops -dna -nmotifs 25 -revcomp -maxw 12 -maxsize 5000000 -oc meme/prembt_active
$MEME meme/mbt_active.fasta    -mod zoops -dna -nmotifs 25 -revcomp -maxw 12 -maxsize 5000000 -oc meme/mbt_active
