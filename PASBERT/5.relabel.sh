tsv=${1}
echo processing $tsv
grep -v "^sequence" $tsv | awk '{seqs=""; for(i = 1; i < 196; i++){seqs=seqs""substr($i,0,1)} print ">"NR"\n"seqs$196}' > $tsv.fa
perl search.motif.pl $tsv.fa $tsv.pas.motif.pos 0 104
Rscript relabel.PAS.r $tsv
