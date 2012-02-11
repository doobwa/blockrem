parallel ./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k {1} -n 1000 -m {2} -d data/eckmann ::: 2 3 5 ::: "baserates" "single" "shared" "full" &
seq 4 | parallel --sshlogin 8/d9,8/d10,8/d11 -j 3 'cd /extra/duboisc0/blockrem; ./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 1000 -m "shared" -d data/eckmann'
parallel --sshlogin 8/d9,8/d10,8/d11 -j 4 'cd /extra/duboisc0/blockrem; ./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k {1} -n 1000 -m {2} -d data/eckmann ::: 2 3 5 ::: "baserates" "single" "shared" "full" &'

./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 1 -n 2000 -m "full" -d "results/eckmann-small"

./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 2000 -m "baserates" -d "results/eckmann-small"
./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 2000 -m "shared" -d "results/eckmann-small"
./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 2000 -m "full" -d "results/eckmann-small"

./R/brem.Rscript -f "data/synthetic.rdata" -k 2 -n 2000 -m "baserates" -d "results/synthetic"
./R/brem.Rscript -f "data/synthetic.rdata" -k 2 -n 2000 -m "shared" -d "results/synthetic"
./R/brem.Rscript -f "data/synthetic.rdata" -k 1 -n 2000 -m "full" -d "results/synthetic"

# Use this now:
./brem.r -f "data/synthetic.rdata" -k 2 -n 1000 -m {} -d "results/synthetic" -g "slow" ::: "full" "slow"

seq 10|./parallel --sshlogin 8/d10 -j +0 'cd /extra/duboisc0/blockrem/; ./brem.r'
seq 5|./parallel --sshlogin 8/d9,8/d10 'module load R/2.13.1; cd /extra/duboisc0/blockrem/; ./tmp.r'
ls | parallel --sshlogin 8/d8,8/d9,8/d10 'echo -n {}" "; ls {}|wc -l'

./parallel --sshlogin 8/d8,8/d9,8/d10 'cd /extra/duboisc0/blockrem;./brem.r -f "data/synthetic.rdata" -k 2 -n 1000 -m {} -d "results/synthetic" -g "slow"'  ::: "full" "shared" "baserates"

cd ~/Documents/blockrem/
rsync -auvz . duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results blockrem/
