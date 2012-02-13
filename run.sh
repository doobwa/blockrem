parallel ./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k {1} -n 1000 -m {2} -d data/eckmann ::: 2 3 5 ::: "baserates" "single" "shared" "full" &
seq 4 | parallel --sshlogin 8/d9,8/d10,8/d11 -j 3 'cd /extra/duboisc0/blockrem; ./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 1000 -m "shared" -d data/eckmann'
seq 3 | parallel --sshlogin 8/d10,8/d11 './tmp.Rscript'
 -f "data/eckmann.dyadic.small.rdata" -k {1} -n 1000 -m {2} -d data/eckmann ::: 2 3 5 ::: "baserates" "single" "shared" "full" &'

./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 1 -n 2000 -m "full" -d "results/eckmann-small"

./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 2000 -m "baserates" -d "results/eckmann-small"
./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 2000 -m "shared" -d "results/eckmann-small"
./R/brem.Rscript -f "data/eckmann.dyadic.small.rdata" -k 2 -n 2000 -m "full" -d "results/eckmann-small"

./R/brem.Rscript -f "data/synthetic.rdata" -k 2 -n 2000 -m "baserates" -d "results/synthetic"
./R/brem.Rscript -f "data/synthetic.rdata" -k 2 -n 2000 -m "shared" -d "results/synthetic"
./R/brem.Rscript -f "data/synthetic.rdata" -k 1 -n 2000 -m "full" -d "results/synthetic"

# Use this now:
./brem.r -f "data/synthetic.rdata" -k 2 -n 500 -m "baserates" -d "results/synthetic" -g "fast"
./brem.r -f "data/synthetic.rdata" -k 2 -n 500 -m "shared" -d "results/synthetic" -g "fast"
./brem.r -f "data/synthetic.rdata" -k 2 -n 500 -m "full" -d "results/synthetic" -g "fast"
./brem.r -f "data/synthetic.rdata" -k 1 -n 500 -m "full" -d "results/synthetic" -g "fast"

./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m "baserates" -d "results/eckmann-small" -g "fast"
./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m "shared" -d "results/eckmann-small" -g "fast"
./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m "full" -d "results/eckmann-small" -g "fast"
./brem.r -f "data/eckmann-small.rdata" -k 1 -n 500 -m "full" -d "results/eckmann-small" -g "fast"

./parallel --sshlogin 8/d10,8/d11 'cd /extra/duboisc0/blockrem;./tmp.r {}' ::: 1 2 3
./parallel --sshlogin 8/d10,8/d11 'cd /extra/duboisc0/blockrem;./brem.r -f \"data/eckmann-small.rdata\" -k 5 -n 500 -m {} -d \"results/eckmann-small\"' ::: "shared" "full"
./parallel --sshlogin 8/d10,8/d'cd /extra/duboisc0/blockrem;./brem.r -f \"data/eckmann-small.rdata\" -k {} -n 500 -m {} -d \"results/eckmann-small\"' ::: 3 5 ::: "shared" "full"


./dashboard.r -d "synthetic"

cd ~/Documents/blockrem/
rsync -auvz . duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
