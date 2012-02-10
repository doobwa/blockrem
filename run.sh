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
./brem.Rscript -f "data/synthetic.rdata" -k 2 -n 2000 -m "full" -d "results/synthetic"
