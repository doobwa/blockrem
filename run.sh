# Synthetic
./brem.r -f "data/synthetic.rdata" -k 2 -n 500 -m "baserates" -d "results/synthetic" -g "fast"
./brem.r -f "data/synthetic.rdata" -k 2 -n 500 -m "shared" -d "results/synthetic" -g "fast"
./brem.r -f "data/synthetic.rdata" -k 2 -n 500 -m "full" -d "results/synthetic" -g "fast"
./brem.r -f "data/synthetic.rdata" -k 1 -n 500 -m "full" -d "results/synthetic" -g "fast"

# Eckmann-small
./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m "baserates" -d "results/eckmann-small" -g "fast"
./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m "shared" -d "results/eckmann-small" -g "fast"
./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m "full" -d "results/eckmann-small" -g "fast"
./brem.r -f "data/eckmann-small.rdata" -k 1 -n 500 -m "full" -d "results/eckmann-small" -g "fast"

opts=list(dataset="synthetic",numclusters=5,model.type="full",niter=10,gibbs=TRUE,numiterations=10)

./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./brem.r -d {} -k 1 -n 500 -m "full"' ::: "synthetic" "eckmann-small" "twitter"

./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./brem.r -d {} -k 2 -n 500 -m "full"' ::: "synthetic" "eckmann-small" 

./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./brem.r -f "data/eckmann-small.rdata" -k 2 -n 500 -m {} -d "results/eckmann-small"' ::: "baserates" "shared" "full"

./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./dashboard.r -d {}' ::: "synthetic" "eckmann-small"


rsync -auvz . duboisc@d1:/extra/duboisc0/blockrem/

cd ~/Documents/blockrem/
rsync -auvz pkg duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz brem.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
