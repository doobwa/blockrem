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

opts=list(dataset="eckmann-small",numclusters=1,model.type="full",niter=10,gibbs="fast",numiterations=100,slice=TRUE)

./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./brem.r -d {} -k 1 -n 500 -m "full"' ::: "synthetic" "eckmann-small" "twitter-small"

./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./brem.r -d {} -k 2 -n 500 -m "shared"' ::: "synthetic" "eckmann-small" 

# Fit models
cd /extra/duboisc0/blockrem;./brem.r -d "twitter-small" -k 1 -n 500 -m "full" 
./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./brem.r -d "synthetic" -k 2 -n 500 -m {} -s FALSE' ::: "baserates" "shared" "full"

./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./brem.r -d "eckmann-small" -k 2 -n 500 -m {}' ::: "baserates" "shared" "full"

./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./brem.r -d "twitter-small" -k 2 -n 500 -m {}' ::: "baserates" "shared" "full"

# Run dashboard
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "synthetic" 
./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "synthetic" "eckmann-small" "twitter-small"


rsync -auvz . duboisc@d1:/extra/duboisc0/blockrem/

cd ~/Documents/blockrem/
rsync -auvz pkg duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz brem.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz data duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
