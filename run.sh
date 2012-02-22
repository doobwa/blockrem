opts=list(dataset="synthetic",numclusters=1,model.type="full",niter=10,gibbs="fast",numiterations=100,slice=TRUE)

# Fit models
./parallel --sshlogin 4/d6 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -m {3} -s TRUE' ::: "irvine" ::: 2 4 ::: "shared" "full"

./parallel --sshlogin 6/d5,6/m 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -m {3} -s TRUE' ::: "synthetic" "eckmann-small" ::: 1 2 3 ::: "shared" "full"


./parallel --sshlogin 3/d7,3/d8 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -m {3} -s TRUE' ::: "twitter-small" ::: 2 3 ::: "shared" "full"

# Get counts
./parallel --sshlogin 3/d2,3/d3 'cd /extra/duboisc0/blockrem;./getcounts.r -d {1} -m {2}' ::: "eckmann-small" ::: "full.1" "full.2" "shared.2" "uniform" "online" "marg"

 "twitter-small" "irvine"

# Make predictions: for twitter, only 1 can put on each datalab

./parallel --sshlogin  7/d10 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -m {2}' ::: "eckmann-small" ::: "full.1" "full.2" "full.3" "shared.2" "shared.3" "uniform" "online" "marg"

./parallel --sshlogin 2/d4,2/d5,2/d6,2/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -m {2}' ::: "irvine" ::: "full.2" "shared.2" "uniform" "online" "marg"

./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -m {2}' ::: "eckmann-small" ::: "full.1" "full.2" "shared.2" "uniform" "online" "marg"

./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -m {2} -i TRUE --fixz={3} -s {4}' ::: "twitter-small" ::: "full.2" ::: TRUE FALSE ::: TRUE

./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -m {2} -i TRUE --fixz={3} -s {4}' ::: "twitter-small" ::: "full.2" ::: TRUE FALSE ::: FALSE


# Run dashboard
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "synthetic" 
./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "twitter-small" 
./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "synthetic" "eckmann-small" "twitter-small"

# Helper commands
cd ~/Documents/blockrem/
rsync -auvz pkg duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz brem.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz predict.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz getcounts.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz data/synthetic.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/data .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.rdata .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.slice.rdata .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.slice.m5.rdata .

