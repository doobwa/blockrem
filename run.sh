#opts=list(dataset="synthetic",numclusters=2,model.type="full",gibbs=TRUE,numiterations=100,slice=TRUE,mh=FALSE,initialize=FALSE,fixz=FALSE,skip.intercept=FALSE)

# Fit models
./parallel --sshlogin 4/d6 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -t {3} -s TRUE' ::: "irvine" ::: 2 4 ::: "shared" "full"

./parallel --sshlogin 6/d5,6/m 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -t {3} -s TRUE' ::: "eckmann-small" ::: 1 2 3 ::: "full"

./parallel --sshlogin 6/d12,6/m 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -t {3} -s TRUE' ::: "synthetic"  ::: 1 2 :::  "full"


./parallel --sshlogin 3/d7,3/d8 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -t {3} -s TRUE' ::: "twitter-small" ::: 2 3 4 ::: "shared" "full"

./parallel --sshlogin 6/d12 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -t {3} -s TRUE' ::: "synthetic" "eckmann-small" "twitter-small" ::: 2 3 ::: "baserates"


# Get counts
./parallel --sshlogin 3/d2,3/d3 'cd /extra/duboisc0/blockrem;./getcounts.r -d {1}' ::: "synthetic" 
"twitter-small" "eckmann-small"
 "twitter-small" "irvine"

# Make predictions: for twitter, only 1 can put on each datalab

./parallel --sshlogin  7/m 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "twitter-small" ::: "full.1" "full.2" "online" "full.3" "shared.2" "shared.3" "uniform" "online" "marg"

./parallel --sshlogin 2/d4,2/d5,2/d6,2/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "synthetic" ::: "full.2" "full.1" "uniform" "online" "marg"

./parallel --sshlogin 2/d4,2/d5,2/d6,2/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "eckmann-small" ::: "full.3" "full.2" "full.1" "uniform" "online" "marg"

./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "synthetic" "eckmann-small" "twitter-small" ::: "baserates.2" "baserates.3"

"full.1" "full.2" "shared.2" "uniform" "online" "marg"


# Run dashboard
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "synthetic" 
./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE' ::: "twitter-small" 
./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE --predictions TRUE' ::: "eckmann-small"

# Helper commands
cd ~/Documents/blockrem/
rsync -auvz pkg duboisc@d1:/extra/duboisc0/blockrem/ --exclude '*.so' --exclude '*.o'
rsync -auvz brem.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz predict.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz getcounts.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz data/synthetic.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/eckmann-small.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz tmp.rdata duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/data .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.rdata .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.slice.rdata .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.slice.m5.rdata .

