opts=list(dataset="irvine",numclusters=1,model.type="full",niter=10,gibbs="fast",numiterations=100,slice=FALSE)

# Fit models
./parallel --sshlogin 4/d5,4/d6,4/d7 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -m {3} -s FALSE' ::: "irvine" ::: 1 2 3 ::: "shared" "full"

# Make predictions
./parallel --sshlogin 8/d10 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -m {2}' ::: "synthetic" "eckmann-small" ::: "full.1" "full.2" "shared.2" "uniform" "online" "marg"

./parallel --sshlogin 2/d4,2/d5,2/d6,2/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -m {2}' ::: "irvine" ::: "full.2" "shared.2" "uniform" "online" "marg"


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
rsync -auvz data/synthetic.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
