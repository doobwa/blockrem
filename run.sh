#opts=list(dataset="synthetic",numclusters=2,model.type="full",gibbs=TRUE,numiterations=100,slice=TRUE,mh=FALSE,initialize=FALSE,fixz=FALSE,skip.intercept=FALSE)
#opts=list(dataset="synthetic",numclusters=2,numiterations=20,splitmerge=FALSE,numextra=2,model.type="full")

./brem.r --dataset "classroom-16" --kinit 2 --kmax 2 --numiterations 500 --splitmerge FALSE --pshifts TRUE --degrees TRUE --negbinom FALSE --collapse TRUE --numextra 5 --force TRUE --sigma.hyper 0.001 --crp.hyper 1


./parallel --sshlogin 8/d5,8/d7,8/d8,4/d9,4/d10,8/d11,4/d12 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --kinit 2 --kmax {2} --numiterations 500 --splitmerge FALSE --pshifts {3} --degrees {4} --negbinom FALSE --collapse TRUE --numextra 5 --force TRUE --sigma.hyper.alpha {5} --sigma.hyper.beta {6} --crp.hyper 1' :::  "synthetic-1" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "enron-small" "realitymining-small" "twitter-small" ::: 1 2 3 10 ::: TRUE FALSE ::: TRUE FALSE ::: 5 1 ::: 1






./parallel --sshlogin 8/d5,8/d7,8/d8,4/d9,4/d10,8/d11,8/d12 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --kinit {2} --kmax {3} --numiterations 500 --splitmerge FALSE --pshifts {4} --degrees {5} --negbinom FALSE --collapse TRUE --numextra 5 --force FALSE --sigma.hyper 1 --crp.hyper 1' :::  "synthetic-1" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "enron-small" "realitymining-small" ::: 1 ::: 1 ::: TRUE FALSE ::: TRUE FALSE

./parallel --sshlogin 8/d12 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --kinit {2} --kmax {3} --numiterations 500 --splitmerge FALSE --pshifts {4} --degrees {5} --negbinom FALSE --collapse TRUE --numextra 5 --force FALSE --sigma.hyper 1 --crp.hyper 1' :::  "twitter-small" "irvine" ::: 1 ::: 1 ::: TRUE FALSE ::: TRUE FALSE


./parallel --sshlogin  8/d11,8/d12 'cd /extra/duboisc0/blockrem;./predict.r --dataset {1} --force {2} --baselines {3} --niters {4}' :::  "synthetic-1" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "enron-small" "realitymining-small" ::: TRUE ::: TRUE FALSE ::: 10

./parallel --sshlogin  8/d11,8/d12 'cd /extra/duboisc0/blockrem;./predict.r --dataset {1} --force {2} --baselines {3} --niters {4}' ::: "classroom-17"  ::: FALSE ::: TRUE FALSE ::: 10


"synthetic-1" "eckmann-small" "classroom-16" "classrom-17" "classroom-27" "classroom-29" "classroom-31" "realitymining-small" "twitter-small" "enron-small" "irvine"

# All
./parallel --sshlogin 8/d11 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --kinit {2} --numiterations 500 --splitmerge {3} --pshifts {4} --degrees {5} --negbinom {6} --collapse {7} --numextra {8} --force {9} --sigma.hyper {10} --crp.hyper {11}' ::: "synthetic-1" "eckmann-small" "classroom-16" "classrom-17" "classroom-27" ::: 2 ::: FALSE ::: FALSE  ::: FALSE  ::: FALSE ::: FALSE  ::: 5 ::: FALSE ::: 
./parallel --sshlogin 8/d11 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --numclusters {2} --numiterations 500 --splitmerge {3} --pshifts {4} --degrees {5} --negbinom {6} --collapse {7} --numextra 5 --force {8}' ::: "realitymining-small" ::: 1 10 ::: FALSE ::: TRUE ::: FALSE TRUE ::: FALSE ::: FALSE  ::: FALSE

# Intercept only models with K=10
./parallel --sshlogin 6/d5,4/d6 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --numclusters {2} --numiterations 500 --splitmerge {3} --pshifts {4} --degrees {5} --negbinom {6} --collapse {7} --numextra 5 --force {8}' :::  "synthetic-1" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "enron-small" "twitter-small" "irvine" ::: 10 ::: FALSE ::: FALSE ::: FALSE ::: FALSE ::: FALSE ::: TRUE

# Intercept and pshift with K=10
./parallel --sshlogin 6/d8,4/d10,6/d11 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --numclusters {2} --numiterations 500 --splitmerge {3} --pshifts {4} --degrees {5} --negbinom {6} --collapse {7} --numextra 5 --force {8}' :::  "synthetic-1" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "enron-small" "twitter-small" ::: 1 10 ::: FALSE ::: TRUE ::: FALSE ::: FALSE ::: FALSE ::: TRUE

# Intercept and pshift and degree with K=1 and K=10
./parallel --sshlogin  6/d8,4/d10,6/d11 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --numclusters {2} --numiterations 500 --splitmerge {3} --pshifts {4} --degrees {5} --negbinom {6} --collapse {7} --numextra 5 --force {8} --sigmahyper {9}' ::: "synthetic-1" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "enron-small" ::: 1 10 ::: FALSE ::: TRUE ::: TRUE ::: FALSE ::: FALSE ::: TRUE ::: .1 .01



./parallel --sshlogin 8/d8 'cd /extra/duboisc0/blockrem;./brem.r --dataset {1} --numclusters {2} --numiterations 500 --splitmerge FALSE --pshifts TRUE --degrees {3} --negbinom {4} --collapse {5} --numextra 5' ::: "eckmann-small" ::: 1 10 ::: FALSE TRUE ::: FALSE ::: FALSE TRUE

./parallel --sshlogin 4/d4,4/d6,4/d7,4/d11,8/d12 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -s {3} -g {4} -b {5} -e 5' ::: "classroom-16" "classrom-17" "classroom-27" "classroom-29" "classroom-31" "enron-small" "irvine" ::: 10 ::: FALSE ::: FALSE ::: FALSE TRUE

./parallel --sshlogin 4/d4,4/d6,4/d7,4/d11,8/d12 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -s {3} -g {4} -b {5} -e 5' ::: "synthetic-1" "eckmann-small" "twitter-small" "realitymining-small" "classroom-16" "classroom-17" "classroom-27" "classroom-29" "classroom-31" "enron-small" "irvine" ::: 1 ::: FALSE ::: FALSE ::: FALSE TRUE

./parallel --sshlogin 8/d11,8/d12 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -s {3} -g {4} -b {5} -e 5' ::: "realitymining-small" "enron-small" ::: 10 ::: FALSE ::: FALSE TRUE ::: FALSE TRUE

./parallel --sshlogin 4/d11 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -s {3} -g {4} -b {5} -e 5' ::: "synthetic-1" ::: 10 ::: FALSE ::: FALSE ::: FALSE TRUE

 "realitymining-small" "eckmann-small" "classroom-16" "classroom-17" "classroom-27" "classroom-29" "classroom-31" "enron-small"

./parallel --sshlogin  4/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1}' ::: "realitymining-small" "enron-small" "twitter-small" "irvine"

./parallel --sshlogin  4/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1}' ::: "realitymining-small" "enron-small" "twitter-small" "irvine"

./parallel --sshlogin 4/d7 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 300 -s TRUE' ::: "eckmann-small" "twitter-small" "realitymining-small" ::: 20

./parallel --sshlogin 6/d5,6/m 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -t {3} -s TRUE' ::: "eckmann-small" ::: 3 ::: "full"

./parallel --sshlogin 6/d5,6/m 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -t {3} -s TRUE --initialize TRUE' ::: "eckmann-small" ::: 2 3 ::: "full"

./parallel --sshlogin 6/d12,6/m 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -t {3} -s TRUE' ::: "synthetic"  ::: 1 2 :::  "full"

./parallel --sshlogin 3/d7,3/d8 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 200 -t {3} -s TRUE' ::: "twitter-small" ::: 2 3 4 ::: "shared" "full"

./parallel --sshlogin 6/d12 'cd /extra/duboisc0/blockrem;./brem.r -d {1} -k {2} -n 500 -t {3} -s TRUE' ::: "synthetic" "eckmann-small" "twitter-small" ::: 2 3 ::: "baserates"


# Get counts
./parallel --sshlogin 3/d2,3/d3 'cd /extra/duboisc0/blockrem;./getcounts.r -d {1}' ::: "synthetic" 
"twitter-small" "eckmann-small"
 "twitter-small" "irvine"

# Make predictions: for twitter, only 1 can put on each datalab

./parallel './predict.r -d {1} -t {2}' ::: "synthetic-1" ::: "full" "uniform" "online" "marg"



./parallel --sshlogin 2/d4,2/d5,2/d6,2/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "synthetic" ::: "full.2" "full.1" "uniform" "online" "marg" "truth"

./parallel --sshlogin 2/d4,2/d5,2/d6,2/d7 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "eckmann-small" ::: "full.3" "full.2" "full.1" "uniform" "online" "marg"

./parallel --sshlogin 8/d5 'cd /extra/duboisc0/blockrem;./predict.r -d {1} -t {2}' ::: "synthetic" "eckmann-small" "twitter-small" ::: "baserates.2" "baserates.3"

"full.1" "full.2" "shared.2" "uniform" "online" "marg"


# Run dashboard
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
./parallel --sshlogin 8/d6 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE --predictions TRUE' ::: "synthetic-1"
./parallel 'cd /extra/duboisc0/blockrem;./dashboard.r -d {} -s TRUE --predictions TRUE' ::: "synthetic-1"

# Helper commands
cd ~/Documents/blockrem/
rsync -auvz pkg duboisc@d1:/extra/duboisc0/blockrem/ --exclude '*.so' --exclude '*.o'
rsync -auvz pkg/R/utils.r duboisc@d1:/extra/duboisc0/blockrem/pkg/R/
rsync -auvz brem.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz predict.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz getcounts.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz data/synthetic-1.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/twitter-small.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/eckmann-small.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/realitymining-small.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/enron-small.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/irvine.rdata duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz data/classroom-* duboisc@d1:/extra/duboisc0/blockrem/data/
rsync -auvz dashboard.r duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz tmp.rdata duboisc@d1:/extra/duboisc0/blockrem/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results . --exclude 'llks' --exclude 'ranks' 
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results/synthetic-1/final results/synethic-1/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results/twitter-small/final results/twitter-small/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results/synthetic-1/final results/synethic-1/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/results/synthetic-1/final results/synethic-1/
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/figs .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/data .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.rdata .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.slice.rdata .
rsync -auvz duboisc@d1:/extra/duboisc0/blockrem/tmp.slice.m5.rdata .

rm -rf results/synthetic-1/llks
rm -rf results/synthetic-1/ranks
rm -rf results/eckmann-small/llks
rm -rf results/eckmann-small/ranks
rm -rf results/enron-small/llks
rm -rf results/enron-small/ranks
rm -rf results/classroom-16/llks
rm -rf results/classroom-16/ranks
rm -rf results/classroom-17/llks
rm -rf results/classroom-17/ranks
rm -rf results/classroom-27/llks
rm -rf results/classroom-27/ranks
rm -rf results/reality-mining/llks
rm -rf results/reality-mining/ranks
rm -rf results/twitter-small/llks
rm -rf results/twitter-small/ranks
