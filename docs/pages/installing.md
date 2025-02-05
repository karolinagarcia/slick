# installing

SLICK depends on versions of DESPOTIC and CAESAR which are not available on pypi, so we'll have to install those ourselves.
Download these as well as SLICK itself.
```
  git clone git@github.com:karolinagarcia/slick.git
  git clone git@bitbucket.org:krumholz/despotic.git
  git clone git@github.com:dnarayanan/caesar.git
```

Make a python environment for slick to exist in:
```
  conda create -y --name slick python=3.10.4
  conda activate slick
```

DESPOTIC depends on libgsl 2.7. 
This can be resolved on HiPerGator by loading the corresponding module:
```
  module load gcc/12.2.0 gsl/2.7
```

We're going to install SLICK first.
It may seem weird to do this before the dependencies, but doing so in this order allows pip to install the dependencies that *are* on pypi (numpy, yt, etc.) for us.
```
  cd slick
  pip install .
  cd ..
```

Install DESPOTIC.
Doing so requires a patch to the makefile which allows the compilers to know where gsl is located on hipergator::
```
  cd despotic
  git checkout 182cd46d
  curl -L https://gist.githubusercontent.com/smsutherland/f12e6dac5bc91c5a227ea349dcce9098/raw/ | git apply
  python setup.py install
  cd ..
```

Install CAESAR:
```
  cd caesar
  git checkout da0dba1e
  python setup.py install
```

If all has gone well, you should be able to run ``slick -h`` and get a help message.
