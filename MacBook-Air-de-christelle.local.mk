#######################################
# darwin.mk
# Options pour MacOS avec Homebrew
#######################################
CC=gcc
LIBSLOCAL=-L/opt/homebrew/opt/lapack/lib -L/opt/homebrew/opt/openblas/lib -llapack -lopenblas -lm
INCLUDEBLASLOCAL=-I/opt/homebrew/opt/lapack/include -I/opt/homebrew/opt/openblas/include
OPTCLOCAL=-fPIC -march=native 