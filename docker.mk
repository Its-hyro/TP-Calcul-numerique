#######################################
# docker.mk
# Options pour l'environnement Docker
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native 