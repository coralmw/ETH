
all:
	icc -mmic -mkl -static dysv.c -o dysv
	scp dysv mic0@/home/tparks/


spins: spins.c kronecker_product.c tools.h
	icc -Wall -O3 -mmic -mkl -static -std=c99 -qopenmp spins.c kronecker_product.c -o spins

upload: spins
	scp spins mic0:/home/tparks/
