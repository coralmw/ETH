
all:
	icc -mmic -mkl -static dysv.c -o dysv
	scp dysv mic0@/home/tparks/


spins: spins.c kronecker_product.c tools.h
	icc -Wall -O3 -mmic -mkl -static -std=c99 -qopenmp spins.c kronecker_product.c -o spins

spinsmac: spins.c kronecker_product.c tools.h
	gcc spins.c kronecker_product.c  -lmkl_intel -lmkl_core -lmkl_intel_thread -o spins


energy: energybasisdensity.c kronecker_product.c tools.h
	icc -Wall -O3 -mmic -mkl -static -std=c99 -qopenmp energybasisdensity.c kronecker_product.c -o energy

upload: energy
	scp energy mic0:/home/tparks/

energylocal: energybasisdensity.c kronecker_product.c tools.h
	icc -std=c11 energybasisdensity.c kronecker_product.c -lmkl_core -lmkl_intel_thread -qopenmp -lpthread -lmkl_intel_lp64 -g  -o energylocal

eigentest: spins.c kronecker_product.c tools.h
	icc -Wall -O3 -mmic -mkl -static -std=c99 -qopenmp test_eigen.c kronecker_product.c -o eigen
	scp eigen mic0:/home/tparks/
