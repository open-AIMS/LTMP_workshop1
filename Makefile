.PHONY: build build_singularity ssh_singularity

build:
	docker build . --tag ltmp

build_singularity:
	docker save ltmp -o ltmp.tar 
	singularity build ltmp.sif docker-archive://ltmp.tar

ssh_singularity:
	scp ltmp.sif mlogan@hpc-l001.aims.gov.au:~/Work/AIMS/LTMP/LTMP-Workshop1
