.PHONY: all gpu cuda clean test

all:
	$(MAKE) -C src all

# Portable OpenCL GPU build (Apple Silicon / AMD / Intel) -> bin/gramsci_cl
gpu:
	$(MAKE) -C src_opencl

# NVIDIA OpenACC GPU build (requires nvfortran) -> bin/gramsci_gpu
cuda:
	$(MAKE) -C src_gpu

test: all
	$(MAKE) -C tests test

clean:
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	@rm -f example/test.2pcf example/test.2pcf.lines example/test.2pcf.md5
