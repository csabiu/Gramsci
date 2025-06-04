.PHONY: all regtest clean

all:
	$(MAKE) -C src all

regtest: all
	./regtest.sh

clean:
	$(MAKE) -C src clean
	@rm -f example/test.2pcf example/test.2pcf.lines example/test.2pcf.md5 \
	example/1.loadnodes.lines example/1.loadnodes.md5
