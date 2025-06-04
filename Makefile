.PHONY: all regtest clean test

all:
	$(MAKE) -C src all
test: all
	$(MAKE) -C tests test


regtest: all
	./regtest.sh

clean:
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	@rm -f example/test.2pcf example/test.2pcf.lines example/test.2pcf.md5 \
	example/1.loadnodes.lines example/1.loadnodes.md5
