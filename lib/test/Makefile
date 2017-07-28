all:
	mkdir -p build/
	cd build/ && cmake ..
	$(MAKE) -C build

test:
	$(MAKE) -C build test

clean:
	rm -rf build/ ./CMakeFiles

.PHONY: all clean test
