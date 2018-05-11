all:
	mkdir -p build/
	cd build/ && cmake .. -DSOPHUS_INCLUDE_DIR=${SOPHUS_ROOT}
	$(MAKE) -C build

debug:
	mkdir -p build/
	cd build/ && cmake .. -DCMAKE_BUILD_TYPE=Debug -DSOPHUS_INCLUDE_DIR=${SOPHUS_ROOT}
	$(MAKE) -C build

test:
	$(MAKE) -C build test

clean:
	rm -rf build/ ./CMakeFiles

.PHONY: all clean test
