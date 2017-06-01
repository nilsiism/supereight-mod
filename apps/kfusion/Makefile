BLOCK_SIZE = 8
SIZE = 512
C=2
#### ICL-NUIM TRAJECTORIES PARAMETERS  ####
0 = -s 5.0 -p 0.34,0.5,0.24 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F
1 = -s 5.0 -p 0.485,0.5,0.55 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F  
2 = -s 5.0 -p 0.34,0.5,0.24 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F
3 = -s 5.0 -p 0.2685,0.5,0.4 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE)  -m 0.1 -F 
live = -s 5.0 -p 0.5,0.5,0.5 -z 1 -c 1 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1
fr2desk = -s 8.0 -p 0.5,0.5,0.5 -z 1 -c 2 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F
# fr2xyz = 	-s 8.0 -p 0.5,0.5,0.5	-z 1 -c 2 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE)  -m 0.1 -F 

ROOT_DIR=$(shell pwd)
TOON_DIR=${ROOT_DIR}/TooN/install_dir
TOON_INCLUDE_DIR=${TOON_DIR}/include/


all : TooN
	mkdir -p build/
	cd build/ && cmake -DCMAKE_BUILD_TYPE=Release -DTOON_INCLUDE_PATH=${TOON_INCLUDE_DIR} $(CMAKE_ARGUMENTS) ..
	$(MAKE) -C build  $(MFLAGS) $(SPECIFIC_TARGET)


debug: TooN
	mkdir -p build/
	mkdir -p build/logs/
	cd build/ && cmake -DCMAKE_BUILD_TYPE=Debug -DTOON_INCLUDE_PATH=${TOON_INCLUDE_DIR} $(CMAKE_ARGUMENTS) ..
	$(MAKE) -C build $(MFLAGS)

stats: 
	mkdir -p build/
	mkdir -p build/logs/
	cd build/ && cmake -DSTATS=ON ..
	$(MAKE) -C build $(MFLAGS)

TooN:
	git clone https://github.com/edrosten/TooN.git
	cd TooN &&  git checkout 92241416d2a4874fd2334e08a5d417dfea6a1a3f
	mkdir -p ${TOON_DIR}
	cd TooN && ./configure --prefix=${TOON_DIR} && make install

#### DATA SET GENERATION ####

living_room_traj%_loop.raw : living_room_traj%_loop
	if test -x ./build/thirdparty/scene2raw ; then echo "..." ; else echo "do make before"; false ; fi
	./build/thirdparty/scene2raw living_room_traj$(*F)_loop living_room_traj$(*F)_loop.raw

living_room_traj%_loop : 
	mkdir $@
	cd $@ ; wget http://www.doc.ic.ac.uk/~ahanda/$@.tgz; tar xzf $@.tgz 

livingRoom%.gt.freiburg : 
	echo  "Download ground truth trajectory..."
	if test -x $@ ; then echo "Done" ; else wget http://www.doc.ic.ac.uk/~ahanda/VaFRIC/$@ ; fi


#### LOG GENERATION ####

%.cpp.log  :  living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/cpp/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	KERNEL_TIMINGS=1 ./build/kfusion-benchmark-cpp $($(*F)) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$@ 2> $(LOGS)/kernels.$@
	./thirdparty/checkPos.py $(LOGS)/benchmark.$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$@
	./thirdparty/checkKernels.py $(LOGS)/kernels.$@   >> $(LOGS)/resume.$@

%.openmp.log  :  living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/openmp/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	KERNEL_TIMINGS=1  ./build/kfusion-benchmark-openmp $($(*F)) -d $(LOGS)/$(*F).openmp.$(SIZE) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$(SIZE).$@ 2> $(LOGS)/kernels.$(SIZE).$@
	./thirdparty/checkPos.py $(LOGS)/benchmark.$(SIZE).$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$(SIZE).$@
	./thirdparty/checkKernels.py $(LOGS)/kernels.$(SIZE).$@   >> $(LOGS)/resume.$(SIZE).$@

fr2desk.openmp.log : fr2_desk.raw 
	KERNEL_TIMINGS=1 ./build/kfusion-benchmark-octree-openmp $(fr2desk) -d octree-openmp-$(BLOCK_SIZE).$(SIZE) -i  fr2_desk.raw -o  benchmark.$(BLOCK_SIZE).$(SIZE).$@ 2> kernels.$(BLOCK_SIZE).$(SIZE).$@
	./thirdparty/checkPos.py benchmark.$(BLOCK_SIZE).$(SIZE).$@  livingRoom2.gt.freiburg > resume.$(BLOCK_SIZE).$(SIZE).$@
	./thirdparty/checkKernels.py kernels.$(BLOCK_SIZE).$(SIZE).$@   >> resume.$(BLOCK_SIZE).$(SIZE).$@

live.log : 
	./build/kfusion-qt-openmp $(live)


#### GENERAL GENERATION ####

clean :
	rm -rf build
cleanall : 
	rm -rf build
	rm -rf living_room_traj*_loop livingRoom*.gt.freiburg living_room_traj*_loop.raw
	rm -f *.log 


.PHONY : clean bench test all validate

.PRECIOUS: living_room_traj%_loop livingRoom%.gt.freiburg living_room_traj%_loop.raw

