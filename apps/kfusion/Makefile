BLOCK_SIZE = 8
SIZE = 512
C=2
#### ICL-NUIM TRAJECTORIES PARAMETERS  ####
0 = -s 5.0 -p 0.34,0.5,0.24 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE) -m 0.1 -F
1 = -s 5.0 -p 0.485,0.5,0.55 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE) -m 0.1 -F  
2 = -s 5.0 -p 0.34,0.5,0.24 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE) -m 0.1 -F
3 = -s 5.0 -p 0.2685,0.5,0.4 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v$(SIZE),$(SIZE),$(SIZE)  -m 0.1 -F 
live = -s 5.0 -p 0.5,0.5,0.5 -z 1 -c 1 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE) -m 0.1
fr2desk = -s 8.0 -p 0.5,0.5,0.5 -z 1 -c 2 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE) -m 0.1 -F
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
	if test -x ./build/kfusion/thirdparty/scene2raw ; then echo "..." ; else echo "do make before"; false ; fi
	./build/kfusion/thirdparty/scene2raw living_room_traj$(*F)_loop living_room_traj$(*F)_loop.raw

living_room_traj%_loop : 
	mkdir $@
	cd $@ ; wget http://www.doc.ic.ac.uk/~ahanda/$@.tgz; tar xzf $@.tgz 

livingRoom%.gt.freiburg : 
	echo  "Download ground truth trajectory..."
	if test -x $@ ; then echo "Done" ; else wget http://www.doc.ic.ac.uk/~ahanda/VaFRIC/$@ ; fi


#### LOG GENERATION ####

%.opencl.log  : living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/opencl/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	LD_PRELOAD=./build/kfusion/thirdparty/liboclwrapper.so ./build/kfusion/kfusion-benchmark-opencl $($(*F)) -i  living_room_traj$(*F)_loop.raw -o $(LOGS)/benchmark.$@ 2> $(LOGS)/oclwrapper.$@
	cat  $(LOGS)/oclwrapper.$@ |grep -E ".+ [0-9]+ [0-9]+ [0-9]+" |cut -d" " -f1,4 >   $(LOGS)/kernels.$@
	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$@
	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$@   >> $(LOGS)/resume.$@

%.cpp.log  :  living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/cpp/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	KERNEL_TIMINGS=1 ./build/kfusion/kfusion-benchmark-cpp $($(*F)) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$@ 2> $(LOGS)/kernels.$@
	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$@
	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$@   >> $(LOGS)/resume.$@

%.octree-cpp.log  :  living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/octree-cpp/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	KERNEL_TIMINGS=1  ./build/kfusion/kfusion-benchmark-octree-cpp $($(*F)) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$(BLOCK_SIZE).$(SIZE).$@ 2> $(LOGS)/kernels.$(BLOCK_SIZE).$(SIZE).$@
	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$(BLOCK_SIZE).$(SIZE).$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$(BLOCK_SIZE).$(SIZE).$@
	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$(BLOCK_SIZE).$(SIZE).$@   >> $(LOGS)/resume.$(BLOCK_SIZE).$(SIZE).$@

%.octree-openmp.log  :  living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/octree-openmp/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	KERNEL_TIMINGS=1  ./build/kfusion/kfusion-benchmark-octree-openmp $($(*F)) -d $(LOGS)/$(*F).octree-openmp-$(BLOCK_SIZE).$(SIZE) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$(BLOCK_SIZE).$(SIZE).$@ 2> $(LOGS)/kernels.$(BLOCK_SIZE).$(SIZE).$@
	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$(BLOCK_SIZE).$(SIZE).$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$(BLOCK_SIZE).$(SIZE).$@
	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$(BLOCK_SIZE).$(SIZE).$@   >> $(LOGS)/resume.$(BLOCK_SIZE).$(SIZE).$@


%.openmp.log  :  living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/openmp/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	KERNEL_TIMINGS=1  ./build/kfusion/kfusion-benchmark-openmp $($(*F)) -d $(LOGS)/$(*F).openmp.$(SIZE) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$(SIZE).$@ 2> $(LOGS)/kernels.$(SIZE).$@
	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$(SIZE).$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$(SIZE).$@
	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$(SIZE).$@   >> $(LOGS)/resume.$(SIZE).$@

%.cuda.log  : living_room_traj%_loop.raw livingRoom%.gt.freiburg
	$(eval LOGS = logs/noiseless/cuda/$(*F)/$(SIZE)/$(shell date +%F--%T))
	mkdir -p $(LOGS)
	/usr/local/cuda-7.5/bin/nvprof --print-gpu-trace ./build/kfusion/kfusion-benchmark-cuda $($(*F)) -i  living_room_traj$(*F)_loop.raw -o  $(LOGS)/benchmark.$@ 2> $(LOGS)/nvprof.$@ || true
	cat  $(LOGS)/nvprof.$@ | kfusion/thirdparty/nvprof2log.py >   $(LOGS)/kernels.$@
	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$@  livingRoom$(*F).gt.freiburg > $(LOGS)/resume.$@
	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$@   >> $(LOGS)/resume.$@

octree-live.log : 
	./build/kfusion/kfusion-qt-octree-openmp $(live)


tum.octree-openmp.log : fr2_desk.raw 
	KERNEL_TIMINGS=1 ./build/kfusion/kfusion-benchmark-octree-openmp $(fr2desk) -d octree-openmp-$(BLOCK_SIZE).$(SIZE) -i  fr2_desk.raw -o  benchmark.$(BLOCK_SIZE).$(SIZE).$@ 2> kernels.$(BLOCK_SIZE).$(SIZE).$@
	./kfusion/thirdparty/checkPos.py benchmark.$(BLOCK_SIZE).$(SIZE).$@  livingRoom2.gt.freiburg > resume.$(BLOCK_SIZE).$(SIZE).$@
	./kfusion/thirdparty/checkKernels.py kernels.$(BLOCK_SIZE).$(SIZE).$@   >> resume.$(BLOCK_SIZE).$(SIZE).$@
# 
# tum.openmp.log : fr2_desk.raw 
# 	KERNEL_TIMINGS=1 ./build/kfusion/kfusion-benchmark-openmp $(fr2desk) -d octree-openmp-$(BLOCK_SIZE).$(SIZE) -i  fr2_desk.raw -o  benchmark.$(BLOCK_SIZE).$(SIZE).$@ 2> kernels.$(BLOCK_SIZE).$(SIZE).$@
# 	./kfusion/thirdparty/checkPos.py benchmark.$(BLOCK_SIZE).$(SIZE).$@  livingRoom2.gt.freiburg > resume.$(BLOCK_SIZE).$(SIZE).$@
# 	./kfusion/thirdparty/checkKernels.py kernels.$(BLOCK_SIZE).$(SIZE).$@   >> resume.$(BLOCK_SIZE).$(SIZE).$@
# 
# fr2desk.cuda.log  : fr2_desk.raw
# 	$(eval LOGS = logs/noiseless/cuda/$(*F)/$(SIZE)/$(shell date +%F--%T))
# 	mkdir -p $(LOGS)
# 	/usr/local/cuda-7.5/bin/nvprof --print-gpu-trace ./build/kfusion/kfusion-benchmark-cuda $($(*F)) -i  fr2_desk.raw -o  $(LOGS)/benchmark.$@ 2> $(LOGS)/nvprof.$@ || true
# 	cat  $(LOGS)/nvprof.$@ | kfusion/thirdparty/nvprof2log.py >   $(LOGS)/kernels.$@
# 	./kfusion/thirdparty/checkPos.py $(LOGS)/benchmark.$@  livingRoom2.gt.freiburg > $(LOGS)/resume.$@
# 	./kfusion/thirdparty/checkKernels.py $(LOGS)/kernels.$@   >> $(LOGS)/resume.$@

#### GENERAL GENERATION ####

clean :
	rm -rf build
cleanall : 
	rm -rf build
	rm -rf living_room_traj*_loop livingRoom*.gt.freiburg living_room_traj*_loop.raw
	rm -f *.log 


.PHONY : clean bench test all validate

.PRECIOUS: living_room_traj%_loop livingRoom%.gt.freiburg living_room_traj%_loop.raw

