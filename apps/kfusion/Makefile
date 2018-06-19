BLOCK_SIZE = 8
SIZE = 512
C=2
#### ICL-NUIM TRAJECTORIES PARAMETERS  ####
0 = -s 5.0 -p 0.34,0.5,0.24 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F
1 = -s 5.0 -p 0.485,0.5,0.55 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F  
2 = -s 5.0 -p 0.34,0.5,0.24 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F
3 = -s 5.0 -p 0.2685,0.5,0.4 	-z 1 -c $C -r 1 -k 481.2,-480,320,240 -B $(BLOCK_SIZE) -v $(SIZE)  -m 0.1 -F 
live = -s 5.0 -p 0.5,0.5,0.5 -z 1 -c 1 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1
fr2desk = -s 5.0 -p 0.5,0.5,0.2 -z 1 -c 2 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE) -m 0.1 -F
# fr2xyz = 	-s 8.0 -p 0.5,0.5,0.5	-z 1 -c 2 -r 1 -k 580.8,581.8,308.8,253.0 -B $(BLOCK_SIZE) -v $(SIZE),$(SIZE),$(SIZE)  -m 0.1 -F 

ROOT_DIR=$(shell pwd)
TOON_DIR=${ROOT_DIR}/TooN/install_dir
TOON_INCLUDE_DIR=${TOON_DIR}/include/


all : TooN
	mkdir -p build/
	cd build/ && cmake -DCMAKE_BUILD_TYPE=Release \
	 	-DTOON_INCLUDE_PATH=${TOON_INCLUDE_DIR} \
	 	-DSOPHUS_INCLUDE_PATH=${SOPHUS_ROOT} \
		$(CMAKE_ARGUMENTS) ..
	$(MAKE) -C build  $(MFLAGS) $(SPECIFIC_TARGET)


debug: TooN
	mkdir -p build/
	mkdir -p build/logs/
	cd build/ && cmake -DCMAKE_BUILD_TYPE=Debug \
	 	-DTOON_INCLUDE_PATH=${TOON_INCLUDE_DIR} \
	 	-DSOPHUS_INCLUDE_PATH=${SOPHUS_ROOT} \
		$(CMAKE_ARGUMENTS) ..
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

live.log : 
	./build/kfusion-qt-openmp $(live)

demo-ofusion:
	./build/kfusion-main-openmp --compute-size-ratio 2 --fps 0 --block-read False --input-file /data/ev314/data/living_room_traj2_frei_png/scene.raw --icp-threshold 1e-05 --mu 0.008 --init-pose 0.34,0.5,0.24 --integration-rate 1 --volume-size 5 -B 8 --tracking-rate 1 --volume-resolution 512 --pyramid-levels 10,5,4 --rendering-rate 1 -k 481.2,-480,320,240

demo-kfusion:
	./build/kfusion-main-openmp --compute-size-ratio 2 --fps 0 --block-read False --input-file /data/ev314/data/living_room_traj2_frei_png/scene.raw --icp-threshold 1e-05 --mu 0.1 --init-pose 0.34,0.5,0.24 --integration-rate 1 --volume-size 5 -B 8 --tracking-rate 1 --volume-resolution 512 --pyramid-levels 10,5,4 --rendering-rate 1 -k 481.2,-480,320,240


#### GENERAL GENERATION ####

clean :
	rm -rf build
cleanall : 
	rm -rf build
	rm -rf living_room_traj*_loop livingRoom*.gt.freiburg living_room_traj*_loop.raw
	rm -f *.log 


.PHONY : clean bench test all validate

.PRECIOUS: living_room_traj%_loop livingRoom%.gt.freiburg living_room_traj%_loop.raw

