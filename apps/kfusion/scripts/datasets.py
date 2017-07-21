import os

from _dataset import *
from systemsettings import *

#
# TUM RGB-D fr2/xyz Settings
#
TUM_RGB_FR2_XYZ = Dataset()
TUM_RGB_FR2_XYZ.dataset_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_xyz/scene.raw')
TUM_RGB_FR2_XYZ.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_xyz/groundtruth.txt')
#TUM_RGB_FR2_XYZ.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_xyz/scene-ground.txt')
TUM_RGB_FR2_XYZ.camera_file = os.path.join(DATASETS_PATH, 'freiburg2.txt')
TUM_RGB_FR2_XYZ.camera = '580.8,581.8,308.8,253.0'
TUM_RGB_FR2_XYZ.quat = '-0.5721,0.6521,-0.3565,0.3469'
TUM_RGB_FR2_XYZ.init_pose = '0.5,0.5,0.5'
TUM_RGB_FR2_XYZ.rgb_image = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_xyz/rgb/')
TUM_RGB_FR2_XYZ.pre_assoc_file_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_xyz/timing.assoc.txt')
TUM_RGB_FR2_XYZ.ate_associate_identity = False
TUM_RGB_FR2_XYZ.descr = 'fr2_xyz'

#
# TUM RGB-D fr1/xyz Settings
#

TUM_RGB_FR1_XYZ = Dataset()
TUM_RGB_FR1_XYZ.dataset_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_xyz/scene.raw')
TUM_RGB_FR1_XYZ.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_xyz/groundtruth.txt')
#TUM_RGB_FR2_XYZ.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_xyz/scene-ground.txt')
TUM_RGB_FR1_XYZ.camera_file = os.path.join(DATASETS_PATH, 'freiburg1.txt')
TUM_RGB_FR1_XYZ.camera = '591.1,590.1,331.0,234.0'
TUM_RGB_FR1_XYZ.quat = '0.6132,0.5962,-0.3311,-0.3986'
TUM_RGB_FR1_XYZ.init_pose = '0.5,0.5,0.5'
TUM_RGB_FR1_XYZ.rgb_image = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_xyz/rgb/')
TUM_RGB_FR1_XYZ.pre_assoc_file_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_xyz/timing.assoc.txt')
TUM_RGB_FR1_XYZ.ate_associate_identity = False
TUM_RGB_FR1_XYZ.descr = 'fr2_xyz'

#
# TUM RGB-D fr2/desk Settings
#
TUM_RGB_FR2_DESK = Dataset()
TUM_RGB_FR2_DESK.dataset_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_desk/scene.raw')
TUM_RGB_FR2_DESK.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_desk/groundtruth.txt')
TUM_RGB_FR2_DESK.camera_file = os.path.join(DATASETS_PATH, 'freiburg2.txt')
TUM_RGB_FR2_DESK.camera = '580.8,581.8,308.8,253.0'
TUM_RGB_FR2_DESK.quat = '0.6529,-0.5483,0.3248,-0.4095'
TUM_RGB_FR2_DESK.init_pose = '0.5,0.5,0'
TUM_RGB_FR2_DESK.pre_assoc_file_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg2_desk/timings.assoc.txt')
TUM_RGB_FR2_DESK.descr = 'fr2_desk'
TUM_RGB_FR2_DESK.ate_associate_identity = False

#
# TUM RGB-D fr1/desk Settings
#
TUM_RGB_FR1_DESK = Dataset()
TUM_RGB_FR1_DESK.dataset_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_desk/scene.raw')
TUM_RGB_FR1_DESK.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_desk/groundtruth.txt')
TUM_RGB_FR1_DESK.camera = '591.1,590.1,331.0,234.0'
TUM_RGB_FR1_DESK.quat = '0.6529,-0.5483,0.3248,-0.4095'
TUM_RGB_FR1_DESK.init_pose = '0.5,0.5,0'
TUM_RGB_FR1_DESK.pre_assoc_file_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg1_desk/timings.assoc.txt')
TUM_RGB_FR1_DESK.descr = 'fr1_desk'
TUM_RGB_FR1_DESK.ate_associate_identity = False

TUM_RGB_FR3_DESK = Dataset()
TUM_RGB_FR3_DESK.dataset_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg3_long_office_household/scene.raw')
TUM_RGB_FR3_DESK.ground_truth = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg3_long_office_household/groundtruth.txt')
TUM_RGB_FR3_DESK.camera = '567.6,570.2,324.7,250.1'
TUM_RGB_FR3_DESK.init_pose = '0.5,0.5,0'
TUM_RGB_FR3_DESK.pre_assoc_file_path = os.path.join(DATASETS_PATH, 'rgbd_dataset_freiburg3_long_office_household/timings.assoc.txt')
TUM_RGB_FR3_DESK.descr = 'fr3_desk'
TUM_RGB_FR3_DESK.ate_associate_identity = False


#
# ICL-NUIM Living Room 0
#
ICL_NUIM_LIV_0 = Dataset()
ICL_NUIM_LIV_0.dataset_path = os.path.join(DATASETS_PATH, 'living_room_traj0_loop/scene.raw')
ICL_NUIM_LIV_0.ground_truth = os.path.join(DATASETS_PATH, 'living_room_traj0_loop/livingRoom0.gt.freiburg.txt')
ICL_NUIM_LIV_0.camera_file = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/camera.txt')
ICL_NUIM_LIV_0.camera = '481.2,-480,320,240'
ICL_NUIM_LIV_0.init_pose = '0.34,0.5,0.24'
ICL_NUIM_LIV_0.ate_associate_identity = True
ICL_NUIM_LIV_0.descr = 'liv_traj_0'


#
# ICL-NUIM Living Room 1
#
ICL_NUIM_LIV_1 = Dataset()
ICL_NUIM_LIV_1.dataset_path = os.path.join(DATASETS_PATH, 'living_room_traj1_loop/scene.raw')
ICL_NUIM_LIV_1.ground_truth = os.path.join(DATASETS_PATH, 'living_room_traj1_loop/livingRoom1.gt.freiburg')
ICL_NUIM_LIV_1.camera_file = os.path.join(DATASETS_PATH, 'living_room_traj1_loop/camera.txt')
ICL_NUIM_LIV_1.camera = '481.2,-480,320,240'
ICL_NUIM_LIV_1.init_pose = '0.485,0.5,0.55'
ICL_NUIM_LIV_1.ate_associate_identity = True
ICL_NUIM_LIV_1.descr = 'liv_traj_1'


#
# ICL-NUIM Living Room 2
#
ICL_NUIM_LIV_2 = Dataset()
ICL_NUIM_LIV_2.dataset_path = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/scene.raw')
ICL_NUIM_LIV_2.ground_truth = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/livingRoom2.gt.freiburg.txt')
ICL_NUIM_LIV_2.camera_file = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/camera.txt')
ICL_NUIM_LIV_2.camera = '481.2,-480,320,240'
ICL_NUIM_LIV_2.init_pose = '0.34,0.5,0.24'
ICL_NUIM_LIV_2.ate_associate_identity = True
ICL_NUIM_LIV_2.descr = 'liv_traj_2'

#
# ICL-NUIM Living Room 3
#
ICL_NUIM_LIV_3 = Dataset()
ICL_NUIM_LIV_3.dataset_path = os.path.join(DATASETS_PATH, 'living_room_traj3_loop/scene.raw')
ICL_NUIM_LIV_3.ground_truth = os.path.join(DATASETS_PATH, 'living_room_traj3_loop/livingRoom3.gt.freiburg')
ICL_NUIM_LIV_3.camera_file = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/camera.txt')
ICL_NUIM_LIV_3.camera = '481.2,-480,320,240'
ICL_NUIM_LIV_3.init_pose = '0.2685,0.5,0.4'
ICL_NUIM_LIV_3.ate_associate_identity = True
ICL_NUIM_LIV_3.descr = 'liv_traj_3'

#
# ICL-NUIM Office 0
#
ICL_NUIM_OFF_0 = Dataset()
ICL_NUIM_OFF_0.dataset_path = os.path.join(DATASETS_PATH, 'office_room_traj0_loop/scene.raw')
ICL_NUIM_OFF_0.ground_truth = os.path.join(DATASETS_PATH, 'office_room_traj0_loop/traj0.gt.freiburg')
ICL_NUIM_OFF_0.camera_file = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/camera.txt')
ICL_NUIM_OFF_0.camera = '481.2,480,320,240'
ICL_NUIM_OFF_0.init_pose = '0.5,0.5,0.5'


#
# ICL-NUIM Office 2
#
ICL_NUIM_OFF_2 = Dataset()
ICL_NUIM_OFF_2.dataset_path = os.path.join(DATASETS_PATH, 'office_room_traj2_loop/scene.raw')
ICL_NUIM_OFF_2.ground_truth = os.path.join(DATASETS_PATH, 'office_room_traj2_loop/traj2.gt.freiburg.txt')
ICL_NUIM_OFF_2.camera_file = os.path.join(DATASETS_PATH, 'living_room_traj2_loop/camera.txt')
ICL_NUIM_OFF_2.camera = '481.2,480,320,240'
ICL_NUIM_OFF_2.init_pose = '0.5,0.5,0.5'

