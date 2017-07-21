#!/usr/bin/python2
from _run import * 
from systemsettings import *
from datasets import *
import numpy as np 


if __name__ == "__main__":
    results_dir = gen_results_dir(RESULTS_PATH)
    algorithm = KinectFusion(BIN_PATH)
    # -q --fps 10 --block-read 1

    # Find the best alignment between the gt and the computed trajectory.
    # It influences results a lot, we should really discuss this.
    algorithm.ate_align = True

    # All must be true for ICL-NUIM
    algorithm.ate_associate_identity = False  # 1to1 association gt-tra
    # When true the trajectory is offset by the first position.
    algorithm.ate_remove_offset = False
    algorithm.voxel_block = '8'
    algorithm.rendering_rate = '1'
    algorithm.bilateralFilter = True
    algorithm.init_pose = '0.5,0.5,0.2'
    min_ate = 100.0
    run_results = {}
    # for mu in [0.1, 0.05]:
    for mu in [0.1]:
        for sequence in [TUM_RGB_FR3_DESK, TUM_RGB_FR2_DESK]:
            for resol in [512]:
                for version in ['openmp']:
                    kernel_data = []
                    algorithm.impl = version
                    algorithm.volume_resolution = str(resol)
                    algorithm.volume_size = '5'
                    algorithm.compute_size_ratio = 2
                    algorithm.integration_rate = 1
                    algorithm.mu = mu
                    res = algorithm.run(sequence)
                    run_results[sequence.descr] = res
                    kernel_data.append(res)

    with open(results_dir + '/resume.log', 'w') as f:
        f.write('{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\n'.format('dataset', 'ATE',
                'preprocessing',
                'integration',
                'raycasting',
                'computation'))
        for k, d in run_results.iteritems():
            data = d['data']
            f.write('{:>10}\t{:>10.4f}\t{:>10.4f}\t{:>10.4f}\t{:>10.4f}\t{:>10.4f}\n'.format(k, float(d['ate_mean']),
                    float(data['preprocessing']['mean']),
                    float(data['integration']['mean']),
                    float(data['raycasting']['mean']),
                    float(data['computation']['mean'])))
