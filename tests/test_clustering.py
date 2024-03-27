# import dadapy
import numpy as np
import sys, os
import dadac

x = np.loadtxt(os.path.join(os.path.dirname(__file__), "gt/Fig2.dat"))
gt_cl_ADP = np.loadtxt(os.path.join(os.path.dirname(__file__), "gt/gt_cl_ADP.dat"))


def test_clustering_ADP():
    dc = dadac.Data(x, verbose=False)
    dc.compute_distances(300, alg="kd")
    dc.compute_id_2NN()
    dc.compute_density_kstarNN()
    dc.compute_clustering_ADP(1.65, use_sparse=True)

    c2 = dc.cluster_assignment

    dc = dadac.Data(x, verbose=False)
    dc.compute_distances(300, alg="kd")
    dc.compute_id_2NN()
    dc.compute_density_kstarNN()
    dc.compute_clustering_ADP(1.65, use_sparse=False)

    c3 = dc.cluster_assignment
    # np.savetxt(os.path.join(os.path.dirname(__file__),"gt/gt_cl_ADP.dat"), c2, fmt="%d")

    assert np.all(c2 == gt_cl_ADP)
    assert np.all(c2 == c3)
