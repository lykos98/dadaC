import numpy as np
import sys, os
#sys.path.append(os.path.join(__file__,"./.."))
#sys.path.append("./..")
from sklearn.neighbors import NearestNeighbors
import dadac
x = np.loadtxt(os.path.join(os.path.dirname(__file__), "gt/Fig2.dat"))
nn = NearestNeighbors(n_neighbors=301, n_jobs=-1)
nn.fit(x)

gt_nn_dist, gt_nn_idx = nn.kneighbors(x) 


def test_ngbh_kdtree():
    dc = dadac.Data(x, verbose=False)
    dc.compute_distances(300, alg="kd")
    n2 = dc.dist_indices
    d2 = dc.distances
    
    assert np.all(n2 == gt_nn_idx)
    assert np.all(np.isclose(d2,gt_nn_dist))

def test_ngbh_vptree():
    dc = dadac.Data(x, verbose=False)
    dc.compute_distances(300, alg="vp")
    n2 = dc.dist_indices
    d2 = dc.distances

    assert np.all(n2 == gt_nn_idx)
    assert np.all(np.isclose(d2,gt_nn_dist))

if __name__ == "__main__":
    test_ngbh_kdtree()
    test_ngbh_vptree()

