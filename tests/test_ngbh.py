import dadapy
import numpy as np
import sys, os
#sys.path.append(os.path.join(__file__,"./.."))
#sys.path.append("./..")
import dadac

def test_ngbh_kdtree():
    n = 10000
    k = 300
    d = np.random.choice([ i for i in range(1,10)])
    x = np.random.rand(n,d)
    dpy = dadapy.Data(x)
    dc = dadac.Data(x)

    dpy.compute_distances(300)
    dc.compute_distances(300, alg="kd")

    n1 = dpy.dist_indices
    dist = dpy.distances
    n2 = dc.dist_indices
    #print(n1.shape)
    #print(n2.shape)
    assert np.all(n1 == n2)

def test_ngbh_vptree():
    n = 10000
    k = 300
    d = np.random.choice([ i for i in range(1,10)])
    x = np.random.rand(n,d)
    dpy = dadapy.Data(x)
    dc = dadac.Data(x)

    dpy.compute_distances(300)
    dc.compute_distances(300, alg="vp")

    n1 = dpy.dist_indices
    dist = dpy.distances
    n2 = dc.dist_indices
    #print(n1.shape)
    #print(n2.shape)
    assert np.all(n1 == n2)

if __name__ == "__main__":
    test_ngbh_kdtree()

