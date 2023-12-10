import dadapy
import numpy as np
import sys
sys.path.append("../")
import dadaC
def test_ngbh_kdtree():
    n = 10000
    k = 300
    d = np.random.choice([ i for i in range(1,10)])
    x = np.random.rand(n,d)
    dpy = dadapy.Data(x)
    dc = dadaC.Data(x)

    dpy.compute_distances(300)
    dc.computeNeighbors(300)

    n1 = dpy.dist_indices
    dist = dpy.distances
    n2 , _ = dc.getNeighbors()
    print(n1.shape)
    print(n2.shape)
    assert np.all(n1 == n2)

if __name__ == "__main__":
    test_ngbh_kdtree()

