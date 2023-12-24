import dadapy
import numpy as np
import sys, os
#sys.path.append("./..")
import dadac

def test_ID_twoNN():
    n = 10000
    k = 300
    d = np.random.choice([ i for i in range(1,10)])
    x = np.random.rand(n,d)
    dpy = dadapy.Data(x)
    dc = dadac.Data(x)

    dpy.compute_distances(300)
    dc.computeNeighbors(300, alg="kd")

    dpy.compute_id_2NN()
    dc.computeIDtwoNN()

    i1 = dpy.intrinsic_dim
    i2 = dc.id 
    assert np.abs(i1 - i2) < 1e-9 

def test_density_kstarNN():
    n = 10000
    k = 300
    d = np.random.choice([ i for i in range(1,10)])
    x = np.random.rand(n,d)
    dpy = dadapy.Data(x)
    dc = dadac.Data(x)

    dpy.compute_distances(300)
    dc.computeNeighbors(300, alg="kd")

    dpy.compute_id_2NN()
    dc.computeIDtwoNN()

    dpy.compute_density_kstarNN()
    dc.computeDensity()

    d1 = dpy.log_den
    d2 = dc.getDensity()

    e1 = dpy.log_den_err
    e2 = dc.getDensityError()

    k1 = dpy.kstar
    k2 = dc.getKstar()
    assert np.all(k1 == k2)
    assert np.all(np.abs(d1 - d2) < 1e-9)
    assert np.all(np.abs(e1 - e2) < 1e-9)




    assert True
