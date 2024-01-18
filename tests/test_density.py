import numpy as np
import sys, os
#sys.path.append("./..")
import dadac

x = np.loadtxt(os.path.join(os.path.dirname(__file__), "gt/Fig2.dat"))
gt_id = 2.024278266766415
gt_density_kstar     = np.loadtxt(os.path.join(os.path.dirname(__file__),"gt/gt_density_kstar.dat"))
gt_density_kstar_err = np.loadtxt(os.path.join(os.path.dirname(__file__),"gt/gt_density_kstar_err.dat"))
gt_kstar  = np.loadtxt(os.path.join(os.path.dirname(__file__),"gt/gt_kstar.dat"))

gt_density_PAk      = np.loadtxt(os.path.join(os.path.dirname(__file__),"gt/gt_density_PAk.dat"))
gt_density_PAk_err  = np.loadtxt(os.path.join(os.path.dirname(__file__),"gt/gt_density_PAk_err.dat"))

def test_ID_twoNN():
    dc = dadac.Data(x,verbose = False)
    dc.compute_distances(300)
    dc.compute_id_2NN()
    i2 = dc.id 

    #dpy = dadapy.Data(x)
    #dpy.compute_distances(300)
    #dpy.compute_id_2NN()
    #i1 = dpy.intrinsic_dim

    #print(i2)
    #assert np.abs(i1 - i2) < 1e-9 
    assert np.isclose(i2,gt_id)

def test_density_kstarNN():
    #dpy = dadapy.Data(x)
    #dpy.compute_distances(300)
    #dpy.compute_id_2NN()
    #dpy.compute_density_kstarNN()
    #d1 = dpy.log_den
    #e1 = dpy.log_den_err
    #k1 = dpy.kstar


    dc = dadac.Data(x,verbose = False)
    dc.compute_distances(300, alg="kd")
    dc.compute_id_2NN()
    dc.compute_density_kstarNN()
    d2 = dc.log_den
    e2 = dc.log_den_err
    k2 = dc.kstar
    



    assert np.all(gt_kstar == k2)
    #assert np.all(np.abs(d1 - d2) < 1e-9)
    #assert np.all(np.abs(e1 - e2) < 1e-9)
    assert np.all(np.isclose(gt_density_kstar, d2))
    assert np.all(np.isclose(gt_density_kstar_err, e2))

def test_density_PAk():
    #
    # /!\   PAk implementation differs from the one in dadapy, 
    # /!\   being slitghtly different and beign the tolerance
    # /!\   on Newton method 10^-3 results need to be tested up to 
    # /!\   that tolerance
    # 

    #dpy = dadapy.Data(x)
    #dpy.compute_distances(300)
    #dpy.compute_id_2NN()
    #dpy.compute_density_PAk()
    #d1 = dpy.log_den
    #e1 = dpy.log_den_err
    #k1 = dpy.kstar

    dc = dadac.Data(x,verbose=False)
    dc.compute_distances(300, alg="kd")
    dc.compute_id_2NN()
    dc.compute_density_PAk()
    d2 = dc.log_den
    e2 = dc.log_den_err




    k2 = dc.kstar
    #assert np.all(np.isclose(d1,d2))
    #assert np.all(np.isclose(e1,e2))
    assert np.all(np.isclose(gt_density_PAk,d2))
    assert np.all(np.isclose(gt_density_PAk_err,e2))


if __name__ == "__main__":
    #test_density_kstarNN()
    #test_ID_twoNN()
    test_density_PAk()


