import dadapy
import dadaC
import numpy as np
import matplotlib.pyplot as plt
import time

dsInfo  = []
results = [] 
## ---------------------------------------------------------------------
## Benchmarks and tests dadaC vs dadapy
## ---------------------------------------------------------------------

## ---------------------------------------------------------------------
## Gaussian mixture 2D
## 
## ---------------------------------------------------------------------
def profileAndRun(dataset,dataset_name,k,Z,results):
    dp = dadapy.Data(dataset,verbose=True)
    dc = dadaC.Data(dataset)

    pres = []

    print(f"**** dadapy on {dataset_name}")
    print("dadapy output\n","-"*(30),"\n")

    pres.append(["Method", "part", "time"])
    
    dsInfo.append(f"{dataset_name} {dataset.shape[0]} points")
    t1 = time.monotonic()
    
    dp.compute_distances(k)
    dp.compute_id_2NN()
    dp.compute_density_kstarNN()

    t2 = time.monotonic()
    pres.append(["py","ngbh and density", f"{t2 - t1: .2}s"]) 

    t1 = time.monotonic()

    dp.compute_clustering_ADP(Z = Z, halo = True)

    t2 = time.monotonic()
    pres.append(["py", "ADP", f"{t2 - t1: .2}s"]) 

    print(f"\n**** dadaC on {dataset_name} ****")
    print("dadapy output\n","-"*(30),"\n")

    t1 = time.monotonic()
    dc.computeNeighbors(k)

    dc.computeIDtwoNN()
    dc.computeDensity()


    t2 = time.monotonic()
    pres.append(["C", "ngbh and density", f"{t2 - t1: .2}s"]) 

    t1 = time.monotonic()

    dc.computeClusteringADP(Z = Z, halo = True)

    t2 = time.monotonic()
    pres.append(["C", "ADP", f"{t2 - t1: .2}s"]) 

    results.append(pres)

    print("\n**** Comparing Results ****")
    c1 = dp.cluster_assignment
    c2 = dc.getClusterAssignment()
    errors = np.where(c1 != c2)[0].shape[0]
    print(f" --> \t Found {errors} errors! \n")



disclaimer = '''
                /!\ NOTE: time measures taken with time.monotonic() /!\ 
             '''

print(disclaimer,"\n")
k  = 300
Z = 3

#n  = 50000
#x1 = np.random.normal([0,2],1,size=(n,2)) 
#x2 = np.random.normal([2,0],1,size=(n,2)) 
#x  = np.concatenate([x1,x2])
#
#profileAndRun(x, "2D Gaussian", k, Z, results)
#
#n  = 50000
#x1 = np.random.normal([0,2,0,0,0],1,size=(n,5)) 
#x2 = np.random.normal([2,0,0,0,0],1,size=(n,5)) 
#x  = np.concatenate([x1,x2])
#
#profileAndRun(x, "5D Gaussian", k, Z, results)


from sklearn.datasets import fetch_openml
x, y = fetch_openml(
    "mnist_784", version=1, return_X_y=True, as_frame=False, parser="pandas"
)

x = np.ascontiguousarray(x)
x = x.astype(np.float64)/255.


profileAndRun(x[:10000], "MNIST subset", k, Z, results)

from tabulate import tabulate

print("Printing results \n")

for header, tab in zip(dsInfo,results):
    print(header)
    print(tabulate(tab, headers = "firstrow", tablefmt='fancy_grid'))


#print(results)






