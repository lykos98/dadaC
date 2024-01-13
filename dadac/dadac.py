import ctypes as ct
import numpy as np  
import os
from sklearn.neighbors import NearestNeighbors
import time
from wurlitzer import sys_pipes


ctFloatType = ct.c_double
ctIdxType = ct.c_uint64

class HeapNode(ct.Structure):
    _fields_ = [
        ("value", ctFloatType),
        ("array_idx", ctIdxType)
    ]


class Heap(ct.Structure):
    _fields_ = [
        ("N", ctIdxType),
        ("count", ctIdxType),
        ("data", ct.POINTER(HeapNode))
    ]

class luDynamicArray(ct.Structure):
    _fields_ = [
        ("data", ct.POINTER(ctIdxType)),
        ("size", ctIdxType),
        ("count", ctIdxType)
    ]

class DatapointInfo(ct.Structure):
    _fields_ = [
        ("g", ctFloatType),
        ("ngbh", Heap),
        ("array_idx", ctIdxType),
        ("log_den", ctFloatType),
        ("log_den_c", ctFloatType),
        ("log_den_err", ctFloatType),
        ("kstar", ctIdxType),
        ("is_center", ct.c_int),
        ("cluster_idx", ct.c_int)
    ]

    def _repr_(self):
        return f"{{ g: {self.g}, ngbh: {self.ngbh}, array_idx: {self.array_idx}, log_den: {self.log_den}, log_den_c: {self.log_den_c}, log_den_err: {self.log_den_err}, kstar: {self.kstar}, isCenter: {self.is_center}, clusterIdx: {self.cluster_idx} }}"

class Border_t(ct.Structure):
    _fields_ = [
        ("idx", ctIdxType),
        ("density", ctFloatType),
        ("error", ctFloatType)
    ]

class SparseBorder_t(ct.Structure):
    _fields_ = [
        ("i", ctIdxType),
        ("j", ctIdxType),
        ("idx", ctIdxType),
        ("density", ctFloatType),
        ("error", ctFloatType)
    ]

class AdjList(ct.Structure):
    _fields_ = [
        ("count", ctIdxType),
        ("size", ctIdxType),
        ("data", ct.POINTER(SparseBorder_t))
    ]

class Clusters(ct.Structure):
    _fields_ = [
        ("UseSparseBorders", ct.c_int),
        ("SparseBorders", ct.POINTER(AdjList)),
        ("centers", luDynamicArray),
        ("borders", ct.POINTER(ct.POINTER(Border_t))),
        ("_borders_data", ct.POINTER(Border_t)),
        ("n", ctIdxType)
    ]
class _dadac_loader():
    def __init__(self):
        path = os.path.join(os.path.dirname(__file__), "bin/libdadac.so")
        print(path)
        if not os.path.exists(path):
            print("dadac is not built yet, calling make for you")
            try:
                os.system(f"make -C {os.path.dirname(__file__)}")
            except:
                print("Cannot build dadac")

        self.lib = ct.CDLL(path)

        global ctFloatType, ctIdxType
        s = self.lib.FloatAndUintSize()

        self._useFloat32 = s < 1
        self._useInt32   = (s == 0) or (s == 2)


        if self._useInt32:
            ctIdxType = ct.c_uint32

        #retrieve function pointers form .so file

        self._NgbhSearch_kdtree = self.lib.NgbhSearch_kdtree_V2
        self._NgbhSearch_kdtree.argtypes = [np.ctypeslib.ndpointer(ctFloatType), ct.c_uint64, ct.c_uint64, ct.c_uint64 ]
        self._NgbhSearch_kdtree.restype  = ct.POINTER(DatapointInfo)

        #Datapoint_info* NgbhSearch_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        self._NgbhSearch_vptree = self.lib.NgbhSearch_vptree_V2
        #actually do not define anything and hope for the best
        #function pointers are not documented
        #METRIC = ct.CFUNCTYPE(x) 
        self._NgbhSearch_vptree.argtypes = [ct.c_void_p, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_void_p ]
        self._NgbhSearch_vptree.restype  = ct.POINTER(DatapointInfo)
        self._eud = self.lib.eud

        self._NgbhSearch_bruteforce = self.lib.NgbhSearch_bruteforce
        #actually do not define anything and hope for the best
        #function pointers are not documented
        #METRIC = ct.CFUNCTYPE(x) 
        self._NgbhSearch_bruteforce.argtypes = [ct.c_void_p, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_void_p ]
        self._NgbhSearch_bruteforce.restype  = ct.POINTER(DatapointInfo)
        self._eud_sq = self.lib.eud_sq

        #void setRhoErrK(Datapoint_info* points, FLOAT_TYPE* rho, FLOAT_TYPE* rhoErr, idx_t* k, size_t n)
        self._setRhoErrK = self.lib.setRhoErrK
        self._setRhoErrK.argtypes = [  ct.POINTER(DatapointInfo),    
                                        np.ctypeslib.ndpointer(ctFloatType), 
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        np.ctypeslib.ndpointer(ctIdxType),
                                        ct.c_uint64 ]


        #Datapoint_info* _allocDatapoints(idx_t* indeces, float_t* distances, idx_t n, idx_t k)

        self._allocateDatapoints = self.lib.allocDatapoints
        self._allocateDatapoints.argtypes = [ctIdxType]
        self._allocateDatapoints.restype  = ct.POINTER(DatapointInfo)


        self._importNeighborsAndDistances = self.lib.importNeighborsAndDistances

        #void _importNeighborsAndDistances(Datapoint_info* points, idx_t* indeces, float_t* distances, idx_t n, idx_t k)
        self._importNeighborsAndDistances.argtypes =   [   ct.POINTER(DatapointInfo),
                                                            np.ctypeslib.ndpointer(ctIdxType),
                                                            np.ctypeslib.ndpointer(ctFloatType),
                                                            ctIdxType,
                                                            ctIdxType
                                                        ]
        #void exportNeighborsAndDistances(Datapoint_info* points, idx_t* dist_indices, float_t* dists, idx_t n, idx_t k)
        self._exportNeighborsAndDistances = self.lib.exportNeighborsAndDistances
        self._exportNeighborsAndDistances.argtypes =   [   ct.POINTER(DatapointInfo),
                                                            np.ctypeslib.ndpointer(ctIdxType),
                                                            np.ctypeslib.ndpointer(ctFloatType),
                                                            ctIdxType,
                                                            ctIdxType
                                                        ]

        #void _importDensity(Datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
        self._importDensity = self.lib.importDensity
        self._importDensity.argtypes = [   ct.POINTER(DatapointInfo),
                                            np.ctypeslib.ndpointer(ctIdxType),
                                            np.ctypeslib.ndpointer(ctFloatType),
                                            np.ctypeslib.ndpointer(ctFloatType),
                                            ctIdxType
                                        ]

        self._exportDensity = self.lib.exportDensity
        self._exportDensity.argtypes = [   ct.POINTER(DatapointInfo),
                                            np.ctypeslib.ndpointer(ctIdxType),
                                            np.ctypeslib.ndpointer(ctFloatType),
                                            np.ctypeslib.ndpointer(ctFloatType),
                                            ctIdxType
                                        ]

        #void computeAvg(Datapoint_info* p, FLOAT_TYPE *va, FLOAT_TYPE* ve, FLOAT_TYPE* vals, FLOAT_TYPE* verr, size_t k, size_t n)
        self._computeAvg = self.lib.computeAvg
        self._computeAvg.argtypes = [  ct.POINTER(DatapointInfo),    
                                        np.ctypeslib.ndpointer(ctFloatType), 
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        ctIdxType,
                                        ct.c_uint64 ]
        #void exportClusterAssignment(Datapoint_info* points, int* labels, idx_t n)

        self._exportClusterAssignment = self.lib.exportClusterAssignment
        self._exportClusterAssignment.argtypes = [ ct.POINTER(DatapointInfo),    
                                                    np.ctypeslib.ndpointer(np.int32), 
                                                    ct.c_uint64 ]


        #void exportBorders(Clusters* clusters, idx_t* border_idx, float_t* border_den, float_t* border_err)
        self._exportBorders = self.lib.exportBorders
        self._exportBorders.argtypes = [ct.POINTER(Clusters),
                                        np.ctypeslib.ndpointer(np.int32),
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        np.ctypeslib.ndpointer(ctFloatType),
                                       ]


        self._idEstimate = self.lib.idEstimate
        self._idEstimate.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64, ctFloatType]
        self._idEstimate.restype  =  ct.c_double


        self._computeRho = self.lib.computeRho
        self._computeRho.argtypes = [ct.POINTER(DatapointInfo), ct.c_double, ct.c_uint64]

        self._PAk = self.lib.PAk
        self._PAk.argtypes = [ct.POINTER(DatapointInfo), ct.c_double, ct.c_uint64]

        self._computeCorrection = self.lib.computeCorrection
        self._computeCorrection.argtypes = [ct.POINTER(DatapointInfo), ctIdxType, ct.c_double]

        self._H1 = self.lib.Heuristic1
        self._H1.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64]
        self._H1.restype = Clusters

        self._ClustersAllocate = self.lib.Clusters_allocate
        self._ClustersAllocate.argtypes = [ct.POINTER(Clusters), ct.c_int]

        self._blas_in_use = self.lib.blas_are_in_use
        self._blas_in_use.restypes = ct.c_int32

        self._H2 = self.lib.Heuristic2
        self._H2.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo)]

        self._H3 = self.lib.Heuristic3
        self._H3.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo), ct.c_double, ct.c_int]
        
        self._freeDatapoints = self.lib.freeDatapointArray
        self._freeDatapoints.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64]

        self._freeClusters = self.lib.Clusters_free
        self._freeClusters.argtypes = [ct.POINTER(Clusters)]


class Data(_dadac_loader):
    def __init__(self, data : np.array):
        """Data object for dadaC library

        Args:
            data (np.array): 2d array to use in searching for clusters 

        Raises:
            TypeError: Raises TypeError if a type different from a matrix is passed 
        """
        super().__init__()
        #initialize class
        if self._useFloat32:
            self._ftype = np.float32
            self.data = data.astype(np.float32)
            self.data = np.ascontiguousarray(self.data, dtype = np.float32)
            ctFloatType = ct.c_float
        else:
            #self.data = np.ascontiguousarray(data, dtype =np.float64)
            self._ftype = np.float64
            self.data = data.astype(np.float64)
            self.data = np.ascontiguousarray(self.data, dtype = np.float64)

        if self._useInt32:
            self._itype = np.uint32
        else:
            self._itype = np.uint64
        self.state = {
                "ngbh" : False,
                "id" : False,
                "density" : False,
                "clustering" : False,
                "useFloat32": self._useFloat32,
                "useInt32": self._useInt32,
                "useSparse": None,
                "computeHalo": None 
                }
        

        if len(self.data.shape) != 2:
            raise TypeError("Please provide a 2d numpy array")


        self._datapoints     = None
        self._clusters       = None
        self.n              = self.data.shape[0]
        self.dims           = self.data.shape[1]
        self.k              = None
        self.id             = None

        self._clusterAssignment = None
        self._distances         = None
        self._dist_indices      = None
        self._log_den_bord      = None
        self._log_den_bord_err  = None
        self._border_indices    = None
        self._N_clusters        = None
        self._cluster_centers   = None

        self.borders            = None
        self._log_den           = None
        self._log_den_err       = None
        self.blas               = self._blas_in_use() != 0 

    def _is_notebook(self) -> bool:
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                print("You are running in a notebook maybe the timing output will break, but everything should be fine ")
                return True   # Jupyter notebook or qtconsole
            elif shell == "google.colab._shell":
                print("You are running in a google colab notebook maybe the timing output will break, but everything should be fine ")
                return True
            elif shell == 'TerminalInteractiveShell':
                return False  # Terminal running IPython
            else:
                return False  # Other type (?)
        except NameError:
            return False

    def computeNeighbors_kdtree(self, k : int):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
            alg (str): default "kd" for kdtree else choose "vp" for vptree
        """
        self.k = k
        #with sys_pipes():
        if self._is_notebook():
            with sys_pipes():
                self._datapoints = self._NgbhSearch_kdtree(self.data, self.n, self.dims, self.k)
        else:
            self._datapoints = self._NgbhSearch_kdtree(self.data, self.n, self.dims, self.k)
        #Datapoint_info* NgbhSearch_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        self.state["ngbh"] = True
        self.neighbors = None

    def computeAvgOverNgbh(self, vals, error,k):
        retVals = np.zeros_like(vals)
        retError = np.zeros_like(error)
        self._computeAvg(self._datapoints, retVals, retError, vals, error, k, self.n)
        return retVals, retError

    def computeNeighbors_vptree(self, k : int, alg="kd"):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
            alg (str): default "kd" for kdtree else choose "vp" for vptree
        """
        self.k = k
        #Datapoint_info* NgbhSearch_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        if self._is_notebook():
            with sys_pipes():
                self._datapoints = self._NgbhSearch_vptree(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, self._eud)
        else:
            self._datapoints = self._NgbhSearch_vptree(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, self._eud)

        self.state["ngbh"] = True
        self.neighbors = None

    def computeNeighbors_bruteforce(self, k : int, alg="kd"):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
            alg (str): default "kd" for kdtree else choose "vp" for vptree
        """
        self.k = k
        #Datapoint_info* NgbhSearch_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        if self._is_notebook():
            with sys_pipes():
                self._datapoints = self._NgbhSearch_bruteforce(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, 0)
        else:
            self._datapoints = self._NgbhSearch_bruteforce(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, 0)
        self.state["ngbh"] = True
        self.neighbors = None

    def compute_distances(self,k : int, alg = "auto"):
        k = k + 1
        if alg == "auto": 
            if (self.data.shape[1] > 15 or k > self.data.shape[0]//2):
                alg = "bf"
            else:
                alg = "kd"
        
        if alg == "kd":
            self.computeNeighbors_kdtree(k)
            return
        if alg == "vp":
            self.computeNeighbors_vptree(k)
            return
        if alg == "bf":
            self.k = k
            if self.blas:
                self.computeNeighbors_bruteforce(k)
            else:
                print("dadac implementation not compiled with blas support for euclidean metric optimization")
                print("if you want to use it consider compiling against a blas library implementation")
                print("--> Falling back to sklearn brute force")
                t1 = time.monotonic() 
                nn = NearestNeighbors(n_neighbors=k, n_jobs=-1, p = 2, algorithm="brute").fit(self.data)
                dist, ngbh = nn.kneighbors(self.data)
                ngbh = ngbh.astype(self._itype)
                print(ngbh.dtype)
                dist = np.ascontiguousarray(dist.astype(self._ftype),dtype = self._ftype)
                ngbh = np.ascontiguousarray(ngbh)
                
                self._datapoints = self._allocateDatapoints(self.n)
                self._importNeighborsAndDistances(self._datapoints,ngbh,dist,self.n, self._itype(k))
                t2 = time.monotonic()

            self.state["ngbh"] = True
            #print(f"\tTotal time: {t2 - t1 : .2f}s")
            #self.computeNeighbors_bruteforce(k)
            return

    def compute_id_2NN(self,fraction = 0.9):

        """ Compute the intrinsic dimension of the dataset via the TWO Nearest Neighbors method.
            Ref. paper 

        Raises:
            ValueError: Raises value error if neighbors are not computed yet, use `Data.computeNeighbors()` 
        """

        if not self.state["ngbh"]:
            raise ValueError("Please compute Neighbors before calling this function")
        if self._is_notebook():
            with sys_pipes():
                self.id = self._idEstimate(self._datapoints,self.n,fraction)
        else:
            self.id = self._idEstimate(self._datapoints,self.n,fraction)
        self.state["id"] = True

    def import_neighbors_and_distances(self,ngbh,dists):
        if self._datapoints is None:
            self._datapoints = self._allocateDatapoints(self.n)
        k = ngbh.shape[1]
        self._importNeighborsAndDistances(self._datapoints,ngbh.astype(self._itype),dists.astype(self._ftype),self.n, np.uint64(k))

    def import_density(self, log_den,log_den_err,kstar):
        if self._datapoints is None:
            self._datapoints = self._allocateDatapoints(self.n)
        self.state["density"] = True
        #void _importDensity(Datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
        self._log_den = log_den.astype(self._ftype)
        self._log_den_err = log_den_err.astype(self._ftype)
        self._kstar = kstar.astype(self._itype)
        self._importDensity(self._datapoints,self.kstar,self.log_den, self.log_den_err,self.n)

    def compute_density_kstarNN(self):

        """Compute density value for each point

        Raises:
            ValueError: Raises value error if ID is not computed, use `Data.computeIDtwoNN` method
        """
        if not self.state["id"]:
            raise ValueError("Please compute ID before calling this function")
        if self._is_notebook():
            with sys_pipes():
                self._computeRho(self._datapoints, self.id, self.n)
        else:
            self._computeRho(self._datapoints, self.id, self.n)

        self.state["density"] = True
        self.density = None
        self.densityError = None

    def compute_density_PAk(self):

        """Compute density value for each point

        Raises:
            ValueError: Raises value error if ID is not computed, use `Data.computeIDtwoNN` method
        """
        if not self.state["id"]:
            raise ValueError("Please compute ID before calling this function")
        if self._is_notebook():
            with sys_pipes():
                self._PAk(self._datapoints, self.id, self.n)
        else:
            self._PAk(self._datapoints, self.id, self.n)

        self.state["density"] = True
        self.density = None
        self.densityError = None

    def compute_clustering_ADP(self,Z : float, halo = False, useSparse = "auto"):

        """Compute clustering via the Advanced Density Peak method

        Args:
            Z (float): Z value for the method 
            useSparse (str): optional [``auto``,`yes`,`no`], use sparse implementation of border storage between clusters. Memory usage for big datsets is significant. 

        Raises:
            ValueError: Raises value error if density is not computed, use `Data.computeDensity()` method
        """
        if not self.state["density"]:
            raise ValueError("Please compute density before calling this function")
        if useSparse == "auto":
            if self.n > 2e6:
                self.state["useSparse"] = True
            else: 
                self.state["useSparse"] = False 
        elif useSparse == "y":
            self.state["useSparse"] = True
        else:
            self.state["useSparse"] = False

        self.state["computeHalo"] = halo 
        self.Z = Z
        if self._is_notebook():
            with sys_pipes():
                self._computeCorrection(self._datapoints, self.n, self.Z)
                self._clusters = self._H1(self._datapoints, self.n)
                self._ClustersAllocate(ct.pointer(self._clusters), 1 if self.state["useSparse"] else 0)
                self._H2(ct.pointer(self._clusters), self._datapoints)
                self._H3(ct.pointer(self._clusters), self._datapoints, self.Z, 1 if halo else 0 )
        else:
            self._computeCorrection(self._datapoints, self.n, self.Z)
            self._clusters = self._H1(self._datapoints, self.n)
            self._ClustersAllocate(ct.pointer(self._clusters), 1 if self.state["useSparse"] else 0)
            self._H2(ct.pointer(self._clusters), self._datapoints)
            self._H3(ct.pointer(self._clusters), self._datapoints, self.Z, 1 if halo else 0 )
        self.state["clustering"] = True
        self.clusterAssignment = None

    def _getClusterAssignment(self):


        """Retrieve cluster assignment

        Raises:
            ValueError: Raises error if clustering is not computed, use `Data.computeClusteringADP(Z)` 

        Returns:
            List of cluster labels
            
        """
        if self.state["clustering"]:
            self._clusterAssignment = np.zeros(self.n, np.int32)
            self._exportClusterAssignment(self._datapoints,self._clusterAssignment,self.n)
        else:
            raise ValueError("Clustering is not computed yet")

    @property
    def cluster_assignment(self):
        if self._clusterAssignment is None:
            self._getClusterAssignment()
        return self._clusterAssignment


    def _getDensity(self):

        """Retrieve list of density values

        Raises:
            ValueError: Raise error if density is not computed,  

        Returns:
            List of density values
            
        """
        if self.state["density"]:
            self._log_den      = np.zeros(self.n,self._ftype)
            self._log_den_err  = np.zeros(self.n,self._ftype)
            self._kstar        = np.zeros(self.n,self._itype)
            self._exportDensity(self._datapoints,self._kstar, self._log_den, self._log_den_err, self.n)
        else:
            raise ValueError("Density is not computed yet use `Data.computeDensity()`")

    def _getNclusters(self):
        if self.state["clustering"]:
            self._N_clusters = len(set(self.cluster_assignment))
        else:
            raise ValueError("Borders not computed yet use `Data.compute_clustering_[DP,ADP,...]()`")


    def _getBorders(self):

        """Retrieve borders 

        Raises:
            ValueError: Raise error if density is not computed,  

        Returns:
            List of density values
            
        """
        if self.state["clustering"]:
            self._log_den_bord      = np.zeros((self.N_clusters,self.N_clusters),self._ftype) 
            self._log_den_bord_err  = np.zeros((self.N_clusters,self.N_clusters),self._ftype)
            self._border_indices    = np.zeros((self.N_clusters,self.N_clusters),np.int32) - 1
            self._exportBorders(self._clusters,self._border_indices, self._log_den_bord, self._log_den_bord_err)
            #ugly thing to retrieve correct thing in dadapy as border density
            mm = np.zeros_like(self._log_den_bord)
            f = np.where(self._border_indices == -1)
            mm[f] = -1
            mm += np.eye(self._N_clusters)
            self._log_den_bord -= mm
        else:
            raise ValueError("Borders not computed yet use `Data.compute_clustering_[DP,ADP,...]()`")

    def _getClusterCenters(self):
        if self.state["clustering"]:
            self._cluster_centers = [int(self._datapoints[i].array_idx)  for i in range(self.n) if bool(self._datapoints[i].is_center)]
        else:
            raise ValueError("Clustering not computed yet use `Data.compute_clustering_[DP,ADP,...]()`")

    @property
    def cluster_centers(self):
        if self._cluster_centers is None:
            self._getClusterCenters()
        return self._cluster_centers

    @property
    def N_clusters(self):
        if self._N_clusters is None:
            self._getNclusters()
        return self._N_clusters


    @property
    def log_den(self):
        if self._log_den is None:
            self._getDensity()
        return self._log_den

    @property
    def log_den_err(self):
        if self._log_den_err is None:
            self._getDensity()
        return self._log_den_err

    @property
    def kstar(self):
        if self._kstar is None:
            self._getDensity()
        return self._kstar

    
    @property
    def distances(self): 
        if self._distances is None:
            self._distances = np.zeros((self.n,self.k), self._ftype)
            self._dist_indices = np.zeros((self.n,self.k), self._itype)
            self._exportNeighborsAndDistances(self._datapoints,self._dist_indices, self._distances, self.n, self.k)
        return self._distances

    @property
    def dist_indices(self): 
        if self._dist_indices is None:
            self._distances = np.zeros((self.n,self.k), self._ftype)
            self._dist_indices = np.zeros((self.n,self.k), self._itype)
            self._exportNeighborsAndDistances(self._datapoints,self._dist_indices,self._distances, self.n, self.k)
        return self._dist_indices

    @property
    def log_den_bord(self):
        if self._log_den_bord is None:
            self._getBorders()
        return self._log_den_bord
    
    @property
    def log_den_bord_err(self):
        if self._log_den_bord_err is None:
            self._getBorders()
        return self._log_den_bord_err

    @property
    def border_indices(self):
        if self._border_indices is None:
            self._getBorders()
        return self._border_indices
    
    #def setRhoErrK(self, rho, err, k):
    #    self._setRhoErrK(self._datapoints, rho, err, k, self.n)
    #    self.state["density"] = True
    #    return

    def __del__(self):
        if not self._datapoints is None:
            self._freeDatapoints(self._datapoints, self.n)
        if not self._clusters is None:
            self._freeClusters(ct.pointer(self._clusters))
