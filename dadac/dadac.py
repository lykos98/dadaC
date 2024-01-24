import ctypes as ct
import numpy as np  
import os
from sklearn.neighbors import NearestNeighbors
import time
from wurlitzer import sys_pipes


ct_float_t = ct.c_double
ct_idx_t = ct.c_uint64

class heap_node(ct.Structure):
    _fields_ = [
        ("value", ct_float_t),
        ("array_idx", ct_idx_t)
    ]


class heap(ct.Structure):
    _fields_ = [
        ("N", ct_idx_t),
        ("count", ct_idx_t),
        ("data", ct.POINTER(heap_node))
    ]

class lu_dynamic_array(ct.Structure):
    _fields_ = [
        ("data", ct.POINTER(ct_idx_t)),
        ("size", ct_idx_t),
        ("count", ct_idx_t)
    ]

class datapoint_info(ct.Structure):
    _fields_ = [
        ("g", ct_float_t),
        ("ngbh", heap),
        ("array_idx", ct_idx_t),
        ("log_den", ct_float_t),
        ("log_den_c", ct_float_t),
        ("log_den_err", ct_float_t),
        ("kstar", ct_idx_t),
        ("is_center", ct.c_int),
        ("cluster_idx", ct.c_int)
    ]

    def _repr_(self):
        return f"{{ g: {self.g}, ngbh: {self.ngbh}, array_idx: {self.array_idx}, log_den: {self.log_den}, log_den_c: {self.log_den_c}, log_den_err: {self.log_den_err}, kstar: {self.kstar}, isCenter: {self.is_center}, clusterIdx: {self.cluster_idx} }}"

class border_t(ct.Structure):
    _fields_ = [
        ("idx", ct_idx_t),
        ("density", ct_float_t),
        ("error", ct_float_t)
    ]

class sparse_border_t(ct.Structure):
    _fields_ = [
        ("i", ct_idx_t),
        ("j", ct_idx_t),
        ("idx", ct_idx_t),
        ("density", ct_float_t),
        ("error", ct_float_t)
    ]

class adj_list(ct.Structure):
    _fields_ = [
        ("count", ct_idx_t),
        ("size", ct_idx_t),
        ("data", ct.POINTER(sparse_border_t))
    ]

class clusters(ct.Structure):
    _fields_ = [
        ("use_sparse_borders", ct.c_int),
        ("sparse_borders", ct.POINTER(adj_list)),
        ("centers", lu_dynamic_array),
        ("borders", ct.POINTER(ct.POINTER(border_t))),
        ("_borders_data", ct.POINTER(border_t)),
        ("n", ct_idx_t)
    ]
class _dadac_loader():
    def __init__(self):
        path = os.path.join(os.path.dirname(__file__), "bin/libdadac.so")
        #print(path)
        if not os.path.exists(path):
            print("dadac is not built yet, calling make for you")
            try:
                os.system(f"make -C {os.path.dirname(__file__)}")
            except:
                print("Cannot build dadac")

        self.lib = ct.CDLL(path)

        global ct_float_t, ct_idx_t
        s = self.lib.float_and_uint_size()

        self._use_float32 = s < 1
        self._use_int32   = (s == 0) or (s == 2)


        if self._use_int32:
            ct_idx_t = ct.c_uint32

        #retrieve function pointers form .so file

        self._ngbh_search_kdtree = self.lib.ngbh_search_kdtree_v2
        self._ngbh_search_kdtree.argtypes = [np.ctypeslib.ndpointer(ct_float_t), ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_int32 ]
        self._ngbh_search_kdtree.restype  = ct.POINTER(datapoint_info)

        #Datapoint_info* ngbh_search_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        self._ngbh_search_vptree = self.lib.ngbh_search_vptree_v2
        #actually do not define anything and hope for the best
        #function pointers are not documented
        #METRIC = ct.CFUNCTYPE(x) 
        self._ngbh_search_vptree.argtypes = [ct.c_void_p, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_void_p, ct.c_int32 ]
        self._ngbh_search_vptree.restype  = ct.POINTER(datapoint_info)
        self._eud = self.lib.eud

        self._ngbh_search_bruteforce = self.lib.ngbh_search_bruteforce
        #actually do not define anything and hope for the best
        #function pointers are not documented
        #METRIC = ct.CFUNCTYPE(x) 
        self._ngbh_search_bruteforce.argtypes = [ct.c_void_p, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_uint64, ct.c_void_p, ct.c_int32 ]
        self._ngbh_search_bruteforce.restype  = ct.POINTER(datapoint_info)
        self._eud_sq = self.lib.eud_sq

        #void set_rho_err_k(Datapoint_info* points, FLOAT_TYPE* rho, FLOAT_TYPE* rhoErr, idx_t* k, size_t n)
        self._set_rho_err_k = self.lib.set_rho_err_k
        self._set_rho_err_k.argtypes = [  ct.POINTER(datapoint_info),    
                                        np.ctypeslib.ndpointer(ct_float_t), 
                                        np.ctypeslib.ndpointer(ct_float_t),
                                        np.ctypeslib.ndpointer(ct_idx_t),
                                        ct.c_uint64 ]

        self._verbose = self.lib.verbose
        self._set_verbose = self.lib.set_verbose_output
        self._set_verbose.argtypes = [ct.c_int32]
        #Datapoint_info* _allocDatapoints(idx_t* indeces, float_t* distances, idx_t n, idx_t k)

        self._allocate_datapoints = self.lib.alloc_datapoints
        self._allocate_datapoints.argtypes = [ct_idx_t]
        self._allocate_datapoints.restype  = ct.POINTER(datapoint_info)


        self._import_neighbors_and_distances = self.lib.import_neighbors_and_distances

        #void _import_neighbors_and_distances(Datapoint_info* points, idx_t* indeces, float_t* distances, idx_t n, idx_t k)
        self._import_neighbors_and_distances.argtypes =   [   ct.POINTER(datapoint_info),
                                                            np.ctypeslib.ndpointer(ct_idx_t),
                                                            np.ctypeslib.ndpointer(ct_float_t),
                                                            ct_idx_t,
                                                            ct_idx_t
                                                        ]
        #void export_neighbors_and_distances(Datapoint_info* points, idx_t* dist_indices, float_t* dists, idx_t n, idx_t k)
        self._export_neighbors_and_distances = self.lib.export_neighbors_and_distances
        self._export_neighbors_and_distances.argtypes =   [   ct.POINTER(datapoint_info),
                                                            np.ctypeslib.ndpointer(ct_idx_t),
                                                            np.ctypeslib.ndpointer(ct_float_t),
                                                            ct_idx_t,
                                                            ct_idx_t
                                                        ]

        #void _import_density(Datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
        self._import_density = self.lib.import_density
        self._import_density.argtypes = [   ct.POINTER(datapoint_info),
                                            np.ctypeslib.ndpointer(ct_idx_t),
                                            np.ctypeslib.ndpointer(ct_float_t),
                                            np.ctypeslib.ndpointer(ct_float_t),
                                            ct_idx_t
                                        ]

        self._export_density = self.lib.export_density
        self._export_density.argtypes = [   ct.POINTER(datapoint_info),
                                            np.ctypeslib.ndpointer(ct_idx_t),
                                            np.ctypeslib.ndpointer(ct_float_t),
                                            np.ctypeslib.ndpointer(ct_float_t),
                                            ct_idx_t
                                        ]

        #void compute_avg(Datapoint_info* p, FLOAT_TYPE *va, FLOAT_TYPE* ve, FLOAT_TYPE* vals, FLOAT_TYPE* verr, size_t k, size_t n)
        self._compute_avg = self.lib.compute_avg
        self._compute_avg.argtypes = [  ct.POINTER(datapoint_info),    
                                        np.ctypeslib.ndpointer(ct_float_t), 
                                        np.ctypeslib.ndpointer(ct_float_t),
                                        np.ctypeslib.ndpointer(ct_float_t),
                                        np.ctypeslib.ndpointer(ct_float_t),
                                        ct_idx_t,
                                        ct.c_uint64 ]
        #void export_cluster_assignment(Datapoint_info* points, int* labels, idx_t n)

        self._export_cluster_assignment = self.lib.export_cluster_assignment
        self._export_cluster_assignment.argtypes = [ ct.POINTER(datapoint_info),    
                                                    np.ctypeslib.ndpointer(np.int32), 
                                                    ct.c_uint64 ]


        #void export_borders(clusters* clusters, idx_t* border_idx, float_t* border_den, float_t* border_err)
        self._export_borders = self.lib.export_borders
        self._export_borders.argtypes = [ct.POINTER(clusters),
                                        np.ctypeslib.ndpointer(np.int32),
                                        np.ctypeslib.ndpointer(ct_float_t),
                                        np.ctypeslib.ndpointer(ct_float_t),
                                       ]


        self._id_estimate = self.lib.id_estimate
        self._id_estimate.argtypes = [ct.POINTER(datapoint_info), ct.c_uint64, ct_float_t, ct.c_int32]
        self._id_estimate.restype  =  ct.c_double


        self._compute_density_kstarnn = self.lib.compute_density_kstarnn
        self._compute_density_kstarnn.argtypes = [ct.POINTER(datapoint_info), ct.c_double, ct.c_uint64, ct.c_int32]

        self._PAk = self.lib.PAk
        self._PAk.argtypes = [ct.POINTER(datapoint_info), ct.c_double, ct.c_uint64, ct.c_int32]

        self._compute_correction = self.lib.compute_correction
        self._compute_correction.argtypes = [ct.POINTER(datapoint_info), ct_idx_t, ct.c_double]

        self._H1 = self.lib.Heuristic1
        self._H1.argtypes = [ct.POINTER(datapoint_info), ct.c_uint64, ct.c_int32]
        self._H1.restype = clusters

        self._clusters_allocate = self.lib.clusters_allocate
        self._clusters_allocate.argtypes = [ct.POINTER(clusters), ct.c_int]

        self._blas_in_use = self.lib.blas_are_in_use
        self._blas_in_use.restypes = ct.c_int32

        self._H2 = self.lib.Heuristic2
        self._H2.argtypes = [ct.POINTER(clusters), ct.POINTER(datapoint_info), ct.c_int32]

        self._H3 = self.lib.Heuristic3
        self._H3.argtypes = [ct.POINTER(clusters), ct.POINTER(datapoint_info), ct.c_double, ct.c_int, ct.c_int32]
        
        self._free_datapoints = self.lib.free_datapoint_array
        self._free_datapoints.argtypes = [ct.POINTER(datapoint_info), ct.c_uint64]

        self._free_clusters = self.lib.clusters_free
        self._free_clusters.argtypes = [ct.POINTER(clusters)]


class Data(_dadac_loader):
    def __init__(self, data : np.array, verbose = True):
        """Data object for dadaC library

        Args:
            data (np.array): 2d array to use in searching for clusters 

        Raises:
            TypeError: Raises TypeError if a type different from a matrix is passed 
        """
        super().__init__()
        #initialize class
        if self._use_float32:
            self._ftype = np.float32
            self.data = data.astype(np.float32)
            self.data = np.ascontiguousarray(self.data, dtype = np.float32)
            ct_float_t = ct.c_float
        else:
            #self.data = np.ascontiguousarray(data, dtype =np.float64)
            self._ftype = np.float64
            self.data = data.astype(np.float64)
            self.data = np.ascontiguousarray(self.data, dtype = np.float64)

        if self._use_int32:
            self._itype = np.uint32
        else:
            self._itype = np.uint64
        self.state = {
                "ngbh" : False,
                "id" : False,
                "density" : False,
                "clustering" : False,
                "use_float32": self._use_float32,
                "use_int32": self._use_int32,
                "use_sparse": None,
                "compute_halo": None 
                }
        

        if len(self.data.shape) != 2:
            raise TypeError("Please provide a 2d numpy array")

        if verbose:
            self._verbose = 1
        else:
            self._verbose = 0


        self._datapoints     = None
        self._clusters       = None
        self.n              = self.data.shape[0]
        self.dims           = self.data.shape[1]
        self.k              = None
        self.id             = None

        self._cluster_assignment = None
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
        self._running_in_notebook = self._is_notebook()

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
                return True  # Other type (?)
        except NameError:
            return False

    def compute_neighbors_kdtree(self, k : int):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
            alg (str): default "kd" for kdtree else choose "vp" for vptree
        """
        self.k = k
        #with sys_pipes():
        if self._running_in_notebook:
            with sys_pipes():
                self._datapoints = self._ngbh_search_kdtree(self.data, self.n, self.dims, self.k, self._verbose)
        else:
            self._datapoints = self._ngbh_search_kdtree(self.data, self.n, self.dims, self.k, self._verbose)
        #Datapoint_info* ngbh_search_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        self.state["ngbh"] = True
        self.neighbors = None

    def compute_avg_over_ngbh(self, vals, error,k):
        retVals = np.zeros_like(vals)
        retError = np.zeros_like(error)
        self._compute_avg(self._datapoints, retVals, retError, vals, error, k, self.n)
        return retVals, retError

    def compute_neighbors_vptree(self, k : int, alg="kd"):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
            alg (str): default "kd" for kdtree else choose "vp" for vptree
        """
        self.k = k
        #Datapoint_info* ngbh_search_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        if self._running_in_notebook:
            with sys_pipes():
                self._datapoints = self._ngbh_search_vptree(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, self._eud, self._verbose)
        else:
            self._datapoints = self._ngbh_search_vptree(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, self._eud, self._verbose)

        self.state["ngbh"] = True
        self.neighbors = None

    def compute_neighbors_bruteforce(self, k : int, alg="kd"):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
            alg (str): default "kd" for kdtree else choose "vp" for vptree
        """
        self.k = k
        #Datapoint_info* ngbh_search_vpTree(void* data, size_t n, size_t byteSize, size_t dims, size_t k, float_t (*metric)(void *, void *));
        if self._running_in_notebook:
            with sys_pipes():
                self._datapoints = self._ngbh_search_bruteforce(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, 0, self._verbose)
        else:
            self._datapoints = self._ngbh_search_bruteforce(self.data.ctypes.data, self.n, self.data.itemsize, self.dims, self.k, 0, self._verbose)
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
            self.compute_neighbors_kdtree(k)
            return
        if alg == "vp":
            self.compute_neighbors_vptree(k)
            return
        if alg == "bf":
            self.k = k
            if self.blas:
                self.compute_neighbors_bruteforce(k)
            else:
                if self._verbose:
                    print("dadac implementation not compiled with blas support for euclidean metric optimization")
                    print("if you want to use it consider compiling against a blas library implementation")
                    print("--> Falling back to sklearn brute force")
                t1 = time.monotonic() 
                nn = NearestNeighbors(n_neighbors=k, n_jobs=-1, p = 2, algorithm="brute").fit(self.data)
                dist, ngbh = nn.kneighbors(self.data)
                ngbh = ngbh.astype(self._itype)

                dist = np.ascontiguousarray(dist.astype(self._ftype),dtype = self._ftype)
                ngbh = np.ascontiguousarray(ngbh)
                
                self._datapoints = self._allocate_datapoints(self.n)
                self._import_neighbors_and_distances(self._datapoints,ngbh,dist,self.n, self._itype(k))

                if self._verbose:
                    t2 = time.monotonic()
                    print(f"\tTotal time: {t2 - t1 : .2f}s")

            self.state["ngbh"] = True
            #self.compute_neighbors_bruteforce(k)
            return

    def compute_id_2NN(self,fraction = 0.9):

        """ Compute the intrinsic dimension of the dataset via the TWO Nearest Neighbors method.
            Ref. paper 

        Raises:
            ValueError: Raises value error if neighbors are not computed yet, use `Data.compute_neighbors()` 
        """

        if not self.state["ngbh"]:
            raise ValueError("Please compute Neighbors before calling this function")
        if self._running_in_notebook:
            with sys_pipes():
                self.id = self._id_estimate(self._datapoints,self.n,fraction, self._verbose)
        else:
            self.id = self._id_estimate(self._datapoints,self.n,fraction, self._verbose)
        self.state["id"] = True

    def import_neighbors_and_distances(self,ngbh,dists,dists_are_sq = False):
        if self._datapoints is None:
            self._datapoints = self._allocate_datapoints(self.n)
        k = ngbh.shape[1]
        if dists_are_sq:
            self._import_neighbors_and_distances(self._datapoints,ngbh.astype(self._itype),dists.astype(self._ftype),self.n, np.uint64(k))
        else:
            dists = dists**2
            self._import_neighbors_and_distances(self._datapoints,ngbh.astype(self._itype),dists.astype(self._ftype),self.n, np.uint64(k))

    def import_density(self, log_den,log_den_err,kstar):
        if self._datapoints is None:
            self._datapoints = self._allocate_datapoints(self.n)
        self.state["density"] = True
        #void _import_density(Datapoint_info* points, idx_t* kstar, float_t* density, float_t* density_err, idx_t n)
        self._log_den = log_den.astype(self._ftype)
        self._log_den_err = log_den_err.astype(self._ftype)
        self._kstar = kstar.astype(self._itype)
        self._import_density(self._datapoints,self.kstar,self.log_den, self.log_den_err,self.n)

    def compute_density_kstarNN(self):

        """Compute density value for each point

        Raises:
            ValueError: Raises value error if ID is not computed, use `Data.computeIDtwoNN` method
        """
        if not self.state["id"]:
            raise ValueError("Please compute ID before calling this function")
        if self._running_in_notebook:
            with sys_pipes():
                self._compute_density_kstarnn(self._datapoints, self.id, self.n, self._verbose)
        else:
            self._compute_density_kstarnn(self._datapoints, self.id, self.n, self._verbose)

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
        if self._running_in_notebook:
            with sys_pipes():
                self._PAk(self._datapoints, self.id, self.n, self._verbose)
        else:
            self._PAk(self._datapoints, self.id, self.n, self._verbose)

        self.state["density"] = True
        self.density = None
        self.densityError = None

    def compute_clustering_ADP(self,Z : float, halo = False, use_sparse = "auto"):

        """Compute clustering via the Advanced Density Peak method

        Args:
            Z (float): Z value for the method 
            use_sparse (str): optional [``auto``,True, False], use sparse implementation of border storage between clusters. Memory usage for big datsets is significant. 

        Raises:
            ValueError: Raises value error if density is not computed, use `Data.computeDensity()` method
        """
        if not self.state["density"]:
            raise ValueError("Please compute density before calling this function")
        if use_sparse == "auto":
            if self.n > 2e6:
                self.state["use_sparse"] = True
            else: 
                self.state["use_sparse"] = False 
        elif use_sparse == True:
            self.state["use_sparse"] = True
        else:
            self.state["use_sparse"] = False

        self.state["compute_halo"] = halo 
        self.Z = Z
        if self.state["clustering"]:
            self._free_clusters(self._clusters)
            self._cluster_assignment = None
            self._border_indices = None
            self._log_den_bord = None
            self._log_den_bord_err = None
        if self._running_in_notebook:
            with sys_pipes():
                self._compute_correction(self._datapoints, self.n, self.Z)
                self._clusters = self._H1(self._datapoints, self.n, self._verbose)
                self._clusters_allocate(ct.pointer(self._clusters), 1 if self.state["use_sparse"] else 0)
                self._H2(ct.pointer(self._clusters), self._datapoints, self._verbose)
                self._H3(ct.pointer(self._clusters), self._datapoints, self.Z, 1 if halo else 0, self._verbose )
        else:
            self._compute_correction(self._datapoints, self.n, self.Z)
            self._clusters = self._H1(self._datapoints, self.n, self._verbose)
            self._clusters_allocate(ct.pointer(self._clusters), 1 if self.state["use_sparse"] else 0)
            self._H2(ct.pointer(self._clusters), self._datapoints, self._verbose)
            self._H3(ct.pointer(self._clusters), self._datapoints, self.Z, 1 if halo else 0, self._verbose )
        self.state["clustering"] = True
        self.clusterAssignment = None

    def _get_cluster_assignment(self):


        """Retrieve cluster assignment

        Raises:
            ValueError: Raises error if clustering is not computed, use `Data.computeClusteringADP(Z)` 

        Returns:
            List of cluster labels
            
        """
        if self.state["clustering"]:
            self._cluster_assignment = np.zeros(self.n, np.int32)
            self._export_cluster_assignment(self._datapoints,self._cluster_assignment,self.n)
        else:
            raise ValueError("Clustering is not computed yet")

    @property
    def cluster_assignment(self):
        if self._cluster_assignment is None:
            self._get_cluster_assignment()
        return self._cluster_assignment


    def _get_density(self):

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
            self._export_density(self._datapoints,self._kstar, self._log_den, self._log_den_err, self.n)
        else:
            raise ValueError("Density is not computed yet use `Data.computeDensity()`")

    def _get_N_clusters(self):
        if self.state["clustering"]:
            self._N_clusters = len(set(self.cluster_assignment))
            if self.state["compute_halo"]:
                self._N_clusters -= 1 
        else:
            raise ValueError("Borders not computed yet use `Data.compute_clustering_[DP,ADP,...]()`")


    def _get_borders(self):

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
            self._export_borders(self._clusters,self._border_indices, self._log_den_bord, self._log_den_bord_err)
            #ugly thing to retrieve correct thing in dadapy as border density
            mm = np.zeros_like(self._log_den_bord)
            f = np.where(self._border_indices == -1)
            mm[f] = -1
            mm += np.eye(self._N_clusters)
            self._log_den_bord -= mm
        else:
            raise ValueError("Borders not computed yet use `Data.compute_clustering_[DP,ADP,...]()`")

    def _get_cluster_centers(self):
        if self.state["clustering"]:
            self._cluster_centers = [int(self._datapoints[i].array_idx)  for i in range(self.n) if bool(self._datapoints[i].is_center)]
        else:
            raise ValueError("Clustering not computed yet use `Data.compute_clustering_[DP,ADP,...]()`")

    @property
    def cluster_centers(self):
        if self._cluster_centers is None:
            self._get_cluster_centers()
        return self._cluster_centers

    @property
    def N_clusters(self):
        if self._N_clusters is None:
            self._get_N_clusters()
        return self._N_clusters


    @property
    def log_den(self):
        if self._log_den is None:
            self._get_density()
        return self._log_den

    @property
    def log_den_err(self):
        if self._log_den_err is None:
            self._get_density()
        return self._log_den_err

    @property
    def kstar(self):
        if self._kstar is None:
            self._get_density()
        return self._kstar

    
    @property
    def distances(self): 
        if self._distances is None:
            self._distances = np.zeros((self.n,self.k), self._ftype)
            self._dist_indices = np.zeros((self.n,self.k), self._itype)
            self._export_neighbors_and_distances(self._datapoints,self._dist_indices, self._distances, self.n, self.k)
            self._distances = np.sqrt(self._distances)
        return self._distances

    @property
    def dist_indices(self): 
        if self._dist_indices is None:
            self._distances = np.zeros((self.n,self.k), self._ftype)
            self._dist_indices = np.zeros((self.n,self.k), self._itype)
            self._export_neighbors_and_distances(self._datapoints,self._dist_indices,self._distances, self.n, self.k)
            self._distances = np.sqrt(self._distances)
        return self._dist_indices

    @property
    def log_den_bord(self):
        if self._log_den_bord is None:
            self._get_borders()
        return self._log_den_bord
    
    @property
    def log_den_bord_err(self):
        if self._log_den_bord_err is None:
            self._get_borders()
        return self._log_den_bord_err

    @property
    def border_indices(self):
        if self._border_indices is None:
            self._get_borders()
        return self._border_indices
    
    #def set_rho_err_k(self, rho, err, k):
    #    self._set_rho_err_k(self._datapoints, rho, err, k, self.n)
    #    self.state["density"] = True
    #    return

    def __del__(self):
        if not self._datapoints is None:
            self._free_datapoints(self._datapoints, self.n)
        if not self._clusters is None:
            self._free_clusters(ct.pointer(self._clusters))
