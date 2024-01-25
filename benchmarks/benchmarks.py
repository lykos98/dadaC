import dadapy
import dadac
import time
from tabulate import tabulate
import numpy as np
from sklearn.datasets import fetch_openml

# change this to run tests on bigger datasets
run_big = False

# change this to affect ngbh search and Z parameter of ADP
k = 300
Z = 3

base_url = (
    "https://raw.githubusercontent.com/sissa-data-science/DADApy/main/examples/datasets"
)
example_names = [
    "Fig1.dat",
    "Fig2.dat",
    "Fig1_mobius.dat",
    "FigS1.dat",
    "FigS2.dat",
    "FigS3.dat",
    "FigS4.dat",
]

dsInfo = []
results = []

# -----------------------------------------
# Change this directory to your saving path
# -----------------------------------------
res_f_name = ""


# ------- UTILITY FUNCTIONS ---------------
def get_example_from_drive():
    import os

    url = "https://drive.usercontent.google.com/download?id=1zL39-F7PZvoYbZihMy5jAGuT3OTPbDPv&export=download&authuser=1&confirm=t&uuid=207b1e31-c1df-46ed-8ea9-fdfae7b3ae77&at=APZUnTVuuHrrOGJaQWYCYGQEDI1i:1701945274700"

    if not os.path.exists("/tmp/example.npy"):
        cmd = f'wget --no-check-certificate "{url}" -O /tmp/example.npy'
        os.system(cmd)

    data = np.fromfile("/tmp/example.npy", np.float32)
    data = data.reshape((data.shape[0] // 5, 5))
    return data


def getFromUrl(base_url, ds_name, sep=" "):
    import urllib

    url = "".join((base_url, "/", ds_name))
    data = []
    print(f"Downloading {ds_name} from dadapy repository")
    for line in urllib.request.urlopen(url):
        data.append(
            [
                float(n.replace("E", "e"))
                for n in line.decode().split(sep)
                if len(n) > 0 and n != "\n"
            ]
        )

    data = np.array(data, dtype=np.float64)
    return data


def stringifyNum(n):
    s = ""
    nn = n

    if n > 10**6:
        s = "G"
        nn = n / (10**9)
        return f"{nn:.1f}{s}"

    if n > 10**3:
        s = "k"
        nn = n / (10**3)
        return f"{nn:.1f}{s}"

    if n > 10**6:
        s = "M"
        nn = n / (10**6)
        return f"{nn:.1f}{s}"


def profileAndRun(dataset, dataset_name, k, Z, results, halo=False):
    dp = dadapy.Data(dataset, verbose=True)
    dc = dadac.Data(dataset)

    pres = []

    print(f"**** dadapy on {dataset_name}")
    print("dadapy output\n", "-" * (30), "\n")

    pres.append(["Method", "part", "time"])

    dsInfo.append(
        f"{dataset_name} N = {stringifyNum(dataset.shape[0])} D = {dataset.shape[1]}"
    )
    t1 = time.monotonic()

    dp.compute_distances(k)
    dp.compute_id_2NN()
    dp.compute_density_kstarNN()

    t2 = time.monotonic()
    pres.append(["py", "ngbh and density", f"{t2 - t1: .2f}s"])

    t1 = time.monotonic()

    dp.compute_clustering_ADP(Z=Z, halo=halo)

    t2 = time.monotonic()
    pres.append(["py", "ADP", f"{t2 - t1: .2f}s"])

    print(f"\n**** dadac on {dataset_name} ****")
    print("dadapy output\n", "-" * (30), "\n")

    t1 = time.monotonic()
    dc.compute_distances(k)

    dc.compute_id_2NN()
    dc.compute_density_kstarNN()

    t2 = time.monotonic()
    pres.append(["C", "ngbh and density", f"{t2 - t1: .2f}s"])

    t1 = time.monotonic()

    dc.compute_clustering_ADP(Z=Z, halo=halo)

    t2 = time.monotonic()
    pres.append(["C", "ADP", f"{t2 - t1: .2f}s"])

    results.append(pres)

    print("\n**** Comparing Results ****")
    c1 = dp.cluster_assignment
    c2 = dc.cluster_assignment
    errors = np.where(c1 != c2)[0].shape[0]
    print(f" --> \t Found {errors} errors! \n")


# ---------------------------------------------------------------------
#
#              Benchmarks and tests dadac vs dadapy
#
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
#
# Gaussian mixture 2D
#
# ---------------------------------------------------------------------


n = 1000
x1 = np.random.normal([0, 2], 1, size=(n, 2))
x2 = np.random.normal([2, 0], 1, size=(n, 2))
x = np.concatenate([x1, x2])

profileAndRun(x, "2D Gaussian", k, Z, results)

# ---------------------------------------------------------------------
#
# Gaussian mixture 5D
#
# ---------------------------------------------------------------------

n = 50000
x1 = np.random.normal([0, 2, 0, 0, 0], 1, size=(n, 5))
x2 = np.random.normal([2, 0, 0, 0, 0], 1, size=(n, 5))
x = np.concatenate([x1, x2])

profileAndRun(x, "5D Gaussian", k, Z, results)

# ---------------------------------------------------------------------
#
# MNIST
#
# ---------------------------------------------------------------------

x, y = fetch_openml(
    "mnist_784", version=1, return_X_y=True, as_frame=False, parser="pandas"
)

x = np.ascontiguousarray(x)
x = x.astype(np.float64) / 255.0


profileAndRun(x, "MNIST", k, Z, results)

# ---------------------------------------------------------------------
#
# Dadapy examples
#
# ---------------------------------------------------------------------

# trying to import from dadapy examples

for ds_name in example_names[:2]:
    data = getFromUrl(base_url, ds_name)
    profileAndRun(data, ds_name, k, Z, results)

# ---------------------------------------------------------------------
#
# Cosmological simulation data
#
# ---------------------------------------------------------------------

data = get_example_from_drive()
profileAndRun(data[:100000], "Astro (sub)Set 1", k, Z, results)

if run_big:
    profileAndRun(data[:500000], "Astro (sub)Set 1", k, Z, results)
    profileAndRun(data, "Astro (sub)Set 1", k, Z, results)


# ---------------------------------------------------------------------
#
# Generate tables
#
# ---------------------------------------------------------------------

print("Printing results \n")

for header, tab in zip(dsInfo, results):
    print(header)
    print(tabulate(tab, headers="firstrow", tablefmt="fancy_grid"))

with open(res_f_name, "w") as f:
    for header, tab in zip(dsInfo, results):
        f.write(f"{header}\n")
        tab = tabulate(tab, headers="firstrow", tablefmt="github", floatfmt=".2f")

        f.write(f"{tab}\n")


# print(results)
