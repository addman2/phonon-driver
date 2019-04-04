import os
import sys
import pkgutil
import pd.examples
from   importlib        import import_module
from   ase.db           import connect
from   pyspglib         import spglib
from   matplotlib       import pyplot           as plt
from   matplotlib       import gridspec         as gridspec

def get_examples():
    ret = []
    package = pd.examples
    for importer, modname, ispkg in pkgutil.walk_packages(path=package.__path__, 
                                                          prefix=package.__name__ + ".",
                                                          onerror=lambda x: None):
        ret.append(modname.split(".")[-1])
    return ret

def make_example(name, folder):
    module_name = "pd.examples.{}".format(name)
    m = import_module(module_name)

    os.makedirs("{}/{}/".format(folder,"settings"), exist_ok = True)

    with open("{}/run.py".format(folder), "wb") as fo:
        fo.write(m.run)
    
    for s in m.settings:
        with open("{}/settings/{}".format(folder, s), "w") as fo:
            fo.write(m.settings[s])

def get_ideal_bandpath(xtal):

    bp = { "cub" : "GXMGRX,MR",
           "bcc" : "GHNGPH,PN",
           #"fcc" : "GXWKGLUWLK,UX",
           "fcc" : "GXMG",
           "hex" : "GMKGALHA,LM,KH",
           "tet" : "GXMGZRAZ,XR,MA" }

    sd = spglib.get_symmetry_dataset(xtal)
    number = sd["number"]
    spacegroup = sd["international"]

    if number in range(195, 230 + 1):
        # Cubic spacegroup
        if spacegroup[0] == "P":
            # Primitive
            return bp["cub"]
        if spacegroup[0] == "I":
            # Base Center Cubic
            return bp["bcc"]
        if spacegroup[0] == "F":
            # Fase Center Cubic
            return bp["fcc"]

    if number in range(168, 194 + 1):
        # Hexagonal spacegroup
        return bp["hex"]

    if number in range(75, 142 + 1):
        # Tetragonal spacegroup
        if spacegroup[0] == "P":
            # Primitive
            return bp["tet"]

    raise Exception("Cannot dermine ideal band path for spacegroup {}".format(spacegroup))

def plot_point(rows, point):

    data = []
    for row in rows:
        p = row.Pressure
        b = row.data["Bands"]
        points = [ pt for pt in b["path"] if pt != "," ]
        ind_point = points.index(point)
        coo_point = b["X"][ind_point]
        ind_coord = list(b["x"]).index(coo_point)
        bands = b["bands"][0][ind_coord]
        data.append([p,bands])

    data = sorted(data, key = lambda x: x[0])

    bands = [ b for p, b in data ]
    ps    = [ p for p, b in data ]

    gs = gridspec.GridSpec(1,10)
    ax1 = plt.subplot(gs[0,:9])

    ax1.plot(ps,
             bands, "b-", lw=1.8)

    #ylim = list(ax1.get_ylim())
    ##ylim[0] = 0.0

    #for X in b["X"]:
    #    ax1.plot( [X]*2,ylim,"k-", lw=2.0)

    #ax1.plot([min(b["x"]),max(b["x"])],[0,0],"r")

    #ax1.set_xlim([min(b["x"]),max(b["x"])])
    #ax1.set_ylim(ylim)
    #ax1.set_xticks(b["X"])
    #points = [ p for p in b["path"] if p != "," ]
    #ax1.set_xticklabels([x.replace("Gamma","$\Gamma$") for x in points])

    #ax1.set_xlabel("q - path")
    #ax1.set_ylabel("[THz]")

    #ax2.plot(d[1],d[0],"b-")

    #plt.tight_layout()
    #os.makedirs("figures",exist_ok=True)

    name = "name"
    try: name = row.name
    except: pass

    ax1.set_title("Energies in {} - point".format(point))

    figname = "{}/{}.{}.png".format("figures",
                                    name,
                                    "point-{}".format(point))
    plt.savefig(figname)

def plot_phonon(row):

    d = row.data["Dos"]
    b = row.data["Bands"]

    gs = gridspec.GridSpec(1,10)
    ax1 = plt.subplot(gs[0,:7])
    ax2 = plt.subplot(gs[0,7:],sharey=ax1)

    ax1.plot(b["x"],
             b["bands"][0], "b-", lw=1.8)

    ylim = list(ax1.get_ylim())
    #ylim[0] = 0.0

    for X in b["X"]:
        ax1.plot( [X]*2,ylim,"k-", lw=2.0)

    ax1.plot([min(b["x"]),max(b["x"])],[0,0],"r")

    ax1.set_xlim([min(b["x"]),max(b["x"])])
    ax1.set_ylim(ylim)
    ax1.set_xticks(b["X"])
    points = [ p for p in b["path"] if p != "," ]
    ax1.set_xticklabels([x.replace("Gamma","$\Gamma$") for x in points])

    ax1.set_xlabel("q - path")
    ax1.set_ylabel("[THz]")

    ax2.plot(d[1],d[0],"b-")

    plt.tight_layout()
    os.makedirs("figures",exist_ok=True)

    name = "name"
    try: name = row.name
    except: pass

    figname = "{}/{}.{}.{:05d}.png".format("figures",
                                           name,
                                           "bands",
                                           row.Pressure)
    print("Saving {}".format(figname))
    plt.savefig(figname)
