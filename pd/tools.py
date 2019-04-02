import os
import sys
from   ase.db           import connect
from   matplotlib       import pyplot           as plt
from   matplotlib       import gridspec         as gridspec

def get_ideal_bandpath(xtal):

    bp = { "cub" : "GXMGRX,MR",
           "bcc" : "GHNGPH,PN",
           "fcc" : "GXWKGLUWLK,UX"
           "hex" : "GMKGALHA,LM,KH"
           "tet" : "GXMGZRAZ,XR,MA" }

    sd = spglib.get_symmetry_dataset(xtal)
    number = sd["number"]
    spacegroup = sd["spacegroup"]

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

    plt.savefig("{}/{:05d}.png".format("figures",row.Pressure))

if __name__ == "__main__":

    cn = connect(sys.argv[1])

    for row in cn.select(Basic = True):
        plot_phonon(row)
