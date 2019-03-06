import os
import sys
from   ase.db           import connect
from   matplotlib       import pyplot           as plt
from   matplotlib       import gridspec         as gridspec

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
