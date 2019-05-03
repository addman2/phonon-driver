#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser
from ase import Atoms
from pd.pd import *

def write(mod, name):
    xtal = Atoms(mod.symbols,
                 cell = mod.cell,
                 positions = mod.get_positions())
    xtal.write(name)

def work():
    parser = ArgumentParser(description = "Tool producing modulated supercells")
    parser.add_argument("Database",
                        metavar="DATABASE",
                        type=str,)
    parser.add_argument("-q",
                        type=float,
                        dest="q",
                        metavar="Q",
                        nargs=3,
                        default=[0,0,0],
                        help="Q-Point")
    parser.add_argument("-d",
                        type=int,
                        dest="d",
                        metavar="Q",
                        nargs=3,
                        default=[2,2,2],
                        help="Supercell dimension")
    parser.add_argument("-i",
                        type=int,
                        dest="index",
                        metavar="I",
                        default=0,
                        help="Mode index. If you want to write all modulations i = -1")
    parser.add_argument("-a",
                        type=float,
                        dest="ampl",
                        metavar="AMP",
                        default=0.5,
                        help="Amplitude (in Angstorm)")
    parser.add_argument("-s",
                        type=float,
                        dest="pressure",
                        metavar="Pressure",
                        default=0.0,
                        help="Pressure")
    parser.add_argument("-p",
                        type=float,
                        dest="phase",
                        metavar="PHASE",
                        default=0.0,
                        help="Phase")
    parser.add_argument("--dest",
                        type=str,
                        dest="dest",
                        metavar="DESTINATION",
                        default="MODULATED_POSCAR",
                        help="Destination file")

    args = parser.parse_args()

    db = args.Database
    if db.endswith(".db"):
        if not os.path.isfile(db):
            print("No such database {}".format(db))
            sys.exit()
        db = db[:-3]
    else:
        if not os.path.isfile(db+".db"):
            print("No such database {}.db".format(db))
            sys.exit()
    p = Pd(db)
    x = p.show_phonon(pressure = args.pressure)

    freqs = x.get_frequencies(args.q)

    if args.index == -1:
        numof = len(freqs)
        mods = []
        for ii in range(numof):
            mods.append([args.q, ii, args.ampl, args.phase])
    else:
        mods = [[ args.q, args.index, args.ampl, args.phase ]]

    x.set_modulations(args.d, mods)
    modulateds = x.get_modulated_supercells()

    if len(modulateds) == 1:
        print("Frequency: {:5.3} THz".format(freqs[args.index]))
        write(modulateds[0], args.dest)
    else:
        for ii, mod in enumerate(modulateds):
            print("Frequency i = {}: {:5.3f} THz".format(ii, freqs[ii]))
            write(mod, args.dest + ".{}".format(ii))


if __name__ == "__main__":
    work()
