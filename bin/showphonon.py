#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser
from ase import Atoms
from pd.pd import *

def work():
    parser = ArgumentParser(description = "Tool producing modulated supercells")
    parser.add_argument("Database",
                        metavar="DATABASE",
                        type=str,)
    parser.add_argument("-q",
                        type=int,
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
                        help="Mode index")
    parser.add_argument("-a",
                        type=float,
                        dest="ampidute",
                        metavar="AMP",
                        default=0.5,
                        help="Amplitude (in Angstorm)")
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
    x = p.show_phonon(None)

    x.set_modulations(args.d, [[ args.q, args.index, args.index, args.phase ]])
    modulated = x.get_modulated_supercells()[0]

    xtal = Atoms(modulated.symbols,
                 cell = modulated.cell,
                 positions = modulated.get_positions())
    xtal.write(args.dest)

if __name__ == "__main__":
    work()
