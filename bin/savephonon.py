#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser
from ase import Atoms
from pd.pd import *

def work():
    parser = ArgumentParser(description = "Saves phonopy data")
    parser.add_argument("Database",
                        metavar="DATABASE",
                        type=str,)
    parser.add_argument("-s",
                        type=float,
                        dest="pressure",
                        metavar="Pressure",
                        default=0.0,
                        help="Pressure")

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

    x.save()

if __name__ == "__main__":
    work()
