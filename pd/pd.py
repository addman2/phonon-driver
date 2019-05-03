import os
import time
import datetime
import textwrap
import numpy                                        as np
from   ase                     import Atoms
from   ase.io                  import read
from   ase.db                  import connect
from   ase.dft.kpoints         import bandpath
from   ase.calculators.vasp    import Vasp
from   phonopy                 import Phonopy
from   pyspglib                import spglib
from   pd.tools                import get_ideal_bandpath

class Pd():

    def __init__(self, name):
        self._name     = name
        self._cn       = connect("{}.db".format(name))
        self._pressure = 0
        self._phonon   = None
        self._load_settings()
        self._prepare_folders()
        try: ret = self._cn.get(Initial = True)
        except:
            self._cn.write(self._xtal, Initial = True)

    def _logthis(self,s):
        with open("pd.log","a") as f:
            ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
            tw = textwrap.TextWrapper()
            for i, line in enumerate(tw.wrap(s)):
                if i == 0:
                    f.write(ts + ": " + line + "\n")
                else:
                    f.write(" " * (len(ts) + 2) + line + "\n")

    def _prepare_folders(self):
        dirs = [ "data",
                 "figures",
                 self._main["vasp_wd"]]
        for dir in dirs:
            os.makedirs( dir, exist_ok = True )

    def _load_settings(self):
        for d, dn, fn in os.walk("settings"):
            sf = [ f for f in fn if f.endswith("_settings.py") ]
            break
        for f in sf:
            global settings
            self._logthis (" Loading {} ...".format(f))
            with open("{}/{}".format(d,f),"rb") as fi:
                exec(fi.read(),globals())
                self.__dict__["_"+f.replace("_settings.py","")] = settings
                del settings

    def _get_similar_crystal(self, pressure):
        xs = list(self._cn.select(Basic = True))
        if len(xs) == 0: return False, None
        xs = sorted(xs, key = lambda x: (x.Pressure - pressure)**2)
        optimized = ( xs[0].Pressure - pressure == 0 )
        return optimized, xs[0].toatoms()

    def _get_crystal(self):
        ret = None

        try: ret = self._cn.get(Initial = True)
        except: pass
        try: ret = self._cn.get(Pressure = self._pressure,
                                Basic    = True)
        except: pass

        if not ret is None:
            ret = ret.toatoms()

        return ret

    def _get_calculator(self, purpose):
        ret = None
        if purpose == "optimize":
            settings = self._common.copy()
            settings.update(self._opt)
            if not self._pressure is None:
                settings["pstress"] = self._pressure
            ret = Vasp(**settings)
        if purpose == "super_optimize":
            settings = self._common.copy()
            settings.update(self._sopt)
            if not self._pressure is None:
                settings["pstress"] = self._pressure
            ret = Vasp(**settings)
        if purpose == "ftest":
            settings = self._common.copy()
            settings.update(self._ftest)
            ret = Vasp(**settings)
        return ret

    def _get_energy(self, crystal, run):

        if run["Aurel"][0] == "X":
            os.environ["VASP_COMMAND"] = self._main["aurrun"]
        if run["Local"][0] == "X":
            os.environ["VASP_COMMAND"] = self._main["locrun"]

        cwd = os.getcwd()
        try: shutil.rmtree(self._main["vasp_wd"]+"/*")
        except: pass
        os.chdir(self._main["vasp_wd"])
        ret = crystal.get_potential_energy()
        os.chdir(cwd)

        return ret

    def _check_forces(self, crystal, fmax, run):
        calculator = self._get_calculator("ftest")
        crystal.set_calculator(calculator)
        E = self._get_energy(crystal, run)
        forces = crystal.get_forces()
        normforces = np.linalg.norm(forces, axis=1)
        #msg = ""
        for ii in range(len(forces)):
            #msg += (3*"{:9.6f}" + " | {:9.6f}\n").format(*(forces[ii]),normforces[ii])
            self._logthis((3*"{:9.6f}" + " | {:9.6f}\n").format(*(forces[ii]),normforces[ii]))
        #self._logthis(msg)
        if normforces.all() < fmax:
            pass
        return normforces.all() < fmax

    def optimize_to_pressure(self, pressure, run):
        self._pressure = pressure
        optimized, x = self._get_similar_crystal(pressure)
        if optimized: return
        if x is None: x = self._get_crystal()
        calculator = self._get_calculator("optimize")
        x.set_calculator(calculator)
        E = self._get_energy(x, run)
        H = E + pressure * x.get_volume() / 1602.176627
        HoA = H/len(x)
        spg = spglib.get_spacegroup(x)
        self._cn.write(x,
                       Basic      = True,
                       Pressure   = pressure,
                       SpaceGroup = spg)

    def prepare_supercell(self, repeat, fmax, run):
        try:
            self._cn.get(Pressure = self._pressure,
                         Supercell = True,
                         Displaced = False)
            return
        except: pass
        x = self._get_crystal().repeat(repeat)
        while not self._check_forces(x, fmax, run):
            self._logthis("Optimizing {}".format(self._pressure))
            calculator = self._get_calculator("super_optimize")
            x.set_calculator(calculator)
            E = self._get_energy(x, run)
        spg = spglib.get_spacegroup(x)
        self._cn.write(x,
                       Supercell  = True,
                       Displaced  = False,
                       Pressure   = self._pressure,
                       SpaceGroup = spg)

    def calculate_phonons(self, run):
        self._prepare_phonons()
        set_of_forces = self._calculate_forces(run)
        self._phonon.produce_force_constants( forces = set_of_forces )

        omega, dos = self._calculate_dos()
        band_data  = self._calculate_bands()

        data = { "Dos"   : [ omega, dos ],
                 "Bands" : band_data }

        row = self._cn.get(Basic    = True,
                           Pressure = self._pressure)

        id = row.id
        self._cn.update(id, data = data)

    def _prepare_phonons(self):
        row = self._cn.get(Pressure  = self._pressure,
                           Supercell = True,
                           Displaced = False)
        x = row.toatoms()
        phonon = Phonopy(x,
                         np.eye(3),
                         primitive_matrix="auto")

        phonon.generate_displacements(distance=self._phpy["h_step"])
        #disps = phonon.get_displacements()
        phonon.supercellrow = row
        self._phonon = phonon

    def _calculate_forces(self, run):
        try: return self._phonon.supercellrow.data["Forces"]
        except: pass
        supercells = self._phonon.get_supercells_with_displacements()
        set_of_forces = []

        for ii, scell in enumerate(supercells):
            x = Atoms(symbols=scell.get_chemical_symbols(),
                         scaled_positions=scell.get_scaled_positions(),
                         cell=scell.get_cell(),
                         pbc=True)

            calculator = self._get_calculator("ftest")
            x.set_calculator(calculator)
            E = self._get_energy(x, run)
            self._cn.write(x,
                           Supercell = True,
                           Displaced = True,
                           Pressure  = self._pressure,
                           Displnumb = ii)

            forces = x.get_forces()

            drift_force = forces.sum(axis=0)

            for force in forces:
                force -= drift_force / forces.shape[0]

            set_of_forces.append(forces)

        self._cn.update(self._phonon.supercellrow.id,
                        data = { "Forces" : set_of_forces })
        return set_of_forces

    def _calculate_dos(self):

        self._phonon.set_mesh(self._dos["kpts"],
                              is_eigenvectors=True)
        qpoints, weights, frequencies, eigvecs = self._phonon.get_mesh()
        self._phonon.set_total_DOS(tetrahedron_method=self._dos["tetrahedron"])

        omega, dos = np.array(self._phonon.get_total_DOS())

        return omega, dos

    def _calculate_bands(self):

        if self._bands["path"].lower() == "auto":
            row = self._cn.get(Pressure  = self._pressure,
                               Supercell = True,
                               Displaced = False)
            x = row.toatoms()
            self._bands["path"] = get_ideal_bandpath(x)


        path = self._bands["path"]
        try: path = self._bands["pathp"]
        except: pass

        kpts, x, X = bandpath(path,
                              np.identity(3),
                              self._bands["npoints"])
        self._phonon.set_band_structure([kpts])
        _, __, bands, ___ = self._phonon.get_band_structure()

        band_data = {"path"  : self._bands["path"],
                     "kpts"  : kpts,
                     "x"     : x,
                     "X"     : X,
                     "bands" : bands }

        return band_data

    def show_phonon(self, pressure = 0):
        self.optimize_to_pressure(pressure, None)
        self.prepare_supercell(None, None, None)
        self.calculate_phonons(None)

        return self._phonon

