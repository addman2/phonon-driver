import os
import numpy                                        as np
from   ase                     import Atoms
from   ase.io                  import read
from   ase.db                  import connect
from   ase.dft.kpoints         import bandpath
from   ase.calculators.vasp    import Vasp
from   phonopy                 import Phonopy

class Pd():

    def __init__(self, name):
        self._name     = name
        self._cn       = connect("{}.db".format(name))
        self._pressure = 0
        self._crystal  = None
        self._phonon   = None
        self._prepare_for_phonon = False
        self._load_settings()
        self._prepare_folders()

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
            print(" Loading {} ...".format(f))
            with open("{}/{}".format(d,f),"rb") as fi:
                exec(fi.read(),globals())
                self.__dict__["_"+f.replace("_settings.py","")] = settings
                del settings

    def _get_crystal(self, pressure = None):
        ret = None
        if pressure is None:
            ret = self._crystal
        else:
            try:
                row = self._cn.get(Pressure = pressure,
                                   Basic    = True )
                print("Reading ...")
                ret = row.toatoms()
            except: pass
        if ret is None:
            if not self._crystal is None:
                return self._crystal
            if type(self._xtal) is str:
                self._crystal = read(self._xtal)
            else:
                self._crystal = Atoms(**self._xtal)
            self._cn.write(self._crystal, Initial = True)
            return self._crystal
        self._crystal = ret
        return ret

    def _get_calculator(self, purpose, pressure = None):
        ret = None
        if purpose == "optimize":
            settings = self._common.copy()
            settings.update(self._opt)
            if not pressure is None:
                settings["pstress"] = pressure
            ret = Vasp(**settings)
        if purpose == "super_optimize":
            settings = self._common.copy()
            settings.update(self._sopt)
            if not pressure is None:
                settings["pstress"] = pressure
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
        forces = np.linalg.norm(forces, axis=1)
        return forces.any() < fmax

    def optimize_to_pressure(self, pressure, run):
        self._pressure = pressure
        calculator = self._get_calculator("optimize")
        x = self._get_crystal(pressure)
        x.set_calculator(calculator)
        E = self._get_energy(x, run)
        H = E + pressure * x.get_volume() / 1602.176627
        HoA = H/len(x)
        self._cn.write(x,
                       Basic    = True,
                       Pressure = pressure)
        self._crystal = x
        self._prepare_for_phonon = False

    def prepare_supercell(self, repeat, fmax, run):
        x = self._crystal.repeat(repeat)
        while not self._check_forces(x, fmax, run):
            print("Optimizing {}".format(self._pressure))
            calculator = self._get_calculator("super_optimize", pressure = self._pressure)
            x.set_calculator(calculator)
            E = self._get_energy(x, run)
        self._cn.write(x,
                       Supercell = True,
                       Displaced = False,
                       Pressure  = self._pressure)

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
        self._phonon = phonon

    def _calculate_forces(self, run):
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
        return set_of_forces

    def _calculate_dos(self):

        self._phonon.set_mesh(self._dos["kpts"],
                              is_eigenvectors=True)
        qpoints, weights, frequencies, eigvecs = self._phonon.get_mesh()
        self._phonon.set_total_DOS(tetrahedron_method=self._dos["tetrahedron"])

        omega, dos = np.array(self._phonon.get_total_DOS())

        return omega, dos

    def _calculate_bands(self):

        kpts, x, X = bandpath(self._bands["path"],
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
