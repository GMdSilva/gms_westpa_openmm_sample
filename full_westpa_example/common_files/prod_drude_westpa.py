from sys import stdout
import json
import pprint

from openmm import *
from openmm.app import *
from openmm.unit import *


def read_configs_json(filename):
    with open(filename) as f:
        prod_configs = json.load(f)
    return prod_configs


def get_cubic_box(psf, rstfile):
    """
    Gets the box dimensions from a OpenMM .rst file and updates the psf topology object
    """
    f = open(rstfile, 'r')
    box = {}

    while True:
        line = f.readline()
        if not line:
            break

        if line.split()[0] == "<A":
            size = line.split()[1].strip('x="')
            box['A'] = float(size)
        elif line.split()[0] == "<B":
            size = line.split()[2].strip('y="')
            box['B'] = float(size)
        elif line.split()[0] == "<C":
            size = line.split()[3].strip('z="').strip('"/>')
            box['C'] = float(size)
        else:
            pass

    boxX = box['A'] * nanometer
    boxY = box['B'] * nanometer
    boxZ = box['C'] * nanometer

    psf.setBox(boxX, boxY, boxZ)

    return psf


def read_toppar(filename, path):
    extlist = ['rtf', 'prm', 'str']

    parFiles = ()
    for line in open(filename, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0:
            ext = parfile.lower().split('.')[-1]
            if not ext in extlist:
                continue
            parFiles += (path + parfile,)

    params = CharmmParameterSet(*parFiles)
    return params


class WEDrudeProd:
    def __init__(self, config_file):
        self.configs = read_configs_json(config_file)
        self.psf = None
        self.pdb = None
        self.rst = None
        self.crd = None
        self.charmm_params = None
        self.temperature = None
        self.pressure = None
        self.psf_box = None
        self.system = None
        self.integrator = None
        self.simulation = None
        self.nsteps = None
        self.nslog = None
        self.nsxyz = None

    def setup_parameters(self):
        self.psf = os.path.join(self.configs["base_dir"], "common_files", str(self.configs["psf_filename"]))
        self.pdb = self.configs["parent_prefix"] + ".pdb"  # westpa_scripts/runseg.sh handles moving to each segment dir
        self.rst = self.configs["parent_prefix"] + ".rst"

        toppar = os.path.join(self.configs["base_dir"], "common_files", str(self.configs["toppar_filename"]))

        self.charmm_params = read_toppar(toppar, os.path.join(self.configs["base_dir"], "common_files"))

        prod_tau = self.configs["prod_tau"] * picosecond
        integration_ts = self.configs['integration_ts'] * picosecond
        log_save_freq = self.configs['log_save_freq'] * picosecond
        xyz_save_freq = self.configs['xyz_save_freq'] * picosecond

        self.temperature = self.configs["temperature"] * kelvin
        self.pressure = self.configs["pressure"] * bar

        self.nsteps = int(
            prod_tau.value_in_unit(picosecond) / integration_ts.value_in_unit(picosecond))
        self.nslog = int(
            log_save_freq.value_in_unit(picosecond) / integration_ts.value_in_unit(picosecond))
        self.nsxyz = int(
            xyz_save_freq.value_in_unit(picosecond) / integration_ts.value_in_unit(picosecond))

        top = CharmmPsfFile(self.psf)
        self.crd = PDBFile(self.pdb)
        self.psf_box = get_cubic_box(top, self.rst)

    def load_system(self):
        system_path = os.path.join(self.configs["base_dir"], "common_files", self.configs['system_path'])
        with open(system_path, 'r') as file:
            xml_content = file.read()
        self.system = XmlSerializer.deserialize(xml_content)

    def create_system(self):
        self.system = self.psf_box.createSystem(self.charmm_params,
                                           nonbondedMethod=PME,
                                           nonbondedCutoff=self.configs["nonbondedCutoff"] * nanometer,
                                           switchDistance=self.configs["switchDistance"] * nanometer,
                                           ewaldErrorTolerance=self.configs["ewaldErrorTolerance"],
                                           constraints=HBonds)

        # Sets the value for the DrudeFF Particles
        nbforce = [self.system.getForce(i) for i in range(self.system.getNumForces()) if
                   isinstance(self.system.getForce(i), NonbondedForce)][0]

        nbforce.setNonbondedMethod(NonbondedForce.PME)
        nbforce.setEwaldErrorTolerance(self.configs["ewaldErrorTolerance"])
        nbforce.setCutoffDistance(self.configs["nonbondedCutoff"] * nanometer)
        nbforce.setUseSwitchingFunction(True)
        nbforce.setSwitchingDistance(self.configs["switchDistance"] * nanometer)

        # not every system has NBFIX terms, so check
        cstnb = [self.system.getForce(i) for i in range(self.system.getNumForces()) if
                 isinstance(self.system.getForce(i), CustomNonbondedForce)]

        if cstnb:
            nbfix = cstnb[0]
            nbfix.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
            nbfix.setCutoffDistance(self.configs["nonbondedCutoff"] * nanometer)
            nbfix.setUseSwitchingFunction(True)
            nbfix.setSwitchingDistance(self.configs["switchDistance"] * nanometer)

        self.system.addForce(MonteCarloBarostat(self.pressure, self.temperature))

    def set_integrator(self):
        self.integrator = DrudeLangevinIntegrator(self.temperature,  # temperature
                                                  5 / picosecond,  # frictionCoeff
                                                  1 * kelvin,  # drudeTemperature
                                                  20 / picosecond,  # drudeFrictionCoeff
                                                  self.configs["integration_ts"] * picoseconds)  # stepSize

        self.integrator.setMaxDrudeDistance(0.02)  # Drude Hardwall

    def build_simulation_obj(self):
        self.integrator.setRandomNumberSeed(RAND) # this gets replaced by the WESTPA random seed after westpa_scripts/runseg.sh is ran
        platform = Platform.getPlatformByName('CUDA')
        process_id = os.environ["WM_PROCESS_INDEX"]  # set by the WESTPA work manager
        properties = {"Precision': 'mixed"}
        properties["DeviceIndex"] = process_id  # use a single GPU

        self.simulation = Simulation(self.psf_box.topology, self.system, self.integrator, platform, properties)
        self.simulation.context.setPositions(self.crd.positions)
        self.simulation.context.computeVirtualSites()

    def simulate(self):
        out_log_file = self.configs["child_prefix"] + ".log"
        out_rst_file = self.configs["child_prefix"] + ".rst"
        out_pdb_file = self.configs["child_prefix"] + ".pdb"
        out_dcd_file = self.configs["child_prefix"] + ".dcd"

        with open(self.rst, 'r') as f:
            self.simulation.context.setState(XmlSerializer.deserialize(f.read()))
            currstep = 0  # need to zero step and time since this is a new iteration
            currtime = 0
            self.simulation.currentStep = currstep
            self.simulation.context.setTime(currtime)

        dcd = DCDReporter(out_dcd_file, self.nslog)
        dcd._dcd = DCDFile(dcd._out, self.simulation.topology, self.simulation.integrator.getStepSize(), 0,
                           self.nsxyz)

        self.simulation.reporters.append(dcd)

        self.simulation.reporters.append(
            StateDataReporter(stdout, self.nslog, step=True, speed=True, progress=True,
                              totalSteps=self.nsteps,
                              remainingTime=True, separator='\t\t'))

        self.simulation.reporters.append(
            StateDataReporter(out_log_file, self.nslog, step=True, kineticEnergy=True, potentialEnergy=True,
                              totalEnergy=True, temperature=True, volume=True, speed=True))

        self.simulation.step(self.nsteps)

        self.simulation.reporters.clear()  # remove all reporters so the next iteration don't trigger them

        state = self.simulation.context.getState(getPositions=True, getVelocities=True)
        with open(out_rst_file, 'w') as f:
            f.write(XmlSerializer.serialize(state)) # save the .rst file for future iterations

        positions = self.simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(self.simulation.topology, positions, open(out_pdb_file, 'w'), keepIds=True) # save the last frame

    def run_prod(self):
        print(f"Reading parameters for WESTPA2 simulation with DrudeFF at {self.configs['base_dir']}\n")
        self.setup_parameters()
        pprint.pprint(self.configs, sort_dicts=False)

        if self.configs["system_path"]: # optionally, load system from disk to save on overhead
            print("\nLoading OpenMM system\n")
            self.load_system()
        else:
            print("\nCreating OpenMM system and adding forces\n")
            self.create_system()

        print("Setting up DrudeLangevinIntegrator\n")
        self.set_integrator()

        print("Building simulation object\n")
        self.build_simulation_obj()

        print(f"Simulating for {self.nsteps} timesteps\n")
        self.simulate()

        print(f"Saved outputs to:"
              f"\nLog - {self.configs['child_prefix']}.log"
              f"\nTrajectory - {self.configs['child_prefix']}.dcd"
              f"\nRestart File - {self.configs['child_prefix']}.rst"
              f"\nLast Frame PDB - {self.configs['child_prefix']}.pdb")


def main():
    drude_prod = WEDrudeProd(
        'abl1_drude_config.json')  # config file gets moved to each iteration dir by westpa_scripts/runseg.sh
    drude_prod.run_prod()


if __name__ == "__main__":
    main()
