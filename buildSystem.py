def buildSystem():
    print('Building system...')

    psf = CharmmPsfFile('./inputFiles/fep-40-i01-1-wb-i.psf')
    pdb = PDBFile('./inputFiles/fep-45-i01-1-wb-i.pdb')
    topology = psf.topology
    positions_init = pdb.positions
    params = CharmmParameterSet('./inputFiles/all_top.rtf', './inputFiles/parameters.prm')
    nonbondedMethod = PME
    nonbondedCutoff = 12*angstroms
    switchDistance=10*angstroms

    ewaldErrorTolerance = 0.0005
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 0.000001

    # Integration Options

    dt = 0.002*picoseconds
    temperature = 276*kelvin
    friction = 1.0/picosecond

    x_PBC_vector_length = 6.05
    y_PBC_vector_length = 5.88
    z_PBC_vector_length = 5.25
  
    psf.setBox(x_PBC_vector_length, y_PBC_vector_length, z_PBC_vector_length)
    platform = Platform.getPlatformByName('CUDA')
    platformProperties = {'DeviceIndex': '0', 'Precision': 'mixed'}

    return psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, switchDistance=switchDistance)
