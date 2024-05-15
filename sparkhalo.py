

class SimulData:
    def __init__(self, 
                 data,
                 boxsize: int = None, # The box size of the simulation
                 mass: float = None # The mass of the particles in the simulation
    ) -> None:
        self.boxsize = boxsize if boxsize is not None else data["boxsize"]  
        self.mass = mass if mass is not None else data["mass"]

    def 

class abacussummit(SimulData):
    def __init__(self,
                 data,
                 boxsize: int = None, # The box size of the simulation
                 mass: float = None # The mass of the particles in the simulation
                 type : str = None,
                 cosmo : str = None,
                 intcont : str = None ) -> None:
        super().__init__(type,data,boxsize,mass)
        self.type = type if type is not None else data["type"]
        self.cosmo = cosmo if cosmo is not None else data["cosmo"]
        self.intcont = intcont if intcont is not None else data["intcont"]
