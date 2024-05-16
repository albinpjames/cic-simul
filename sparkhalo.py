from scipy.stats import binned_statistic_dd 

class SimulData(object):
    def __init__(self, 
                 params,
                 boxsize: int = None, # The box size of the simulation
                 mass: float = None # The mass of the particles in the simulation
    ) -> None:
        self.boxsize = boxsize if boxsize is not None else params["boxsize"]  
        self.mass = mass if mass is not None else params["mass"]

    def cic(self,nw_boxsize):
        cutside = int(self.boxsize/nw_boxsize)

        x_bins_dd = pos_start[0] + np.linspace(0,nw_boxsize*cutside,cutside+1)
        y_bins_dd = pos_start[1] + np.linspace(0,nw_boxsize*cutside,cutside+1)
        z_bins_dd = pos_start[2] + np.linspace(0,nw_boxsize*cutside,cutside+1)
        try:
            boxdata = binned_statistic_dd(data
                ,values = None, statistic = 'count', 
                bins =[x_bins_dd, y_bins_dd, z_bins_dd]).statistic
        except Exception as __e:
            print(x_bins_dd)
            print(y_bins_dd)
            print(z_bins_dd)
            raise __e

        boxdata = boxdata.ravel()

        return boxdata

class abacussummit(SimulData):
    def __init__(self,
                 params,
                 boxsize: int = None, # The box size of the simulation
                 mass: float = None # The mass of the particles in the simulation
                 type : str = None,
                 cosmo : str = None,
                 intcont : str = None ) -> None:
        super().__init__(type,params,boxsize,mass)
        self.type = type if type is not None else params["type"]
        self.cosmo = cosmo if cosmo is not None else params["cosmo"]
        self.intcont = intcont if intcont is not None else params["intcont"]
        self.dataloc = dataloc if dataloc is not None else params["dataloc"]
    
    def readdata(self, redshift: float = None, clms: list = None):

    

    

if __init__ == '__main__':
    
    params = {
        "dataloc" : "C:/"
        "name": "abacussummit",  # The name of the simulation
        "type": "base",  # The type of simulation (base, small for abacussumit)
        "cosmo": "c000",  # The cosmology used in the simulation
        "intcont": "ph000",  # The initial condition used for the simulation
        "boxsize": 2000,  # The box size of the simulation
        "mass": 2.109081520453063 * 10**9,  # The mass of the particles in the simulation
        }
    
    simulation = abacussummit(params)
    """ Redshifts & boxsizes to be computed """
    # redshifts = ["3.000","2.500","2.000"]
    redshifts = ["3.000"]
    # nw_boxsizes = [30]
    nw_boxsizes = [50,15,10,5]
   
    for redshift in redshifts:
        print(f"redshift being computed: {redshift}")
        clms = ["N","SO_central_particle"]
        simulation.readdata(redshift,clms)

