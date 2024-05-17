
import os
from pathlib import Path
import numpy as np
import math
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from scipy.stats import binned_statistic_dd 

class SimulData(object):
    def __init__(self, 
                 params,
                 boxsize: int = None, # The box size of the simulation
                 mass: float = None # The mass of the particles in the simulation
    ) -> None:
        self.boxsize = boxsize if boxsize is not None else params["boxsize"]  
        self.mass = mass if mass is not None else params["mass"]

    def cic(self,nw_boxsize,massrange):
        
        cutside = int(self.boxsize/nw_boxsize)
        totalboxes = int(cutside**3) 
        self.cicdata = np.zeros(totalboxes)

        x_bins_dd = np.linspace(0,nw_boxsize*cutside,cutside+1)
        y_bins_dd = np.linspace(0,nw_boxsize*cutside,cutside+1)
        z_bins_dd = np.linspace(0,nw_boxsize*cutside,cutside+1)

        data = self.cat[["xpos","ypos","zpos"]][(massrange[0] <= self.cat["halomass"]) & (self.cat["halomas"] < massrange[1])]
        data = np.array([data['xpos'], data['ypos'], data['zpos']]).T
        boxdata = binned_statistic_dd(data
                ,values = None, statistic = 'count', 
                bins =[x_bins_dd, y_bins_dd, z_bins_dd]).statistic

        self.cicdata = boxdata.ravel()

    def getcicdata(self):
        return self.cicdata

class abacussummit(SimulData):
    def __init__(self,
                 params,
                 boxsize: int = None, # The box size of the simulation
                 mass: float = None, # The mass of the particles in the simulation
                 type : str = None,
                 cosmo : str = None,
                 intcont : str = None,
                 dataloc: str = None  ) -> None:
        super().__init__(params,boxsize,mass)
        self.type = type if type is not None else params["type"]
        self.cosmo = cosmo if cosmo is not None else params["cosmo"]
        self.intcont = intcont if intcont is not None else params["intcont"]
        self.dataloc = dataloc if dataloc is not None else params["dataloc"]
    
    def readdata(self, redshift: float = None, clms: list = None):
        print("Reading the data")

        # Location of the data
        file = os.path.join(
            self.dataloc,
            "AbacusSummit_Public_Data_Access/AbacusSummit_"
            + self.type
            + "_"
            + self.cosmo
            + "_"
            + self.intcont,
            "halos/z" + redshift,
            "halo_info/",
        )

        self.cat = CompaSOHaloCatalog(file, cleaned=False)  # Reads the data
        self.cat = self.cat.halos[clms]  # Reads the given column and saves it to an array

        if "N" in clms:
            print("converting to mass")
            self.cat["halomass"] = self.cat["N"] * self.mass
        
        if "SO_central_particle" in clms:
            self.cat["SO_central_particle"] += self.boxsize / 2

            self.cat.add_columns(
                [
                    self.cat["SO_central_particle"][:, 0],
                    self.cat["SO_central_particle"][:, 1],
                    self.cat["SO_central_particle"][:, 2],
                ],
                names=["xpos", "ypos", "zpos"],
            )

            del self.cat["SO_central_particle"]

    

if __name__ == '__main__':
    
    params = {
        "dataloc" : "/mnt/data/DATA/Simulations",
        "name": "abacussummit",  # The name of the simulation
        "type": "base",  # The type of simulation (base, small for abacussumit)
        "cosmo": "c000",  # The cosmology used in the simulation
        "intcont": "ph000",  # The initial condition used for the simulation
        "boxsize": 2000,  # The box size of the simulation
        "mass": 2.109081520453063 * 10**9,  # The mass of the particles in the simulation
        }
    
    simul = abacussummit(params)

    """ Redshifts & boxsizes to be computed """
    redshifts = ["3.000"]
    nw_boxsizes = [50,15,10,5]
    massrange = [10**5,10**15]
   
    for redshift in redshifts:
        print(f"redshift being computed: {redshift}")
        clms = ["N","SO_central_particle"]
        simul.readdata(redshift,clms)
        print(simul.cat)

        # for nw_boxsize in nw_boxsizes: 
        #     print(f"Boxsize being computed: {nw_boxsize}")
        #     simul.cic(nw_boxsize,massrange)
        #     print(simul.getcicdata())








