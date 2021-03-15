#!/usr/bin/env python3
from MDAnalysis import Universe
import mdtraj as md
import numpy as np
np.seterr(divide = 'ignore')
import argparse

def parser():
    p = argparse.ArgumentParser(description='This is Cryptic-site predictor')
    p.add_argument('-r' , '--ref'      , required = True  , help = '')
    p.add_argument('-i' , '--trj'      , required = True  , help = '')
    p.add_argument('-o' , '--out'      , required = False , help = 'output suffix'         , type=str  , default = 'cryptic')    
    p.add_argument('-a' , '--alpha'    , required = False , help = 'Upper bound of Delta F', type=float, default = 2.0 )
    p.add_argument('-b' , '--beta'     , required = False , help = 'Lower bound of sigma'  , type=float, default = 0.0)
    p.add_argument('-th', '--threshold', required = False , help = 'threshold'             , type=float, default = 0.75)
    p.add_argument('-tr', '--trj_range', required = False , help = 'traj range'            , type=int  , default = [0,-1], nargs='+')    
    args = p.parse_args()
    return args

class CrypticSitePredictor():

    def __init__(self):
        args = parser()
        self.__ref      , self.__traj     = args.ref, args.trj
        self.__alpha    , self.__beta     = args.alpha, args.beta
        self.__upper_lim_of_buried_state, self.__out_suffix = args.threshold, args.out
        self.__trj_range                  = args.trj_range
        self.__traj_data = md.load(self.__traj, top = self.__ref)

    def print_info(self):
        info = f''' 
        Reference      : {self.__ref}
        Traj           : {self.__traj}
        alpha          : {self.__alpha}
        beta           : {self.__beta}
        shreshold      : {self.__upper_lim_of_buried_state}
        Output suffix  : {self.__out_suffix}
        Trj Range      : {self.__trj_range} (Currently Inactive)
        '''
        print(info)

    def run_SASA(self):
        import pandas as pd
        PROBE_RADIUS, N_SPHERE_POINTS = 0.14, 100
        print(self.__traj_data)

        def get_residue_tag():
            """
            Todo:
                I want to do this function by mdtraj instead of MDAnalysis
            Return:
                residue_tags: List; e.g, [PHE1A, ARG2A, ...]
            """
            from MDAnalysis import Universe 
            u = Universe(self.__ref)
            residue_tags = []
            nchains = len(set(u.select_atoms('protein').segids))
            print(f'n chains = {nchains}')
            for chain, resn, resi in zip(u.residues.segids, u.residues.resnames, u.residues.resids):
            
                if chain == 'SYSTEM': self.chain = chain.replace('SYSTEM', 'A') # self.chain should have default vaule in __init__
                 
                residue_tags.append(resn + str(resi) + self.chain)
            
            return residue_tags

        keys = get_residue_tag()
        print(keys)
        
        sasa_matrix = md.shrake_rupley(self.__traj_data,
                                    probe_radius    = PROBE_RADIUS   ,
                                    n_sphere_points = N_SPHERE_POINTS,
                                    mode='residue')

        print('sasa data shape', sasa_matrix.shape)
        print(len(sasa_matrix[0,:]))
        
        ### Save a table storing sasa values
        df_sasa            = pd.DataFrame(sasa_matrix).round(4)
        df_sasa.columns    = keys    
        df_sasa.index.name = 'Frame No'
        df_sasa.to_csv(f'sasa_{self.__out_suffix}.csv')
        return df_sasa

    def to_rSASA(self, df_sasa):
        # The values were computed by mdtraj's sasa class.
        # unit, nm^2
        # Note: Relative SASAs of non-aromatic residues are NOT computed.

        # These values were caluculated for an ALA-X-ALA tripeptide.
        sasa_values_of_naked_aa = {'PHE': 2.15,
                                   'TYR': 2.30,
                                   'TRP': 2.57,
                                   'HIS': 1.90,
                                   'ARG': 2.39}

        def to_rSASA_each_series(ser_sasa):
            key = ser_sasa.name
            
            if key[0:3] == 'PHE':
                print(key)
                rsasa= ser_sasa / sasa_values_of_naked_aa['PHE']

            elif key[0:3] == 'TYR':
                print(key)
                rsasa= ser_sasa / sasa_values_of_naked_aa['TYR']

            elif key[0:3] == 'TRP':
                print(key)
                rsasa= ser_sasa / sasa_values_of_naked_aa['TRP']

            elif key[0:3] == 'HIS':
                print(key)
                rsasa= ser_sasa / sasa_values_of_naked_aa['HIS']

            elif key[0:3] == 'ARG':
                print(key)
                rsasa= ser_sasa / sasa_values_of_naked_aa['ARG']

            else:
                print(f'{key} was not converted, because non-aromatic aa are not handled by this program.')        
                rsasa= ser_sasa
 
            return rsasa
 
        df_rsasa = df_sasa.apply(to_rSASA_each_series)
        df_rsasa.to_csv(f'rsasa_{self.__out_suffix}.csv')        
        return df_rsasa

    def ratio_RbRe(self, rSASA, buried_upper_limit):
        """
        Args:
            rSASA: [dict] relative SASA like {PHE105:ASA, ..., TYR200:ASA}
            buried_upper_limit: [float] upper limit of buriedness (standardized SASA) [no unit]
        Return:
            Rbe: Ratio of exposed probability to buried one. 
                 Small value indicates buried states are more dominant than exposed.
        """
        bin_start, bin_end, binwidth = 0, 1, 0.03 # Note 'binwidth' does not affect the result, because ratio of Pe to Pb is taken and thus binwidth is offseted.
        hist = np.histogram(rSASA, bins=[i for i in np.arange(bin_start,bin_end,binwidth)])
        probs, bin_edges = hist[0], hist[1]
        half_width       = 0.5 * abs(bin_edges[0] - bin_edges[1])
        hist_xpoints     = np.array([edge+half_width for edge in bin_edges[:-1]])

        index_Pe = np.where(hist_xpoints > buried_upper_limit)[0]
        index_Pb = np.where(hist_xpoints <= buried_upper_limit)[0]
        Pe       = [probs[i] for i in index_Pe]
        xe       = hist_xpoints[hist_xpoints > buried_upper_limit]
        Pb       = [probs[j] for j in index_Pb]
        xb       = hist_xpoints[hist_xpoints <= buried_upper_limit]

        if np.sum(Pe) == 0.0:
            Rbe = np.inf

        else:
            Rbe = np.sum(Pb) / np.sum(Pe) # R_eb, the more it is, the more exposed and vice versa.

        return Rbe
    
    def run(self):
        df_sasa  = self.run_SASA()
        df_rsasa = self.to_rSASA(df_sasa)
        self.generate_cryptic_site_index(df_rsasa)

    def generate_cryptic_site_index(self, df_rsasa):
        """
            df_rsasa: dataframe (e.g., {PHE1:0.4, TYR2:0.2,...})
            self.cryptic_index: [dict] delta F and sigma are assigned to each aromatic residue
        Note: this function works only for the aromatic residues.
        """
        S_aromatic_resname = set(['PHE','TRP','TYR','HIS'])
        self.cryptic_index = {}
        for key in df_rsasa:

            if key[0:3] in S_aromatic_resname:

                # -- Calculate the ratio (Rbe) of buried to exposed states
                Rbe = self.ratio_RbRe(df_rsasa[key], self.__upper_lim_of_buried_state)
                RT = 0.59 # at 300 K
                
                if Rbe == np.inf:
                    dF = - np.inf

                else:
                    dF = - RT * np.log(Rbe)

                std_sasa = np.std(df_rsasa[key])
                
                self.cryptic_index[key] = (dF, std_sasa)

        #self.cryptic_index = sorted(self.cryptic_index.items(), key=lambda x:int(x[0][3:-1])) # sorted by residue number
        #self.cryptic_index = sorted(self.cryptic_index.items(), key=lambda x:x[1]) # sorted by sigma 

    def output_pdb_w_index(self):
        #This scales sigma. The reason for this is because PDB files accepts few significant digits/
        # sigma is usually 10^2 ~ 10^3 order, so if sigma was 0.011, then the sigma value to be written would be 0.01 in the PDB. I want to avoid this. 
        scale_factor = 100.0 

        u = Universe(self.__ref)        
        #initialize the b-factor column
        u.atoms.tempfactors = 0
        for icalpha in u.atoms.select_atoms('name CA'):
            if icalpha.resname in ['PHE','TRP','TYR','HIS']:
                 key = icalpha.resname + str(icalpha.resid) + self.chain
                 DF    = self.cryptic_index[key][0]
                 sigma = self.cryptic_index[key][1]
                 print(key, DF, sigma)
                 if np.abs(DF) < self.__alpha:
#                     print(key, DF, sigma)
                     icalpha.tempfactor = sigma * scale_factor

        u.select_atoms('protein').write(f'index_{self.__out_suffix}.pdb')

def main():
    CSP = CrypticSitePredictor()
    CSP.print_info()
    CSP.run()
    CSP.output_pdb_w_index()
    
#    CSP.df_rsasa.to_csv('rsasa.csv')

if __name__ == '__main__':
    main()
