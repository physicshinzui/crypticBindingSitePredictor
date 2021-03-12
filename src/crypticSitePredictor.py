from MDAnalysis import Universe
import mdtraj as md
import numpy as np
np.seterr(divide = 'ignore')
import argparse

def parser():
    p = argparse.ArgumentParser(description='This is Cryptic-site predictor')
    p.add_argument('-r' , '--ref'      , required = True  , help = '')
    p.add_argument('-i' , '--trj'      , required = True  , help = '')
    p.add_argument('-o' , '--out'      , required = False , help = 'output name'           , type=str  , default = 'cryptic_site_index.csv')    
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
        self.__threshold, self.__out_name = args.threshold, args.out
        self.__trj_range                  = args.trj_range
        self.__traj_data = md.load(self.__traj, top = self.__ref)

    def print_info(self):
        info = f''' 
        Reference: {self.__ref}
        Traj     : {self.__traj}
        alpha    : {self.__alpha}
        beta     : {self.__beta}
        shreshold: {self.__threshold}
        Output   : {self.__out_name}
        Trj Range: {self.__trj_range}
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
                if chain == 'SYSTEM': chain = chain.replace('SYSTEM', 'A')
                residue_tags.append(resn+str(resi)+chain)
            return residue_tags

        keys = get_residue_tag()
        print(keys)
        
        self.sasa_matrix = md.shrake_rupley(self.__traj_data,
                                    probe_radius    = PROBE_RADIUS   ,
                                    n_sphere_points = N_SPHERE_POINTS,
                                    mode='residue')

        print('sasa data shape', self.sasa_matrix.shape)
        print(len(self.sasa_matrix[0,:]))
        
        ### Save a table storing sasa values
        df = pd.DataFrame(self.sasa_matrix).round(4)
        df.columns = keys    
        df.index.name = 'Frame No'
        df.to_csv('sasa.csv')

    def to_rSASA(self):
        # The values were computed by mdtraj's sasa class.
        # unit, nm^2
        # Note: Relative SASAs of non-aromatic residues are NOT computed.

        # These values were caluculated for an ALA-X-ALA tripeptide.
        sasa_values_of_naked_aa = {'PHE': 2.15,
                                   'TYR': 2.30,
                                   'TRP': 2.57,
                                   'HIS': 1.90,
                                   'ARG': 2.39}
        key = self.sasa_matrix.name
        
        if key[0:3] == 'PHE':
            print(key)
            rSASA = self.sasa_matrix / sasa_values_of_naked_aa['PHE']

        elif key[0:3] == 'TYR':
            print(key)
            rSASA = self.sasa_matrix / sasa_values_of_naked_aa['TYR']

        elif key[0:3] == 'TRP':
            print(key)
            rSASA = self.sasa_matrix / sasa_values_of_naked_aa['TRP']

        elif key[0:3] == 'HIS':
            print(key)
            rSASA = self.sasa_matrix / sasa_values_of_naked_aa['HIS']

        elif key[0:3] == 'ARG':
            print(key)
            rSASA = self.sasa_matrix / sasa_values_of_naked_aa['ARG']

        else:
            print(f'{key} was not converted, because non-aromatic aa are not handled by this program.')        
            rSASA = self.sasa_matrix

        return rSASA


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

    def predict(self, verbose=False):
        '''
        From RSASA, it predicts high score aromatic residues
        '''
        pass


def main():
    CSP = CrypticSitePredictor()
    CSP.print_info()
    CSP.run_SASA()

if __name__ == '__main__':
    main()