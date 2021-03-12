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

    def info(self):
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

    def rSASA_calculation(self):
        traj = md.load(__traj, top = __ref)
        pass

    def run(self, verbose=False):
        pass

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


        pass

def main():
    CSP = CrypticSitePredictor()
    CSP.info()
if __name__ == '__main__':
    main()