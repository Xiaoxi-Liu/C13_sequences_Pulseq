# C13 Metabolite-Specific 3D bSSFP spiral 

Author: Xiaoxi Liu, Di Cui, 03/24,
        xiaoxi.liu@ucsf.edu, di.cui@ucsf.edu
        
The MS-3DbSSFP sequence was designed for C13 imaging with spiral readout. This sequence framework supports shifting different metabolic-specific RF pulses in the same image series. The examples include C1-pyruvate, C1-lactat[1], bicarbonate[2] and urea[3].

Different RF pulses can be set via 'set sequence parameters'

1. set metabolites number for acquisition, name and corresponding frequencies.
2. set resolution for each metabolite (support variable resolution)
3. set flip angle for each metabolite
4. Set canalization pulse flip angle list 
5. support loading customized readout gradient and built-in spiral readout generation


We have updated the corresponding functions in the toppe GitHub repository (develop branch), please check whether you are using the updated version.

Reference: 

[1] Tang, Shuyu, et al. "A metabolite‐specific 3D stack‐of‐spiral bSSFP sequence for improved lactate imaging in hyperpolarized [1‐13C] pyruvate studies on a 3T clinical scanner." Magnetic resonance in medicine 84.3 (2020): 1113-1125.

[2] Liu, Xiaoxi, et al. "A metabolite specific 3D stack-of-spirals bSSFP sequence for improved bicarbonate imaging in hyperpolarized [1-13C] Pyruvate MRI." Journal of Magnetic Resonance 353 (2023): 107518.

[3] Liu, Xiaoxi, et al. "Development of specialized magnetic resonance acquisition techniques for human hyperpolarized [13C, 15N2] urea+[1‐13C] pyruvate simultaneous perfusion and metabolic imaging." Magnetic resonance in medicine 88.3 (2022): 1039-1054.
