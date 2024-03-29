# C13 2D GRE spiral (single-echo or multi-echo)

Author: Xiaoxi Liu, Di Cui, 03/17,
        xiaoxi.liu@ucsf.edu, di.cui@ucsf.edu
        
The 2D GRE sequence was designed for C13 imaging with spiral readout, supporting single-echo and multi-echo mode.

Different modes can be set via 'set sequence parameters'

1. set metabolites number for acquisition and corresponding frequencies.
2. set resolution for each metabolite (support variable resolution)
3. set echo numbers (support arbitrary >=1 echo for each metabolite)
4. decide echo spacing (only valid when there is a >1 echo number)
5. set TR for each metabolite (considering echo numbers)
6. set flip angle for each metabolite (support arbitrary flip angle for each metabolite)
7. only support loading customized RF pulse
8. support loading customized readout gradient and built-in spiral readout generation

Current readout trajectory is spiral. If you want to use EPI/cartesian or other trajectories, please refer to examples in pulseq GitHub or contact us (xiaoxi.liu@ucsf.edu, di.cui@ucsf.edu) for further support.

We have updated the corresponding functions in the toppe GitHub repository (develop branch), please check whether you are using the updated version.

Reference: Liu, X., Cui, D., Xu, D., Bok, R., Wang, Z.J., Vigneron, D.B., Larson, P.E. and Gordon, J.W., 2024. Dynamic T2* relaxometry of hyperpolarized [1‐13C] pyruvate MRI in the human brain and kidneys. Magnetic Resonance in Medicine, 91(3), pp.1030-1042.
