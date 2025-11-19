"""
Stripper to CSS lattice definition.
Reference beam: 64Zn28.5+ at 19.8500 MeV/u
"""

import numpy as np
import pyJuTrack as jt

def create_stripper2CSS_lattice():
    """
    Create pyJuTrack lattice for stripper to CSS section.
    
    Returns:
        lattice: jt.Lattice object (apertures defined in element RApertures)
    """
    # Define all unique elements
    FS1_STRL_STRIP_D2237 = jt.MARKER("FS1_STRL_STRIP_D2237")
    trg_drift_D2248 = jt.DRIFT("trg_drift_D2248", 1.112338, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_BPM_D2248 = jt.MARKER("FS1_CSS_BPM_D2248")
    trg_drift_D2254 = jt.DRIFT("trg_drift_D2254", 0.457162, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QV_D2254 = jt.KQUAD("FS1_CSS_QV_D2254", 0.261, k1=3.8357202079, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2257 = jt.DRIFT("trg_drift_D2257", 0.187006, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_DCH_D2257 = jt.CORRECTOR("FS1_CSS_DCH_D2257", 0.0, xkick=0.0020071345)
    FS1_CSS_DCV_D2257 = jt.CORRECTOR("FS1_CSS_DCV_D2257", 0.0, ykick=-0.00011507654399999999)
    trg_drift_D2260 = jt.DRIFT("trg_drift_D2260", 0.201994, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QH_D2260 = jt.KQUAD("FS1_CSS_QH_D2260", 0.261, k1=-4.8150216628, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2264 = jt.DRIFT("trg_drift_D2264", 0.310501, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_BCM_D2264 = jt.DRIFT("FS1_CSS_BCM_D2264", 0.0, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2272 = jt.DRIFT("trg_drift_D2272", 0.628499, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QV_D2272 = jt.KQUAD("FS1_CSS_QV_D2272", 0.261, k1=2.6654394091, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2276 = jt.DRIFT("trg_drift_D2276", 0.233156, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_DCH_D2276 = jt.CORRECTOR("FS1_CSS_DCH_D2276", 0.0, xkick=-0.0005800467999999999)
    FS1_CSS_DCV_D2276 = jt.CORRECTOR("FS1_CSS_DCV_D2276", 0.0, ykick=-0.000132846712)
    trg_drift_D2278 = jt.DRIFT("trg_drift_D2278", 0.239394, RApertures=[0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    FS1_CSS_BPM_D2278 = jt.MARKER("FS1_CSS_BPM_D2278")
    trg_drift_D2280 = jt.DRIFT("trg_drift_D2280", 0.0543, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QH_D2280 = jt.KQUAD("FS1_CSS_QH_D2280", 0.261, k1=-0.3437955768, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2290a = jt.DRIFT("trg_drift_D2290a", 0.574094, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2290b = jt.DRIFT("trg_drift_D2290b", 0.006, RApertures=[0.015, 0.015, 0.015, 0.015, 0.015, 0.015])
    trg_drift_D2290c = jt.DRIFT("trg_drift_D2290c", 0.18, RApertures=[0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    FS1_CSS_DH_D2290_0 = jt.DRIFT("FS1_CSS_DH_D2290_0", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2290_1 = jt.DRIFT("FS1_CSS_DH_D2290_1", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2290_2 = jt.DRIFT("FS1_CSS_DH_D2290_2", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2290_3 = jt.DRIFT("FS1_CSS_DH_D2290_3", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2290_4 = jt.DRIFT("FS1_CSS_DH_D2290_4", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    trg_drift_D2296 = jt.DRIFT("trg_drift_D2296", 0.25, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_DH_D2296_0 = jt.DRIFT("FS1_CSS_DH_D2296_0", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2296_1 = jt.DRIFT("FS1_CSS_DH_D2296_1", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2296_2 = jt.DRIFT("FS1_CSS_DH_D2296_2", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2296_3 = jt.DRIFT("FS1_CSS_DH_D2296_3", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2296_4 = jt.DRIFT("FS1_CSS_DH_D2296_4", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    trg_drift_D2302 = jt.DRIFT("trg_drift_D2302", 0.35, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_DH_D2302_0 = jt.DRIFT("FS1_CSS_DH_D2302_0", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2302_1 = jt.DRIFT("FS1_CSS_DH_D2302_1", 0.06, RApertures=[0.15, 0.15,0.15,0.15,0.15,0.15])
    FS1_CSS_DH_D2302_2 = jt.DRIFT("FS1_CSS_DH_D2302_2", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2302_3 = jt.DRIFT("FS1_CSS_DH_D2302_3", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2302_4 = jt.DRIFT("FS1_CSS_DH_D2302_4", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    trg_drift_D2308 = jt.DRIFT("trg_drift_D2308", 0.25, RApertures=[0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    FS1_CSS_DH_D2308_0 = jt.DRIFT("FS1_CSS_DH_D2308_0", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2308_1 = jt.DRIFT("FS1_CSS_DH_D2308_1", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2308_2 = jt.DRIFT("FS1_CSS_DH_D2308_2", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2308_3 = jt.DRIFT("FS1_CSS_DH_D2308_3", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    FS1_CSS_DH_D2308_4 = jt.DRIFT("FS1_CSS_DH_D2308_4", 0.06, RApertures=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    trg_drift_D2313 = jt.DRIFT("trg_drift_D2313", 0.371042, RApertures=[0.02, 0.02,0.02,0.02,0.02,0.02])
    FS1_CSS_BPM_D2313 = jt.MARKER("FS1_CSS_BPM_D2313")
    trg_drift_D2325c = jt.DRIFT("trg_drift_D2325c", 0.194, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2325b = jt.DRIFT("trg_drift_D2325b", 0.006, RApertures=[0.015, 0.015,0.015,0.015,0.015,0.015])
    trg_drift_D2325d = jt.DRIFT("trg_drift_D2325d", 0.124, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_MGB01_CAV_D2325 = jt.DRIFT("FS1_MGB01_CAV_D2325", 1.14, RApertures=[0.018, 0.018,0.018,0.018,0.018,0.018])
    trg_drift_D2351a = jt.DRIFT("trg_drift_D2351a", 1.336322, RApertures=[0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    trg_drift_D2351b = jt.DRIFT("trg_drift_D2351b", 0.006, RApertures=[0.015, 0.015, 0.015, 0.015, 0.015, 0.015])
    trg_drift_D2351c = jt.DRIFT("trg_drift_D2351c", 0.6936, RApertures=[0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    FS1_CSS_DCH_D2351 = jt.CORRECTOR("FS1_CSS_DCH_D2351", 0.0, xkick=0.00019520926999999999)
    FS1_CSS_DCV_D2351 = jt.CORRECTOR("FS1_CSS_DCV_D2351", 0.0, ykick=-0.00037992356799999995)
    trg_drift_D2353 = jt.DRIFT("trg_drift_D2353", 0.206701, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_BCM_D2353 = jt.DRIFT("FS1_CSS_BCM_D2353", 0.0, RApertures=[0.015, 0.015, 0.015, 0.015, 0.015, 0.015])
    trg_drift_D2356 = jt.DRIFT("trg_drift_D2356", 0.137799, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QH_D2356 = jt.KQUAD("FS1_CSS_QH_D2356", 0.261, k1=0.0110558368, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2362a = jt.DRIFT("trg_drift_D2362a", 0.2312, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2362b = jt.DRIFT("trg_drift_D2362b", 0.006, RApertures=[0.015, 0.015,0.015,0.015,0.015,0.015])
    trg_drift_D2362c = jt.DRIFT("trg_drift_D2362c", 0.1268, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QH_D2362 = jt.KQUAD("FS1_CSS_QH_D2362", 0.261, k1=-0.3347871172, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2367 = jt.DRIFT("trg_drift_D2367", 0.3945, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_DCH_D2367 = jt.CORRECTOR("FS1_CSS_DCH_D2367", 0.0)
    FS1_CSS_DCV_D2367 = jt.CORRECTOR("FS1_CSS_DCV_D2367", 0.0)
    trg_drift_D2369 = jt.DRIFT("trg_drift_D2369", 0.1685, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_BPM_D2369 = jt.MARKER("FS1_CSS_BPM_D2369")
    trg_drift_D2372 = jt.DRIFT("trg_drift_D2372", 0.18815, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QV_D2372 = jt.KQUAD("FS1_CSS_QV_D2372", 0.261, k1=4.1021249271, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2377 = jt.DRIFT("trg_drift_D2377", 0.239, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_QH_D2377 = jt.KQUAD("FS1_CSS_QH_D2377", 0.261, k1=-3.5350014470, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    trg_drift_D2381 = jt.DRIFT("trg_drift_D2381", 0.2395, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_DCH_D2381 = jt.CORRECTOR("FS1_CSS_DCH_D2381", 0.0, xkick=-0.00024034590999999999)
    FS1_CSS_DCV_D2381 = jt.CORRECTOR("FS1_CSS_DCV_D2381", 0.0, ykick=0.000217389376)
    trg_drift_D2383 = jt.DRIFT("trg_drift_D2383", 0.26299, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_BPM_D2383 = jt.MARKER("FS1_CSS_BPM_D2383")
    trg_drift_D2385 = jt.DRIFT("trg_drift_D2385", 0.145282, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_CSS_PM_D2385 = jt.MARKER("FS1_CSS_PM_D2385")
    trg_drift_D2394 = jt.DRIFT("trg_drift_D2394", 0.381604, RApertures=[0.025, 0.025,0.025,0.025,0.025,0.025])
    FS1_BBS_DH_D2394_0 = jt.SBEND("FS1_BBS_DH_D2394_0", 0.104065, angle=-0.0785398163, e1=-0.1221730476, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_1 = jt.SBEND("FS1_BBS_DH_D2394_1", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_2 = jt.SBEND("FS1_BBS_DH_D2394_2", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_3 = jt.SBEND("FS1_BBS_DH_D2394_3", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_4 = jt.SBEND("FS1_BBS_DH_D2394_4", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_5 = jt.SBEND("FS1_BBS_DH_D2394_5", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_6 = jt.SBEND("FS1_BBS_DH_D2394_6", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_7 = jt.SBEND("FS1_BBS_DH_D2394_7", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_8 = jt.SBEND("FS1_BBS_DH_D2394_8", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=0.0000000000, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    FS1_BBS_DH_D2394_9 = jt.SBEND("FS1_BBS_DH_D2394_9", 0.104065, angle=-0.0785398163, e1=0.0000000000, e2=-0.1221730476, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    trg_drift_D2412a = jt.DRIFT("trg_drift_D2412a", 0.03, RApertures=[0.04, 0.04,0.04,0.04,0.04,0.04])
    trg_drift_D2412b = jt.DRIFT("trg_drift_D2412b", 0.59, RApertures=[0.06, 0.06,0.06,0.06,0.06,0.06])
    trg_drift_D2412c = jt.DRIFT("trg_drift_D2412c", 0.025, RApertures=[0.026, 0.026,0.026,0.026,0.026,0.026])
    trg_drift_D2412d = jt.DRIFT("trg_drift_D2412d", 0.0001, RApertures=[0.07, 0.07,0.07,0.07,0.07,0.07])

    # Combine elements into lattice sequence
    elements = [
        FS1_STRL_STRIP_D2237,
        trg_drift_D2248, FS1_CSS_BPM_D2248, trg_drift_D2254, FS1_CSS_QV_D2254, trg_drift_D2257,
        FS1_CSS_DCH_D2257, FS1_CSS_DCV_D2257, trg_drift_D2260, FS1_CSS_QH_D2260, trg_drift_D2264,
        FS1_CSS_BCM_D2264, trg_drift_D2272, FS1_CSS_QV_D2272, trg_drift_D2276, FS1_CSS_DCH_D2276,
        FS1_CSS_DCV_D2276, trg_drift_D2278, FS1_CSS_BPM_D2278, trg_drift_D2280, FS1_CSS_QH_D2280,
        trg_drift_D2290a,trg_drift_D2290b,trg_drift_D2290c, FS1_CSS_DH_D2290_0, FS1_CSS_DH_D2290_1, FS1_CSS_DH_D2290_2, FS1_CSS_DH_D2290_3,
        FS1_CSS_DH_D2290_4, trg_drift_D2296, FS1_CSS_DH_D2296_0, FS1_CSS_DH_D2296_1, FS1_CSS_DH_D2296_2,
        FS1_CSS_DH_D2296_3, FS1_CSS_DH_D2296_4, trg_drift_D2302, FS1_CSS_DH_D2302_0, FS1_CSS_DH_D2302_1,
        FS1_CSS_DH_D2302_2, FS1_CSS_DH_D2302_3, FS1_CSS_DH_D2302_4, trg_drift_D2308, FS1_CSS_DH_D2308_0,
        FS1_CSS_DH_D2308_1, FS1_CSS_DH_D2308_2, FS1_CSS_DH_D2308_3, FS1_CSS_DH_D2308_4, trg_drift_D2313,
        FS1_CSS_BPM_D2313, trg_drift_D2325c,trg_drift_D2325b,trg_drift_D2325c,trg_drift_D2325d, FS1_MGB01_CAV_D2325, trg_drift_D2351a,trg_drift_D2351b,trg_drift_D2351c, FS1_CSS_DCH_D2351,
        FS1_CSS_DCV_D2351, trg_drift_D2353, FS1_CSS_BCM_D2353, trg_drift_D2356, FS1_CSS_QH_D2356,
        trg_drift_D2362a,trg_drift_D2362b,trg_drift_D2362c, FS1_CSS_QH_D2362, trg_drift_D2367, FS1_CSS_DCH_D2367, FS1_CSS_DCV_D2367,
        trg_drift_D2369, FS1_CSS_BPM_D2369, trg_drift_D2372, FS1_CSS_QV_D2372, trg_drift_D2377,
        FS1_CSS_QH_D2377, trg_drift_D2381, FS1_CSS_DCH_D2381, FS1_CSS_DCV_D2381, trg_drift_D2383,
        FS1_CSS_BPM_D2383, trg_drift_D2385, FS1_CSS_PM_D2385, trg_drift_D2394,
        FS1_BBS_DH_D2394_0,FS1_BBS_DH_D2394_1,FS1_BBS_DH_D2394_2,FS1_BBS_DH_D2394_3,FS1_BBS_DH_D2394_4,FS1_BBS_DH_D2394_5,
        FS1_BBS_DH_D2394_6,FS1_BBS_DH_D2394_7,FS1_BBS_DH_D2394_8,FS1_BBS_DH_D2394_9,
        trg_drift_D2412a,trg_drift_D2412b,trg_drift_D2412c,trg_drift_D2412d
    ]

    # Create the lattice
    lattice = jt.Lattice(elements)

    return lattice