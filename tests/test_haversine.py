import numpy as np
from acrg.grid.haversine import distance

def test_acrg_haversine():
    '''
    Test if the shape of the regriddded 2D array is correct
    '''
    test_distance = distance([18.,25.],[19.,26.])
    test_distance2 = distance([18.,29.],[19.,33.])
    test_distance3 = distance([-28.,29.],[-30.,33.])
    test_distance4 = distance([-28.,29.],[28.,33.])
    assert np.round(test_distance,6) == 153.242785
    assert np.round(test_distance2,6) == 436.190463
    assert np.round(test_distance3,6) == 448.045143
    assert np.round(test_distance4,6) == 6241.498555


