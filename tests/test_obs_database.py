from acrg_obs.utils import obs_database
from acrg_config.paths import paths


def test_obs_database():
    """Create obs_database
    """
    
    obs_database()

    # This isn't a very good test...
    assert (paths.obs / "obs.db").exists()

