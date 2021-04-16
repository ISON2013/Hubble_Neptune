def get_ephemerides(hdr=None, return_eph=True, loc=None, epoch=None, target=None,
                    cache=False, **kwargs):
    """
    Loads ephemeris data from JPL horizons for provided FITS header object.

    Telescope location, target and observation datetime can be found from
    header object. See astroquery documentation for more details.

    Parameters
    ----------
    hdr
        FITS header object for observation requiring ephemeris data. Set to
        None to use data from loc, epoch and target instead.

    return_eph : bool
        Toggle returning full ephemeris object or returning default mapping
        parameters.

    loc : dict
        Location of observer.

    epoch : dict
        Observation timings.

    target : str
        Observation target.

    cache : bool
        Toggle saving/using locally cached results.

    Returns
    -------
    ephemerides
    """
    # Get info from FITS header
    if target is None:
        target = hdr['OBJECT'].casefold()
    target_dict = {
            'io': 501,
            'europa': 502,
            'jupiter': 599,
            'uranus' : 799,
            'neptune': 899,
            }
    if target in target_dict:
        target = target_dict[target] # Correct names to unambiguous ID if necessary
    if epoch is None:
        epoch = hdr['EXPSTART']
    if loc is None:
        loc = {'lon': hdr['ESO TEL GEOLON'],
               'lat': hdr['ESO TEL GEOLAT'],
               'elevation': hdr['ESO TEL GEOELEV']/1e3  # m -> km
               }
    SpaceTelescope_dict = {
            'hst': '@hst',
            'jwst': '@-170',
            }
    if loc in SpaceTelescope_dict:
        loc = SpaceTelescope_dict[loc]
    # Load ephemeris data
    conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'
    obj = Horizons(id=target, location=loc, id_type='majorbody', epochs=epoch)
    # Astroquery automatically caches query result locally, so should still work if no internet
    # connection is available and exact query has been performed before on the same machine.

    if cache:
        # Use cached results from previous runs etc. to avout bottleneck
        if cache is True:
            cache = 'rw'
        cache_path = tools.path.code('cache', 'ephemerides.pickle.gz')
        try:
            with gzip.open(cache_path, 'rb') as f:
                cache_db = pickle.load(f)
        except FileNotFoundError:
            cache_db = {}

        key = str(obj)
        if key in cache_db and 'r' in cache:
            eph = cache_db[key]
        else:
            eph = obj.ephemerides(**kwargs)

        if 'w' in cache and key not in cache_db:
            cache_db[key] = eph
            tools.file.check_path(cache_path)
            with gzip.open(cache_path, 'wb') as f:
                pickle.dump(cache_db, f)
    else:
        eph = obj.ephemerides(**kwargs)

    if return_eph:
        return eph
    # Return data relevant for mapping
    obs_long = -float(eph['PDObsLon'])
    obs_lat = float(eph['PDObsLat'])
    sun_long = -float(eph['PDSunLon'])
    sun_lat = float(eph['PDSunLat'])
    np_ang = float(eph['NPole_ang'])
    return obs_long, obs_lat, sun_long, sun_lat, np_ang
