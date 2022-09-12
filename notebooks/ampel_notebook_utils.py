import requests

from ampel.log.AmpelLogger import AmpelLogger
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper
from ampel.ztf.ingest.ZiDataPointShaper import ZiDataPointShaper
from ampel.content.T1Document import T1Document
from ampel.view.LightCurve import LightCurve


def api_name2alerts(name, token):
    """
    Find all alerts belonging to ZTF transient
    """

    # Setup connection
    endpoint = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/object/{}/alerts?".format(name)
    header = {"Authorization": "bearer "+token}
    response = requests.get(endpoint, headers=header )

    if not response.ok:
        print('... failed to get alert for {}'.format(name))
        return None

    # Find the candidate ID of the final alert from this event
    return response.json()


def api_name2candids(name, token):
    """
    Find all candids of alerts belonging to ZTF transient
    """

    # Setup connection
    endpoint = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/object/{}/alerts?".format(name)
    header = {"Authorization": "bearer "+token}
    response = requests.get(endpoint, headers=header )

    if not response.ok:
        print('... failed to get alert for {}'.format(name))
        return None

    # Find the candidate ID of the final alert from this event
    alerts = response.json()

    return [a['candid']for a in alerts]

def api_name2candid(name, token):
    """
    Find candid of latest alert belonging to ZTF transient
    """

    # Setup connection
    endpoint = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/object/{}/alerts?".format(name)
    header = {"Authorization": "bearer "+token}
    response = requests.get(endpoint, headers=header )

    if not response.ok:
        print('... failed to get alert for {}'.format(name))
        return None

    # Find the candidate ID of the final alert from this event
    alerts = response.json()

    return sorted(alerts, key=lambda i: (i['candidate']['jd']), reverse=True)[0]['candid']


def api_get_lightcurve(name, token, shaper=None):
    """
    Retrieve the alert history of a SN and convert to a LightCurve object.

    An alert token is needed for access. Archive token can be retrieved from
    https://ampel.zeuthen.desy.de/live/dashboard/tokens
    once the user has registered with the github AmpelProject.

    :param name: str, ZTF name of transient.
    :param token: str, AMPEL archive token:
    :param shaper: optional, AMPEL data point shaper (otherwise using default)

    return LightCurve object

    """

    # Setup connection
    endpoint = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/object/{}/photopoints".format(name)
    header = {"Authorization": "bearer "+token}

    response = requests.get(endpoint, headers=header )

    if not response.ok:
        print('... failed to get alert')
        return None

    # Convert
    alert = response.json()
    if alert is None:
        print(' ... no alert content')
        return None

    pps = [alert['candidate']]
    pps.extend( [prv_cand for prv_cand in alert['prv_candidates'] ] )

    if shaper is None:
        tmplog = AmpelLogger.get_logger()
        shaper = ZiDataPointShaper(logger=tmplog)

    stockId = ZTFIdMapper.to_ampel_id(name)
    dps = shaper.process( pps, stockId)
    t1d = T1Document(stock=stockId, link=0)
    return LightCurve.build(t1d, dps)


def api_get_lightcurves(name, token, shaper=None):
    """

    Retrieve alert history of a SN and convert to a list of LightCurve objects.

    An alert token is needed for access. Archive token can be retrieved from
    https://ampel.zeuthen.desy.de/live/dashboard/tokens
    once the user has registered with the github AmpelProject.

    :param name: str, ZTF name of transient.
    :param token: str, AMPEL archive token:
    :param shaper: optional, AMPEL data point shaper (otherwise using default)

    return List[LightCurve object]

    """

    if shaper is None:
        shaper = ZiDataPointShaper(logger=AmpelLogger.get_logger())

    # Setup connection
    endpoint = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/object/{}/alerts?with_history=true&with_cutouts=false".format(name)
    header = {"Authorization": "bearer "+token}

    response = requests.get(endpoint, headers=header )

    if not response.ok:
        print('... failed to get alert')
        return None

    # Convert
    alerts = response.json()
    if alerts is None:
        print(' ... no alert content')
        return None

    lcs = []

    for alert in alerts:
        pps = [alert['candidate']]
        pps.extend( [prv_cand for prv_cand in alert['prv_candidates'] ] )

        stockId = ZTFIdMapper.to_ampel_id(name)
        dps = shaper.process( pps, stockId)
        t1d = T1Document(stock=stockId, link=0)
        lcs.append( LightCurve.build(t1d, dps) )

    return lcs


DecentFilterSingle = {
    'min_ndet': 1,
    'min_tspan': -99,
    'max_tspan': 10,
    'min_archive_tspan': -99,
    'max_archive_tspan': 10,
    'min_rb': 0.2,
    'min_drb': 0.995,   # Default for PossisHQ
    'max_fwhm': 5.5,
    'min_gal_lat': 0,
    'ps1_sgveto_rad': 1,
    'ps1_sgveto_th': 0.8,
    'ps1_confusion_rad': 3,
    'ps1_confusion_sg_tol': 0.1,
    'gaia_rs': 20,
    'gaia_pm_signif': 3,
    'gaia_plx_signif': 3,
    'gaia_veto_gmag_min': 9,
    'gaia_veto_gmag_max': 20,
    'gaia_excessnoise_sig_max': 999,
    'max_elong': 1.4,
    'max_magdiff': 1,
    'max_nbad': 2,
    'min_sso_dist': 20,
}
T2CatalogMatchMultipleFirstPPS = {
    'catalogs' : {
        'SDSS_spec' : {
            'use' : 'extcats',
            'rs_arcsec' : 10.0,
            'keys_to_append' : ['z', 'bptclass', 'subclass'],
            'all': False,
        },
        'NEDz' : {
            'use' : 'catsHTM',
            'rs_arcsec' : 10.0,
            'keys_to_append' : ['ObjType', 'Velocity', 'z'],
        },
        'GLADEv23' : {
            'use' : 'extcats',
            'rs_arcsec' : 10,
            'keys_to_append' : ['z', 'dist', 'dist_err', 'flag1', 'flag2', 'flag3'],
        },
        'LSPhotoZZou' : {
            'use' : 'extcats',
            'rs_arcsec' : 10.0,
            'keys_to_append' : ['photoz','ra','dec','e_photoz','specz','_6','logMassBest','logMassInf','logMassSup'],
            'pre_filter' : None,
            'post_filter' : None,
            'all': False,
        },
        'wiseScosPhotoz' : {
            'use' : 'extcats',
            'rs_arcsec' : 10.0,
            'keys_to_append' : ['zPhoto_Corr','ra','dec','wiseID','w1mCorr','w2mCorr'],
            'pre_filter' : None,
            'post_filter' : None,
        },
        'twoMPZ' : {
            'use' : 'extcats',
            'rs_arcsec' : 10.0,
            'keys_to_append' : ['zPhoto','ra','dec','zSpec'],
            'pre_filter' : None,
            'post_filter' : None,
        },

    }
}

T2DigestRedshiftsMultipleFirstPPSPhotoZ = {
    "max_redshift_category" : 7,
    "t2_dependency": [
        {
            "unit": "T2CatalogMatch",
            "config": T2CatalogMatchMultipleFirstPPS,
            "link_override": {
                'filter': 'PPSFilter', 'sort': 'jd', "select": "first"
                }
        },
    ]
}

T2PropagateStockInfoHealpix = {
    'prop_paths': {'explosion_time':['journal','healpix','time'] }
}

T2RunPossisStockTriggerPhotoZ = possis_config = {
#    "possis_dir": "/home/jnordin/data/models/kilonova_models",
#    "explosion_time_jd": Time(trigger_time, scale="utc").jd,
    "explosion_time_jd": "StockTriggerTime",
    "redshift_kind" : 'T2DigestRedshifts',
    "max_ampelz_group" : 7,      # For this purpose we use any available redshift
    "t2_dependency": [
        {
            "unit": "T2DigestRedshifts",
            "config": T2DigestRedshiftsMultipleFirstPPSPhotoZ,
        },
        {
            "unit": "T2PropagateStockInfo",
            "config": T2PropagateStockInfoHealpix,
        },
    ],
    "plot_props" : {
        "tags":["POSSIS"],
        "file_name":{
            "format_str": "%s.svg",
            "arg_keys": ["stock"]
        },
        "title":{
            "format_str": "%s",
            "arg_keys": ["stock"]
        },
        "fig_text":{
            "format_str": "%s %s \nz-source %s \nchisq %.2f ndof %s",
            "arg_keys": ["stock", "redshift_kind", "chisq", "ndof"]
        },
        "width":10,
        "height":6,
        "compress":2,
        "id_mapper":"ZTFIdMapper",
        "disk_save":"/home/jnordin/tmp/possis",
    }

}



def standard_configs(unit: str, mode: str = None):
    """
    Return a standard configuration for use with a particular unit.
    The optional mode is used e.g. for chained units.

    Current options: (unit, mode)
    - DecentFilter, Single
    - T2CatalogMatch, MultipleFirstPPS
    - T2DigestRedshifts, MultipleFirstPPSPhotoZ
    - T2RunPossis, StockTriggerPhotoZ
    - T2PropagateStockInfo, Healpix

    return dict
    """

    if unit=='DecentFilter':
        if mode=='Single':
            return DecentFilterSingle
    elif unit=='T2CatalogMatch':
        if mode=='MultipleFirstPPS':
            return T2CatalogMatchMultipleFirstPPS
    elif unit=='T2DigestRedshifts':
        if mode=='MultipleFirstPPSPhotoZ':
            return T2DigestRedshiftsMultipleFirstPPSPhotoZ
    elif unit=='T2RunPossis':
        if mode=='StockTriggerPhotoZ':
            return T2RunPossisStockTriggerPhotoZ
    elif unit=='T2PropagateStockInfo':
        if mode=='Healpix':
            return T2PropagateStockInfoHealpix
