import numpy as np


def mean_position(
    ra: list[float], dec: list[float]
) -> tuple[float, float, float, float]:
    ra = np.radians(ra)
    dec = np.radians(dec)
    ref_ra = ra[0]
    ref_dec = dec[0]
    dx = (ra - ref_ra) * np.cos(ref_dec)
    dy = dec - ref_dec
    mean_dx = np.mean(dx)
    mean_dy = np.mean(dy)
    mean_ra = ref_ra + mean_dx / np.cos(ref_dec)
    mean_dec = ref_dec + mean_dy
    dra = ref_ra + dx / np.cos(ref_dec)
    ddec = ref_dec + dy
    std_ra = np.std(dra)
    std_dec = np.std(ddec)
    return (
        np.degrees(mean_ra),
        np.degrees(mean_dec),
        np.degrees(std_ra),
        np.degrees(std_dec),
    )
