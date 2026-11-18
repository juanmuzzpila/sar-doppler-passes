"""Geodesy helpers: WGS84 ellipsoid, ECEF/geodetic conversions, TEME -> ECEF."""

from __future__ import annotations

import numpy as np

# WGS84
A = 6378137.0  # semi-major axis [m]
F = 1.0 / 298.257223563  # flattening
E2 = F * (2.0 - F)  # first eccentricity squared
OMEGA_E = 7.292115146706979e-5  # Earth rotation rate [rad/s]
C = 299_792_458.0  # speed of light [m/s]


def geodetic_to_ecef(lat_deg, lon_deg, h=0.0):
    """Geodetic (deg, deg, m) -> ECEF [m]. Accepts scalars or arrays."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    n = A / np.sqrt(1.0 - E2 * np.sin(lat) ** 2)
    x = (n + h) * np.cos(lat) * np.cos(lon)
    y = (n + h) * np.cos(lat) * np.sin(lon)
    z = (n * (1.0 - E2) + h) * np.sin(lat)
    return np.stack([x, y, z], axis=-1)


def ecef_to_geodetic(r):
    """ECEF [m] (..., 3) -> (lat_deg, lon_deg, h_m). Iterative, converges in a few steps."""
    r = np.asarray(r, dtype=float)
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1.0 - E2))
    h = np.zeros_like(lat)
    for _ in range(8):
        n = A / np.sqrt(1.0 - E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1.0 - E2 * n / (n + h)))
    return np.degrees(lat), np.degrees(lon), h


def up_vector(lat_deg, lon_deg):
    """Local ellipsoid normal (unit 'up') in ECEF."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) * np.ones_like(lon)],
        axis=-1,
    )


def gmst(jd_ut1):
    """Greenwich Mean Sidereal Time [rad] (IAU-82), vectorised over a Julian date array."""
    t = (np.asarray(jd_ut1) - 2451545.0) / 36525.0
    sec = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * t
        + 0.093104 * t**2
        - 6.2e-6 * t**3
    )
    return np.mod(np.radians(sec / 240.0), 2.0 * np.pi)


def teme_to_ecef(r_teme, v_teme, jd):
    """Rotate TEME position/velocity to ECEF (GMST only; polar motion ignored).

    Velocity is the velocity *relative to the rotating Earth*, i.e. what a
    stationary ground target sees: v_ecef = R (v_teme - omega x r_teme).
    """
    th = gmst(jd)
    c, s = np.cos(th), np.sin(th)
    x, y, z = r_teme[:, 0], r_teme[:, 1], r_teme[:, 2]
    vx, vy, vz = v_teme[:, 0], v_teme[:, 1], v_teme[:, 2]
    # omega x r with omega = (0, 0, w)
    vx_rel = vx + OMEGA_E * y
    vy_rel = vy - OMEGA_E * x
    r = np.stack([c * x + s * y, -s * x + c * y, z], axis=-1)
    v = np.stack([c * vx_rel + s * vy_rel, -s * vx_rel + c * vy_rel, vz], axis=-1)
    return r, v
