import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris

from poliastro.twobody.propagation import propagate, cowell
from poliastro.ephem import build_ephem_interpolant
from poliastro.core.elements import rv2coe
from poliastro.atmosphere import COESA76

from poliastro.constants import rho0_earth, H0_earth
from poliastro.core.perturbations import atmospheric_drag_exponential, atmospheric_drag_model, third_body, J2_perturbation
from poliastro.bodies import Earth, Moon
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter3D
import plotly.io as pio
pio.renderers.default = "notebook_connected"
from poliastro.twobody.events import LithobrakeEvent

R = Earth.R.to(u.km).value

orbit = Orbit.circular(Earth, 300 * u.km)
t_decay = 7.17 * u.d

# parameters of a body
C_D = 2.2  # dimentionless (any value would do)
#A_over_m = ((np.pi / 4.0) * (u.m ** 2) / (100 * u.kg)).to_value(
#u.km ** 2 / u.kg
#)  # km^2/kg
A_over_m = (((4 * np.pi * 0.0144) / 4.0) * (u.m ** 2) / (48 * u.kg)).to_value(u.km ** 2 / u.kg)

tofs = [365] * u.d

lithobrake_event = LithobrakeEvent(R)
events = [lithobrake_event]

coesa76 = COESA76()

rr, _ = cowell(
    Earth.k,
    orbit.r,
    orbit.v,
    tofs,
    ad=atmospheric_drag_model,
    R=R,
    C_D=C_D,
    A_over_m=A_over_m,
    model=coesa76,
    events=events,
)

print('{0} {1:.2f} {2}'.format('orbital decay seen after', lithobrake_event.last_t.to(u.d).value, 'days'))
