import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris

from poliastro.twobody.propagation import propagate, cowell
from poliastro.ephem import build_ephem_interpolant
from poliastro.core.elements import rv2coe

from poliastro.constants import rho0_earth, H0_earth
from poliastro.core.perturbations import atmospheric_drag_exponential, third_body, J2_perturbation
from poliastro.bodies import Earth, Moon
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter3D
import plotly.io as pio
pio.renderers.default = "notebook_connected"


R = Earth.R.to(u.km).value
k = Earth.k.to(u.km ** 3 / u.s ** 2).value



C_D = 2.2
A_over_m = ((np.pi / 4.0) * (u.m ** 2) / (1000 * u.kg)).to_value(u.km ** 2 / u.kg)
B = C_D * A_over_m
rho0 = rho0_earth.to(u.kg / u.km ** 3).value
H0 = H0_earth.to(u.km).value



orbit = Orbit.circular(Earth, 200 * u.km, epoch=Time(0.0, format="jd", scale="tdb"))
tofs = TimeDelta(np.linspace(0 * u.h, 1000000 * u.d, num=200))

from poliastro.twobody.events import LithobrakeEvent
lithobrake_event = LithobrakeEvent(R)
events = [lithobrake_event]

rr = propagate(
        orbit,
        tofs,
        method=cowell,
        ad=atmospheric_drag_exponential,
        R=R,
        C_D=C_D,
        A_over_m=A_over_m,
        H0=H0,
        rho0=rho0,
        events=events,
    )

print('{0} {1:.2f} {2}'.format('orbital decay seen after', lithobrake_event.last_t.to(u.d).value, 'days'))


plt.ylabel("h(t)")
plt.xlabel("t, days")
plt.plot(tofs.value, rr.norm() - Earth.R)

