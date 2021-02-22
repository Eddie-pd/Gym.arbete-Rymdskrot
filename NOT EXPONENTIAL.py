#The MIT License (MIT)

#Copyright (c) 2012-2019 Juan Luis Cano Rodríguez and the poliastro development team

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#"copies or substantial portions of the Software.


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


R = Earth.R.to(u.km).value
k = Earth.k.to(u.km ** 3 / u.s ** 2).value



C_D = 2.2
A_over_m = (((4 * np.pi * 0.0144) / 4.0) * (u.m ** 2) / (48 * u.kg)).to_value(u.km ** 2 / u.kg)
B = C_D * A_over_m
#rho0 = rho0_earth.to(u.kg / u.km ** 3).value
#H0 = H0_earth.to(u.km).value



orbit = Orbit.circular(Earth, 200 * u.km, epoch=Time(0.0, format="jd", scale="tdb"))
tofs = TimeDelta(np.linspace(0 * u.h, 100000000 * u.d, num=200))

from poliastro.twobody.events import LithobrakeEvent
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
        model = coesa76,
        events=events,
    )

print('{0} {1:.2f} {2}'.format('orbital decay seen after', lithobrake_event.last_t.to(u.d).value, 'days'))


#plt.ylabel("h(t)")
#plt.xlabel("t, days")
#plt.plot(tofs.value, rr, _.norm() - Earth.R)
