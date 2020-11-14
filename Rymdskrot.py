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
#Atmospheric drag
#The poliastro package now has several commonly used natural perturbations. One of them is atmospheric drag! See how one can monitor decay of the near-Earth orbit over time using our new module poliastro.twobody.perturbations!


R = Earth.R.to(u.km).value
k = Earth.k.to(u.km ** 3 / u.s ** 2).value


# parameters of a body
#C_D = 2.2  # dimentionless (any value would do)
#Min kommentar: ändrade från den betämda luftmotståndskonstanten C_D = 2 till att testa olika värden(värdena 1.02 till 2.0) i en vektor.
#Min kommentar: Dessutom lade jag till en for-loop för att programmet skulle ge mig svar för olika värden på C_D.
#C_D = np.array([1.22, 1.44, 1.66 ])
C_D = np.arange(1.50, 1.6, 0.1) #Min kommentar: Jag bytade från np.array till np.arange(börjar:1.5,slutar2.0,hoppar0.1)
for x in C_D:
    A_over_m = ((np.pi / 4.0) * (u.m ** 2) / (100 * u.kg)).to_value(u.km ** 2 / u.kg)  # km^2/kg #Min kommentar: Här tror jag vi antar att föremålet väger 100kg, men jag förstår inte area-beräkningen?
    B = x * A_over_m #Min kommentar: Skrev in x istället för C_D
    # parameters of the atmosphere
    rho0 = rho0_earth.to(u.kg / u.km ** 3).value  # kg/km^3
    H0 = H0_earth.to(u.km).value #Min kommentar: Ändras atmosfärens egenskaper ju högre radien är i det här programmet?



#[<matplotlib.lines.Line2D at 0x7ffb5f511048>]

#Orbital Decay
#If atmospheric drag causes the orbit to fully decay, additional code is needed to stop the integration when the satellite reaches the surface.

#Please note that you will likely want to use a more sophisticated atmosphere model than the one in atmospheric_drag for these sorts of computations.


    orbit = Orbit.circular(Earth, 230 * u.km, epoch=Time(0.0, format="jd", scale="tdb")) #Min kommentar: Här tror jag att jag kan ändra på värdet hos radien på omloppsbanan
    tofs = TimeDelta(np.linspace(0 * u.h, 100 * u.d, num=200))#Min kommentar: Här tror jag deltatiden beskrivs, att den mäter mellan noll timmar till ett max antal dagar, samt att den mäter 2000 gånger.

    from poliastro.twobody.events import LithobrakeEvent
    lithobrake_event = LithobrakeEvent(R)
    events = [lithobrake_event]

    rr = propagate(
        orbit,
        tofs,
        method=cowell,
        ad=atmospheric_drag_exponential,
        R=R,
        C_D=x, #Min kommentar: Här ändrade jag från C_D till x på högra sidan likhetstecknet
        A_over_m=A_over_m,
        H0=H0,
        rho0=rho0,
        events=events,
    )

    print('{0} {1:.2f} {2}'.format('orbital decay seen after', lithobrake_event.last_t.to(u.d).value, 'days'))

    plt.ylabel("h(t)")
    plt.xlabel("t, days")
    plt.plot(tofs.value, rr.norm() - Earth.R)

#[<matplotlib.lines.Line2D at 0x7ffb5f53f898>




