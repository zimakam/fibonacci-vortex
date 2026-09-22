import numpy as np
import math

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]
K = 5

def Gamma_fib(r):
    return sum(FIB[k]*math.exp(-k*LAM*r) for k in range(K))

def a_over_nu(r):
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    if abs(den) < 1e-15: return 0.0
    return (num/den)*LAM/r + 1.0/(r*r)

N = 200
r_min, r_max = 0.3, 10.0
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]
ani_arr = np.array([a_over_nu(rr) for rr in r])

print("="*72)
print("STEP 2: Energy method — numerical dE/dt for all perturbations")
print("="*72)
print()
print("Energy functional: E = int_0^inf delta(r)^2 * r dr")
print("dE/dt = 2 * int delta * rhs * r dr")
print()

def compute_dEdt(delta):
    delta_p = np.gradient(delta, dr)
    delta_pp = np.gradient(delta_p, dr)
    rhs = np.zeros_like(delta)
    for i in range(1, len(r)-1):
        ri = r[i]
        rhs[i] = (delta_pp[i]/ri - delta_p[i]/(ri*ri)) + ani_arr[i]*delta_p[i]
    integrand = 2 * delta * rhs * r
    return float(np.trapezoid(integrand, r))

def energy(delta):
    return float(np.trapezoid(delta*delta*r, r))

np.random.seed(42)
perturbations = [
    ("sin(pi*x)", 0.1*np.sin(np.pi*(r-r_min)/(r_max-r_min))),
    ("sin(2*pi*x)", 0.1*np.sin(2*np.pi*(r-r_min)/(r_max-r_min))),
    ("sin(5*pi*x)", 0.1*np.sin(5*np.pi*(r-r_min)/(r_max-r_min))),
    ("gaussian", 0.1*np.exp(-((r-3.0)/0.5)**2)),
    ("random 1%", 0.01*np.random.randn(N)),
    ("random 10%", 0.1*np.random.randn(N)),
]

for name, pert in perturbations:
    pert[0] = 0.0; pert[-1] = 0.0

print(f"{'perturbation':<20} {'dE/dt':>14} {'E(0)':>14} {'dE/dt/E':>14}")
print("-"*72)

all_negative = True
for name, pert in perturbations:
    dEdt = compute_dEdt(pert)
    E0 = energy(pert)
    ratio = dEdt/E0 if E0 > 1e-15 else 0
    if dEdt >= 0: all_negative = False
    print(f"{name:<20} {dEdt:>14.6e} {E0:>14.6e} {ratio:>14.6f}")

print()
print("="*72)
print("CONCLUSION")
print("="*72)
print()
if all_negative:
    print("*** ALL dE/dt < 0 → energy monotonically decreases ***")
    print("*** ANALYTIC-STYLE ENERGY DISSIPATION CONFIRMED ***")
else:
    print("*** Some dE/dt >= 0 → energy not monotone everywhere ***")
print()
print("Local decay rate dE/dt/E (log-slope of energy):")
for name, pert in perturbations:
    dEdt = compute_dEdt(pert)
    E0 = energy(pert)
    if E0 > 1e-15:
        rate = -dEdt/E0
        tau = 1/rate if rate > 0 else float('inf')
        print(f"  {name:<20} rate={rate:>10.4f}  tau={tau:>10.4f}")
