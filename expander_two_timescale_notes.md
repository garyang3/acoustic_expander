# Two-timescale analysis of the reed-valve acoustic expander

Model: Adams (2024), Eqs. 3.34–3.42. Goal: separate the simulated pressure trace into a slow mean and a fast oscillation,

$$P_{cv}(t) = \bar P(t) + \delta(t)\cos\phi(t)$$

and derive, in order, the evolution laws for the mean $\bar P(t)$, the amplitude $\delta(t)$, and the phase rate $\dot\phi = \omega(\bar P)$. Everything below is derived before it is used. Numerical checks refer to the validation run in `two_timescale_validation.png` (`full_model.py`, `averaged_model.py`).

Notation: overbars are cycle means, primes are fast fluctuations ($P = \bar P + P'$), hats are complex amplitudes ($P' = \mathrm{Re}[\hat p e^{i\phi}]$, $\delta = \lvert\hat p\rvert$, $a = \lvert\hat x\rvert$). Gas: $R = c_p - c_v$, $\gamma = c_p/c_v$. Angle brackets $\langle q \rangle$ denote the average over one oscillation cycle.

---

## 1. The model

Control volume (CV) with inlet/outlet reed valves; gas piston of mass $m_{gas}$, face area $A$, damping $b_{gas}$, riding between the CV and a sealed passive gas volume.

$$\frac{dM}{dt} = F_{in} - F_{out} \tag{1}$$

$$\frac{dE}{dt} = -P\frac{dV}{dt} + F_{in}c_pT_{in} - F_{out}c_pT \tag{2}$$

$$E = c_v M T \tag{3}$$

$$PV = MRT \tag{4}$$

$$m_{gas}\ddot x + b_{gas}\dot x = (P - P_{pass})A \tag{5}$$

$$V = V_{0c} + Ax \tag{6}$$

$$P_{pass}(x) = P_{p0}\left(\frac{V_{0p}}{V_{0p} - Ax}\right)^{\gamma} \tag{7}$$

$$F(\Delta P) = c_1 \dot m_{max}(1-u)\sqrt{u} \quad \text{for } 0 \lt u \lt 1, \quad \text{else } 0, \qquad u = \frac{\Delta P}{\Delta P_{max}}, \qquad c_1 = \frac{3\sqrt{3}}{2} \tag{8}$$

Here (1)–(4) are thesis 3.38, 3.34–3.36 (adiabatic CV, work $P dV/dt$), and 3.42; (5)–(6) are 3.39 and 3.41; (8) is the reed flow curve 3.37, applied as $F_{in} = F(P_{in}-P)$ and $F_{out} = F(P-P_{out})$; (7) treats the passive spring as a sealed adiabatic gas volume (its linearization is the thesis $k_{passive}$).

Note the two distinct "$P dV$" appearances in (2): the boundary work on the moving piston face, and the flow work hidden inside the port enthalpies — $h = u + Pv = c_pT$ for an ideal gas referenced to $T=0$, which is why $c_p$ (not $c_v$) multiplies the port terms.

---

## 2. First law + mass conservation reduced to a pressure equation

### 2.1 Two identities

From (3) and (4), eliminating $MT = PV/R$:

$$E = c_vMT = \frac{c_v}{R}PV \tag{9}$$

$$H = E + PV = Mc_pT = \frac{c_p}{R}PV \tag{10}$$

Equation (10) is the enthalpy-content identity you wanted justified. Dimension check: $c_p/R \approx 3.43$ for the thesis air properties ($c_p = 1050$, $R = 306$ in J per kg per K) is a pure number, and $PV$ is joules; the $1/R$ enters through $MT = PV/R$, so $H = c_p P V$ alone would not be dimensionally consistent. $H$ here is a state quantity of the CV; the quantities that flow through the ports in (2) are the specific enthalpies $c_pT_{in}$ and $c_pT$. Its cycle average needs the fast dynamics and is taken in section 8.

### 2.2 The pressure equation

Differentiate (9) and insert into (2):

$$\frac{c_v}{R}\dot P V + \frac{c_v}{R}P\dot V = -P\dot V + c_p(F_{in}T_{in} - F_{out}T)$$

Collect the $P\dot V$ terms:

$$\frac{c_v}{R}\dot P V = -\left(1 + \frac{c_v}{R}\right)P\dot V + c_p(F_{in}T_{in} - F_{out}T)$$

Use $1 + c_v/R = c_p/R$ and multiply by $R/(c_v V)$, using $(R/c_v)c_p = \gamma R$:

$$\frac{dP}{dt} = \frac{\gamma R}{V}(F_{in}T_{in} - F_{out}T) - \gamma\frac{P}{V}\frac{dV}{dt} \tag{11}$$

This single ODE is the first law for this system, with the enthalpy fluxes and the $P dV$ work explicit. Integrating the state $(M,P,x,\dot x)$ with (11) reproduces the thesis-literal $(M,E,x,\dot x)$ integration of (2) to better than one part in $10^{8}$ (`full_model.py` runs both).

### 2.3 The temperature equation

From $T = PV/MR$, logarithmic differentiation gives $\dot T/T = \dot P/P + \dot V/V - \dot M/M$. Insert (11) and (1), and use $1/M = RT/PV$:

$$\frac{dT}{dt} = \frac{RT}{PV}\left[F_{in}(\gamma T_{in} - T) - F_{out}(\gamma - 1)T\right] - (\gamma-1)\frac{T}{V}\frac{dV}{dt} \tag{12}$$

### 2.4 Limits worth memorizing

Sealed CV ($F_{in} = F_{out} = 0$): (11) and (12) reduce to $\dot P/P = -\gamma\dot V/V$ and $\dot T/T = -(\gamma-1)\dot V/V$, which integrate to

$$PV^{\gamma} = \mathrm{const}, \qquad TV^{\gamma-1} = \mathrm{const}, \qquad T \propto P^{(\gamma-1)/\gamma} \tag{13}$$

So the isentrope on your list is not an added assumption — it emerges whenever the valves are shut. Filling at $T_{in} = T$ heats the CV (flow work); blowdown cools it, $\dot T = -(\gamma-1)TF_{out}/M$ — the effect the expander exists to exploit.

Second law, briefly: $\dot S = F_{in}s_{in} - F_{out}s + \dot S_{gen}$ with $\dot S_{gen} \ge 0$. The sealed fast motion is isentropic by (13); entropy is generated only in the reed jets (throttling with the jet kinetic energy lost to turbulence, about $R\ln(P_{up}/P_{down})$ per unit mass), in blow-by (thesis 3.44), and in the piston damping. None of it is needed to close what follows; it is why the machine is sub-isentropic.

---

## 3. The two-timescale ansatz

### 3.1 Scale separation

Fast rate: $\omega \approx 820$ rad/s ($f \approx 130$ Hz). Slow rates in the same run: envelope growth or decay 2 to 4 per second, fractional mean drift below 1 per second, valve mass-exchange rate $F/M \approx 0.4$ per second (the valves pass about $F/(Mf) \approx 0.3$ percent of the CV mass per cycle). The separation parameter is

$$\varepsilon = \frac{\text{slow rate}}{\omega} \approx 10^{-2} \tag{14}$$

### 3.2 The ansatz (your step 2)

$$P(t) = \bar P(t) + \mathrm{Re}\left[\hat p(t)e^{i\phi(t)}\right], \qquad \dot\phi = \omega(\bar P) \tag{15}$$

subject to the slowness conditions that make the split well-defined — the envelope and the mean may change only by a fraction of order $\varepsilon$ per radian:

$$\left|\frac{1}{\hat p}\frac{d\hat p}{dt}\right| = O(\varepsilon\omega) \qquad \text{and} \qquad \left|\frac{1}{\bar P}\frac{d\bar P}{dt}\right| = O(\varepsilon\omega) \tag{16}$$

The same decomposition, with the same conditions, applies to $T$, $M$, $x$, $V$. Without (16) the split is empty (any signal can be written as "mean plus modulated cosine" in infinitely many ways); with it, $\bar P$ absorbs what varies slowly and $\hat p$ carries the envelope. The machinery is the method of averaging: project every model equation onto its cycle mean and onto its first harmonic. The valve nonlinearity (8) also generates harmonics $e^{2i\phi}$ and higher; the resonator responds to those off-resonance, suppressed by order $\varepsilon$, so first-harmonic closure is consistent — the 2 to 3 percent envelope error in the validation is their footprint.

### 3.3 Why there is no $e^{-ikx}$ factor

Your ansatz carried $e^{i\omega t - ikx}$. The device is acoustically compact: $c = \sqrt{\gamma R T} \approx 360$ m/s, so the wavelength at 130 Hz is $\lambda_{ac} \approx 2.8$ m against a 0.15 m neck, giving $kL \approx 0.34$ rad. Across any single element $e^{-ikx}$ is constant to that accuracy: the volumes carry spatially uniform $P'$ (this is exactly the thesis uniformity assumption) and the neck gas moves as a rigid slug — the Helmholtz gas piston. The spatial wave structure is replaced by the lumped mode shape, and $k$ survives only implicitly inside $\omega(\bar P)$. If the expander is later coupled to a long line or resonator tube, that element keeps its $e^{i(\omega t - kx)}$ and this lumped model becomes its boundary condition.

---

## 4. The fast subsystem

### 4.1 Within one cycle the CV is a sealed gas spring

Per cycle the valves exchange only about 0.3 percent of the CV mass (section 3.1), so on the fast timescale the valve terms in (11) and (12) are $O(\varepsilon)$ perturbations. At leading order the fast motion obeys the sealed relations (13). Linearize $PV^{\gamma} = \mathrm{const}$ about the cycle mean — expand $(\bar P + P')(\bar V + V')^{\gamma}$ to first order in the fluctuations:

$$\frac{P'}{\bar P} + \gamma\frac{V'}{\bar V} = 0 \qquad \Longrightarrow \qquad V' = -\frac{\bar V_{cv}}{\gamma\bar P}P' \tag{17}$$

This is the adiabatic gas-spring (compliance) relation; the $\gamma$ rather than 1 is compression heating, the same $\gamma P$ that appears in $c^2 = \gamma P/\rho$. Similarly $TV^{\gamma-1} = \mathrm{const}$ gives your step 3, the temperature decomposition:

$$\frac{T'}{\bar T} = \frac{\gamma-1}{\gamma}\frac{P'}{\bar P} \tag{18}$$

in phase with $P'$. Check against thesis Fig. 3.19: an amplitude of 0.25 bar on 5.7 bar at 300 K gives a 3.8 K swing — the plotted 294.5 to 302 K. With $V' = Ax'$, (17) converts piston to pressure amplitude:

$$\hat x = -\frac{\bar V_{cv}}{\gamma\bar P A}\hat p, \qquad \delta = \frac{k_{cv}a}{A}, \qquad k_{cv} = \frac{\gamma\bar P A^2}{\bar V_{cv}} \tag{19}$$

where $k_{cv}$ is the CV gas-spring constant (force $AP' = -k_{cv}x'$). Linearizing the passive adiabat (7) the same way gives the other spring:

$$k_{pass} = \frac{\gamma\bar P_{pass}A^2}{\bar V_{pass}} \tag{20}$$

### 4.2 Valve feedback on the fast pressure

Now keep the $O(\varepsilon)$ valve terms. Expand the flows about the slow pressure differences $\Delta\bar P_{in} = P_{in} - \bar P$ and $\Delta\bar P_{out} = \bar P - P_{out}$:

$$F_{in} = F(\Delta\bar P_{in} - P') \approx \bar F_{in} - F'_{in}P', \qquad F_{out} = F(\Delta\bar P_{out} + P') \approx \bar F_{out} + F'_{out}P' \tag{21}$$

with $F' = dF/d\Delta P$ — note the opposite signs: raising $P$ starves the inlet and feeds the outlet. Insert (21) into (11) and keep the first-harmonic terms:

$$\frac{dP'}{dt} = -\frac{k_{cv}}{A}\frac{dx'}{dt} + g P' \tag{22}$$

$$g = -\frac{\gamma R}{\bar V_{cv}}\left(T_{in}F'_{in} + \bar T F'_{out}\right) \tag{23}$$

The first term of (22) is the gas spring of (17) written with (19); the second is the valve feedback. Three small terms are dropped from (22), with their sizes in the validation run: the $T'$ feedback on the outflow enthalpy (about 0.2 per second, against $g$ of about 25 per second), the $V'$ modulation of the $1/V$ prefactors (relative size $\delta/\gamma\bar P$, about 1 percent), and slow-fast cross terms (order $\varepsilon^2$).

Differentiating the flow curve (8):

$$F'(\Delta P) = \frac{c_1\dot m_{max}}{\Delta P_{max}}\cdot\frac{1-3u}{2\sqrt{u}} \tag{24}$$

so $F' \gt 0$ for $u \lt 1/3$ (positive resistance: this valve damps, it makes $g$ more negative), $F' \lt 0$ for $1/3 \lt u \lt 1$ (negative resistance: it pumps, $g \gt 0$), and $F' = 0$ closed. This is thesis Fig. 3.13 and 3.15 in one formula.

### 4.3 Eigenvalues: frequency, growth rate, and a third root

The fast system is (22) plus the fluctuating part of the piston equation (5):

$$m_{gas}\ddot x' + b_{gas}\dot x' + k_{pass}x' = AP' \tag{25}$$

Eliminate $P'$: differentiate (25) once, substitute (22) for $\dot P'$, then use (25) again to replace $AP'$:

$$m_{gas}\dddot x' + (b_{gas} - g m_{gas})\ddot x' + (k_{cv} + k_{pass} - g b_{gas})\dot x' - g k_{pass}x' = 0$$

so $x' \propto e^{st}$ requires, writing $k = k_{cv} + k_{pass}$ and dropping subscripts on $m$, $b$:

$$m s^3 + (b - gm)s^2 + (k - gb)s - g k_{pass} = 0 \tag{26}$$

Look for the oscillatory pair perturbatively: $s = i\omega_0 + \lambda$ with $\omega_0^2 = k/m$ and $\lambda$, $b/m$, $g$ all of order $\varepsilon\omega_0$. To first order $s^2 = -\omega_0^2 + 2i\omega_0\lambda$ and $s^3 = -i\omega_0^3 - 3\omega_0^2\lambda$. Substitute into (26), discarding order $\varepsilon^2$ products ($gb$, $\lambda b$, $\lambda g$, $\lambda^2$):

imaginary part: $-m\omega_0^3 + k\omega_0 = 0$ — satisfied identically, so the frequency is unshifted at first order;

real part: $-3m\omega_0^2\lambda - (b - gm)\omega_0^2 + k\lambda - g k_{pass} = 0$, and with $k = m\omega_0^2$ this is $-2m\omega_0^2\lambda = (b-gm)\omega_0^2 + gk_{pass}$, giving

$$\omega = \omega_0(\bar P) = \sqrt{\frac{k_{cv}+k_{pass}}{m_{gas}}} = A\sqrt{\frac{\gamma\bar P}{m_{gas}}\left(\frac{1}{\bar V_{cv}} + \frac{1}{\bar V_{pass}}\right)} \tag{27}$$

$$\lambda = -\frac{b_{gas}}{2m_{gas}} + \frac{g}{2}\cdot\frac{k_{cv}}{k_{cv}+k_{pass}} \tag{28}$$

(the second form of (27) uses $\bar P_{pass} = \bar P$, derived at (36); it is the square-root-of-pressure chirp law — validation: 129.9 to 132.4 Hz predicted as the mean charges, 130.0 to 132.3 Hz simulated, 0.16 Hz rms). In (28), $b/2m = \omega_0/2Q$ is the resonator loss; the valve pumping is diluted by the lever ratio $k_{cv}/(k_{cv}+k_{pass})$ because the valves act only on the CV-side spring — one half for your equal volumes.

### 4.4 Cross-check of (28) by energy balance

Solve (22) for the harmonic amplitudes: $i\omega\hat p = -(k_{cv}/A)i\omega\hat x + g\hat p$, so

$$\hat p = -\frac{k_{cv}}{A}\hat x + i\frac{g k_{cv}}{\omega A}\hat x + O(g^2/\omega^2) \tag{29}$$

The in-phase part is the spring (19); the small quadrature part is the valves'. With $\langle q_1 q_2\rangle = \frac{1}{2}\mathrm{Re}[\hat q_1\hat q_2^{*}]$ and $\hat{\dot x} = i\omega\hat x$, the mean acoustic power the CV gas delivers to the piston is

$$\dot W_{ac} = A\langle P'\dot x'\rangle = \frac{A}{2}\mathrm{Re}\left[\hat p(i\omega\hat x)^{*}\right] = \frac{g k_{cv}}{2}a^2 \tag{30}$$

(the spring part is 90 degrees out of phase and averages to zero; only the quadrature part survives). The oscillator energy is $E_{osc} = \frac{1}{2}ka^2$, and its budget $dE_{osc}/dt = \dot W_{ac} - b\langle\dot x'^2\rangle$ with $\langle\dot x'^2\rangle = \frac{1}{2}\omega_0^2a^2$ gives $k a^2 \lambda = \frac{1}{2}g k_{cv}a^2 - \frac{1}{2}b\omega_0^2a^2$, i.e. exactly (28). So the growth rate is (valve acoustic power minus damping loss) over twice the stored energy.

### 4.5 The third root

The product of the roots of (26) is $g k_{pass}/m$; the oscillatory pair contributes $\omega_0^2$, so the remaining real root is

$$s_3 = \frac{g k_{pass}}{m\omega_0^2} = g\cdot\frac{k_{pass}}{k_{cv}+k_{pass}} = g\cdot\frac{\bar V_{cv}}{\bar V_{cv}+\bar V_{pass}} \tag{31}$$

This is not an acoustic mode: it is the drift rate of the mean pressure, and it reappears independently from the averaged equations (section 6.4). One gain, two instabilities: the same $g \gt 0$ that pumps the oscillation also tips the DC state, which is why the transient is simultaneously an amplitude transient and a mean ramp.

### 4.6 Equivalence with the thesis admittance criterion (3.32)

Onset, $\lambda \gt 0$, is $g \gt b\omega_0^2/k_{cv}$ by (28). Convert $g$ to a volumetric reed admittance ($\rho = \bar P/R\bar T$, $T_{in} \approx \bar T$): $g = -(\gamma\bar P/\bar V_{cv})Y_{reed}$ with $Y_{reed} = (F'_{in}+F'_{out})/\rho$, negative on the pumping branch. The input admittance of the resonator seen from the CV port is the CV compliance in parallel with the neck branch; its imaginary parts cancel precisely at (27), where for light damping its real part is $A^2\omega_0^2 b/k_{cv}^2$. Substituting both into $\lambda \gt 0$ collapses it to

$$-\mathrm{Re}[Y_{reed}] \gt \mathrm{Re}[Y_{res}(\omega_0)]$$

which is thesis Eq. 3.32 — here upgraded from a yes/no criterion to the quantitative rate (28).

---

## 5. Finite amplitude: describing functions for the valves

At your amplitudes the valves are swept across their whole curve and clip shut for part of each cycle, so the tangent slopes of (21) are not enough. The systematic replacement: for a valve seeing $\Delta P(t) = \Delta\bar P + \delta\cos\psi$, Fourier-project its actual flow waveform onto mean and first harmonic,

$$F_0(\Delta\bar P, \delta) = \frac{1}{2\pi}\int_0^{2\pi} F(\Delta\bar P + \delta\cos\psi)d\psi \tag{32}$$

$$F_1(\Delta\bar P, \delta) = \frac{1}{\pi}\int_0^{2\pi} F(\Delta\bar P + \delta\cos\psi)\cos\psi~d\psi, \qquad N = \frac{F_1}{\delta} \tag{33}$$

(the inlet sees the modulation with opposite sign, but both projections are parity-even, so the same functions serve both valves). Higher harmonics are non-resonant (section 3.2) and dropped. Then everywhere above and below:

$$F(\Delta\bar P) \to F_0 \quad \text{(slow equations)}, \qquad F'(\Delta\bar P) \to N \quad \text{(inside } g\text{)} \tag{34}$$

$$g(\bar P, \delta) = -\frac{\gamma R}{\bar V_{cv}}\left(T_{in}N_{in} + \bar T N_{out}\right) \tag{35}$$

Small-amplitude limit, by Taylor expansion under the integrals (odd powers of $\cos\psi$ average out; the mean of $\cos^2$ is one half and of $\cos^4$ is three eighths):

$$F_0 = F + \frac{1}{4}F''\delta^2 + O(\delta^4), \qquad N = F' + \frac{1}{8}F'''\delta^2 + O(\delta^4) \tag{36}$$

recovering section 4.2 as $\delta \to 0$. Because $F$ vanishes outside $0 \lt u \lt 1$, the integrals automatically cover both operating modes: blow-by bias ($\zeta \lt 1$, valves open at the mean) and peak-cracking bias ($\zeta \gt 1$, valves closed at the mean — $F$ and $F'$ both zero there — yet $F_0$ and $N$ nonzero because the swing reaches the open band near the pressure extremes). Your run is the second mode: the drive is a purely finite-amplitude effect that no tangent linearization could produce. ($\zeta$ is the thesis reed bias factor, Eq. 3.43.) These two integrals are the only non-elementary objects in the whole solution; `averaged_model.F0F1` evaluates them in a few lines.

---

## 6. The slow subsystem

Now average (1), (11), (5) over a cycle.

### 6.1 Mass

$$\frac{d\bar M}{dt} = F_{0,in} - F_{0,out} \tag{37}$$

The rectified flows (32), not the bare curve, drive the drift — at large $\delta$ they differ qualitatively (in your run the bare curve gives exactly zero).

### 6.2 The averaged work term

With $P = \bar P + P'$ and $\dot V = d\bar V/dt + \dot V'$:

$$\langle P\dot V\rangle = \bar P\frac{d\bar V}{dt} + \langle P'\dot V'\rangle = \bar P\frac{d\bar V}{dt} + \dot W_{ac} \tag{38}$$

the displacement work of the slowly drifting piston plus the acoustic power (30). This is the promised split of the $P dV$ work; converting (30) with (19), $\dot W_{ac} = g\bar V_{cv}\delta^2/(2\gamma\bar P)$. When $g \gt 0$ the mean state pays for the growing oscillation; during collapse ($g \lt 0$) the oscillation pays it back.

### 6.3 The quasi-static piston

Average the piston equation (5): the mean inertia and drag terms are of order $\varepsilon$ and $\varepsilon^2$ against the pressure scale, so at leading order the mean position continuously re-equilibrates, $\bar P_{pass}(\bar x) = \bar P$. Inverting the passive adiabat (7):

$$\bar V_{pass} = V_{0p}\left(\frac{P_{p0}}{\bar P}\right)^{1/\gamma}, \qquad \bar x = \frac{V_{0p} - \bar V_{pass}}{A}, \qquad \bar V_{cv} = V_{0c} + A\bar x \tag{39}$$

This drifting equilibrium is what your rolling-median baseline tracks. Since $\bar V_{cv} + \bar V_{pass} = V_{0c} + V_{0p}$ is constant, differentiating (39) gives

$$\frac{d\bar V_{pass}}{d\bar P} = -\frac{\bar V_{pass}}{\gamma\bar P}, \qquad \frac{d\bar V_{cv}}{d\bar P} = +\frac{\bar V_{pass}}{\gamma\bar P} \tag{40}$$

### 6.4 The averaged pressure equation

Average (11). The flux terms give $\gamma R(T_{in}F_{0,in} - \bar T F_{0,out})/\bar V_{cv}$ plus covariance corrections of relative order $\delta/\gamma\bar P$ (about 1 percent — the fluctuating parts of $1/V$, $T$, and the flows correlating with each other); the work term gives (38) times $\gamma/\bar V_{cv}$:

$$\frac{d\bar P}{dt} = \frac{\gamma R}{\bar V_{cv}}\left(T_{in}F_{0,in} - \bar T F_{0,out}\right) - \frac{\gamma\bar P}{\bar V_{cv}}\frac{d\bar V_{cv}}{dt} - \frac{g\delta^2}{2\bar P} \tag{41}$$

The middle term is not independent: by (40), $d\bar V_{cv}/dt = (\bar V_{pass}/\gamma\bar P)(d\bar P/dt)$, so it equals $(\bar V_{pass}/\bar V_{cv})(d\bar P/dt)$. Move it to the left and divide:

$$\frac{d\bar P}{dt} = \frac{\gamma R\left(T_{in}F_{0,in} - \bar T F_{0,out}\right) - \bar V_{cv}g\delta^2/(2\bar P)}{\bar V_{cv} + \bar V_{pass}} \tag{42}$$

The effective volume is the sum: the quasi-static piston transmits pressure, so slow charging fills both springs — with your equal volumes, half the naive rate. Linearizing (42) at small $\delta$ about the balance point $T_{in}F_{in} = \bar T F_{out}$ reproduces the drift rate $s_3$ of (31), as promised there. The set closes with the mean gas law $\bar T = \bar P\bar V_{cv}/(\bar M R)$.

### 6.5 The envelope equation

Two independent first-order effects add. The valve-and-damping balance gives $\dot a = \lambda a$ with $\lambda$ from (28), now evaluated with the finite-amplitude $g$ of (35). Separately, because $\omega(\bar P)$ drifts, the adiabatic invariant of a slowly varying oscillator ($E_{osc}/\omega$ constant, so $a \propto \omega^{-1/2}$) contributes $\dot a/a = -\dot\omega/2\omega$:

$$\frac{da}{dt} = \lambda(\bar P, \delta)a - \frac{\dot\omega}{2\omega}a, \qquad \delta = \frac{k_{cv}a}{A} \tag{43}$$

With a 2 percent chirp the second term is a 1 percent correction; the exponential of the integrated $\lambda$ carries the physics.

---

## 7. The assembled solution

Slow system — three ODEs, everything on the right explicit in the state $(\bar M, \bar P, a)$: equations (37), (42), (43), with the algebra (39), (27), (28), (35), (32)–(33), and $\bar T = \bar P\bar V_{cv}/\bar M R$. Fast reconstruction:

$$P_{cv}(t) = \bar P(t) + \delta(t)\cos\phi(t), \qquad \phi(t) = \phi_0 + \int_0^t\omega(\bar P(t'))dt' \tag{44}$$

with $T$, $x$, $V$ slaved by (18) and (19). Formally the envelope integrates to

$$\delta(t) = \delta(0)\cdot\frac{k_{cv}(t)}{k_{cv}(0)}\cdot\sqrt{\frac{\omega(0)}{\omega(t)}}\cdot\exp\left(\int_0^t\lambda(\bar P, \delta)dt'\right) \tag{45}$$

This is analytical in the same sense as a WKB or averaged van der Pol solution: the fast oscillation has been integrated out exactly to first order in $\varepsilon$, and what remains is smooth, non-stiff, and about a thousand times cheaper than the original system. In the operating phases the remaining integrals collapse to elementary functions:

**Ignition** (small $\delta$, coefficients frozen at $\bar P_0$): $\delta = \delta_0 e^{\lambda_0 t}$ and the mean departs its balance point as $e^{s_3 t}$, with $\lambda_0$ from (28) and $s_3$ from (31). Validation: predicted growth 4.0 per second, measured 3.5.

**Turnover**: the amplitude peaks where $\lambda(\bar P, \delta) = 0$ — the finite-amplitude form of thesis 3.32. In a saturating run this root is the limit cycle; in your charging run the mean rises through it (the inlet's open sliver slides from the negative- to the positive-resistance branch of (24)), $g$ changes sign, and collapse begins. Validation: predicted at $t = 0.126$ s and 5.90 bar; simulated peak at 0.12 s.

**Collapse**: $\delta \approx \delta_{pk}e^{-(t-t_{pk})/\tau_d}$ with $\tau_d = 1/\lvert\lambda\rvert$ — shorter than the bare ring-down $2Q/\omega$, because the inlet valve now actively absorbs the oscillation.

**Final approach of the mean** (outlet shut, $\delta$ small, opening small so $F_{0,in} \to c_1\dot m_{max}\sqrt{\Delta\bar P_{in}/\Delta P_{max}}$ by (8)): write $\Delta\bar P_{in} = P_{in} - \bar P$ and $V_{eff} = \bar V_{cv} + \bar V_{pass}$; then (42) reads

$$\frac{d\Delta\bar P_{in}}{dt} = -K\sqrt{\Delta\bar P_{in}}, \qquad K = \frac{\gamma R T_{in}c_1\dot m_{max}}{V_{eff}\sqrt{\Delta P_{max}}}$$

Substitute $y = \sqrt{\Delta\bar P_{in}}$: then $2y(dy/dt) = -Ky$, so $dy/dt = -K/2$, a linear decay of $y$, and

$$\Delta\bar P_{in}(t) = \Delta\bar P_{in}(t_1)\left(1 - \frac{t - t_1}{\tau_c}\right)^2, \qquad \tau_c = \frac{2V_{eff}\sqrt{\Delta\bar P_{in}(t_1)\Delta P_{max}}}{\gamma R T_{in}c_1\dot m_{max}} \tag{46}$$

a parabolic, finite-time landing on $P_{in}$ at $t_1 + \tau_c$ — not an exponential tail; the square-root orifice law is what makes the contact finite-time. Validation: from 0.021 bar at $t_1 = 0.25$ s, $\tau_c = 31$ ms, flat by about 0.28 s versus about 0.30 s simulated. (Convention caveat: with the check-valve clamp, $F = 0$ for $\Delta P \le 0$, a residual oscillation strands the mean a couple of hundredths of a bar above $P_{in}$; a valve that relieves backward lands exactly on $P_{in}$, as your notebook's trace does.)

**Ring-down** (both valves out of play): $g = 0$, so $\delta \propto e^{-\omega t/2Q}$ at constant $\bar P$. Validation: measured decay 2.03 per second against $\omega/2Q = 2.05$.

---

## 8. Cycle-averaged enthalpy (your step 4, completed)

The instantaneous identity is (10). Averaging it with $P = \bar P + P'$, $V = \bar V + V'$:

$$\langle H\rangle = \frac{c_p}{R}\left(\bar P\bar V + \langle P'V'\rangle\right) \tag{47}$$

Now use the fast relation (17), which makes $P'$ and $V'$ perfectly anti-correlated, and the average of $\cos^2$ equal to one half:

$$\langle P'V'\rangle = -\frac{\bar V_{cv}}{\gamma\bar P}\langle P'^2\rangle = -\frac{\bar V_{cv}\delta^2}{2\gamma\bar P} \qquad \Longrightarrow \qquad \frac{\langle H\rangle}{(c_p/R)\bar P\bar V_{cv}} = 1 - \frac{\delta^2}{2\gamma\bar P^2} \tag{48}$$

(The valves' order-$\varepsilon$ correction to (17) is the quadrature term of (29), which is 90 degrees out of phase with $V'$ and averages to zero here — the result is more robust than the leading-order derivation suggests.) The correction is about $10^{-4}$ at your amplitudes; the validation run measured a per-cycle residual of $-7.7$ against a predicted $-7.3$ (both times $10^{-5}$) at peak amplitude, correlated 1.0000 across all cycles. Practical consequence for your scatter test: regress $\langle H\rangle$ against the per-cycle product $\langle P\rangle\langle V\rangle$ — slope $c_p/R$, intercept zero, residuals $10^{-4}$. Plotted against $\langle P\rangle$ alone, the few-percent drift of $\bar V_{cv}$ from (39) masquerades as curvature; that, not the covariance, is the dominant deviation.

---

## 9. Validation summary and assumptions

Validation parameters: thesis Table 3.1 fluid (air, $c_p = 1050$, $c_v = 744$, $T_{in} = 300$ K, $P_{in} = 6$ bar, $P_{out} = 5.45$ bar); $V_{0c} = V_{0p} = 2\times10^{-4}$ m³, $A = 2.85\times10^{-4}$ m², $m_{gas} = 1$ g, $Q = 200$, $\Delta P_{max} = 0.20$ bar, $\dot m_{max} = 0.5$ g/s — a peak-cracking, charge-and-stall case like your run. Slow system versus full simulation of (1)–(8), no fitted parameters:

| Quantity | Two-timescale prediction | Full simulation |
|---|---|---|
| mean pressure, whole run | tracks | within 0.003 bar |
| envelope, growth through ring-down | tracks | within 2 to 3 percent |
| frequency chirp | 129.9 to 132.4 Hz | 130.0 to 132.3 Hz (0.16 Hz rms) |
| ignition growth rate | 4.0 per second | 3.5 per second |
| turnover (amplitude peak) | 0.126 s at 5.90 bar | 0.12 s |
| landing time | flat by 0.28 s | flat by 0.30 s |
| ring-down rate | 2.05 per second | 2.03 per second |
| enthalpy covariance (48) | $-7.3\times10^{-5}$ | $-7.7\times10^{-5}$, correlation 1.0000 |

Assumptions inherited from the thesis model: massless reeds, adiabatic uniform CV, ideal gas, fixed plenum pressures, amplitude-independent $Q$. Added by this analysis: scale separation (16) (holds at about $10^{-2}$ throughout; would degrade only if the valves were strong enough that $g$ approached $\omega$); first-harmonic closure (the 2 to 3 percent envelope error); quasi-static mean piston (order $\varepsilon^2$); and the small terms itemized after (22) and in section 6.4, all of relative order $\delta/\gamma\bar P$ or smaller. The check-valve convention $F = 0$ for $\Delta P \le 0$ matters only for the final few hundredths of a bar.
