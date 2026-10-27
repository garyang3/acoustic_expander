# Section 5.1, expanded: substituting the split into the pressure equation

This note spells out every line that section 5.1 of `expander_packet_derivation.md` compresses into "substitute everything into (7) and multiply out". Equation numbers in square brackets refer to that file; plain numbers are this note's own. Notation follows the packet file: the flows are $\dot m_{in}$, $\dot m_{out}$, the slope of the valve curve is $F'$, and $\bar V$ is short for $\bar V_{cv}$. "The resonator fills" means the slow rise of $\bar P$.

**Where we start.** The pressure equation [7],

$$\frac{dP}{dt} = \frac{\gamma R}{V}\left(\dot m_{in}T_{in} - \dot m_{out}T\right) - \gamma\frac{P}{V}\frac{dV}{dt} \tag{1}$$

and the split [9] of every variable into a cycle mean plus a first harmonic with slowly varying amplitude:

$$P = \bar P + P', \qquad T = \bar T + T', \qquad V = \bar V + V', \qquad V' = Ax', \qquad \frac{dV}{dt} = \frac{d\bar V}{dt} + A\dot x' \tag{2}$$

The goal is two equations: one for the slow part $d\bar P/dt$ (handed to section 7) and one for the oscillating part $dP'/dt$, which is [14].

**The three rates in the problem.** Everything below sorts terms by size, so it helps to name the sizes first. The oscillation rate is $\omega \approx 800$ per second. The valve feedback rate $g$ (defined in step 5) is about 10 per second. The rates at which the mean state moves — the mean flows divided by the mass in the control volume, $\bar{\dot m}/\bar M$, and the relative rate of change of the mean volume, $(d\bar V/dt)/\bar V$ — are a few tenths per second. So $\omega \gg g \gg \bar{\dot m}/\bar M$, each step being a factor of 50 to 100 at your numbers. Section 5.1 keeps the first two rates in the oscillating equation and drops the third.

---

## Step 1: the ingredients, written out

**The volume prefactor.** Since $V' = Ax'$ is small compared with $\bar V$ (a few percent),

$$\frac{1}{V} = \frac{1}{\bar V + V'} = \frac{1}{\bar V}\cdot\frac{1}{1 + V'/\bar V} = \frac{1}{\bar V}\left(1 - \frac{V'}{\bar V} + \frac{V'^2}{\bar V^2} - \ldots\right) \tag{3}$$

**The flows.** Each valve sees a pressure drop that is its mean drop plus or minus the fluctuation: $\Delta P_{in} = P_{in} - P = \Delta\bar P_{in} - P'$ and $\Delta P_{out} = P - P_{out} = \Delta\bar P_{out} + P'$, with $\Delta\bar P_{in} = P_{in} - \bar P$ and $\Delta\bar P_{out} = \bar P - P_{out}$ as in the packet file. Expanding the valve curve $\dot m(\Delta P)$ about each mean drop to first order in $P'$ is [13]; with $F'_{in} \equiv F'(\Delta\bar P_{in})$ and $F'_{out} \equiv F'(\Delta\bar P_{out})$ the slopes of the valve curve at the two mean drops,

$$\dot m_{in} = \bar{\dot m}_{in} + \dot m'_{in}, \qquad \dot m'_{in} = -F'_{in}P' \tag{4}$$

$$\dot m_{out} = \bar{\dot m}_{out} + \dot m'_{out}, \qquad \dot m'_{out} = +F'_{out}P' \tag{5}$$

The signs are the point of (4) and (5): a rise in $P$ narrows the inlet's drop and widens the outlet's, so in the notation of [13], $d\dot m_{in}/dP' = -F'_{in}$ and $d\dot m_{out}/dP' = +F'_{out}$.

**The temperature fluctuation is not independent.** Within one cycle the gas is sealed to leading order, so $T$ follows the isentrope $T \propto P^{(\gamma-1)/\gamma}$ (from [3] and [8]); linearizing,

$$\frac{T'}{\bar T} = \frac{\gamma - 1}{\gamma}\cdot\frac{P'}{\bar P} \tag{6}$$

in phase with $P'$. This is used only to size a dropped term in step 5.

---

## Step 2: the flux term, multiplied out

Insert (3), (4), (5) and $T = \bar T + T'$ into the first term of (1). First the bracket:

$$\dot m_{in}T_{in} - \dot m_{out}T = \bar S + \left(\dot m'_{in}T_{in} - \dot m'_{out}\bar T\right) - \bar{\dot m}_{out}T' - \dot m'_{out}T', \qquad \bar S \equiv \bar{\dot m}_{in}T_{in} - \bar{\dot m}_{out}\bar T \tag{7}$$

where $\bar S$ is the mean enthalpy flux (per $c_p$). Now multiply by the prefactor (3) and sort by how many fluctuating factors each term contains:

$$\text{zeroth order:}\qquad \frac{\gamma R}{\bar V}\bar S \tag{8}$$

$$\text{first order:}\qquad \frac{\gamma R}{\bar V}\left(\dot m'_{in}T_{in} - \dot m'_{out}\bar T\right) - \frac{\gamma R}{\bar V}\bar{\dot m}_{out}T' - \frac{\gamma R}{\bar V^2}\bar S V' \tag{9}$$

$$\text{second order:}\qquad -\frac{\gamma R}{\bar V}\dot m'_{out}T' - \frac{\gamma R}{\bar V^2}\left(\dot m'_{in}T_{in} - \dot m'_{out}\bar T - \bar{\dot m}_{out}T'\right)V' + \frac{\gamma R}{\bar V^3}\bar S V'^2 \tag{10}$$

Substituting (4) and (5) into the first group of (9) turns it into a term proportional to $P'$:

$$\frac{\gamma R}{\bar V}\left(\dot m'_{in}T_{in} - \dot m'_{out}\bar T\right) = -\frac{\gamma R}{\bar V}\left(T_{in}F'_{in} + \bar TF'_{out}\right)P' \tag{11}$$

This is the valve feedback: the flows respond to the pressure fluctuation, and that response feeds back on the pressure fluctuation itself.

---

## Step 3: the work term, multiplied out

Insert (2) and (3) into the second term of (1):

$$-\gamma\frac{P}{V}\frac{dV}{dt} = -\frac{\gamma}{\bar V}\left(\bar P + P'\right)\left(1 - \frac{V'}{\bar V} + \ldots\right)\left(\frac{d\bar V}{dt} + A\dot x'\right) \tag{12}$$

Sorted the same way:

$$\text{zeroth order:}\qquad -\frac{\gamma\bar P}{\bar V}\frac{d\bar V}{dt} \tag{13}$$

$$\text{first order:}\qquad -\frac{\gamma\bar PA}{\bar V}\dot x' - \frac{\gamma}{\bar V}\frac{d\bar V}{dt}P' + \frac{\gamma\bar P}{\bar V^2}\frac{d\bar V}{dt}V' \tag{14}$$

$$\text{second order:}\qquad -\frac{\gamma A}{\bar V}P'\dot x' + \frac{\gamma\bar PA}{\bar V^2}V'\dot x' + \frac{\gamma}{\bar V^2}\frac{d\bar V}{dt}P'V' + \ldots \tag{15}$$

The first term of (14) is the gas spring: the piston's velocity changes the volume, and the sealed gas answers with a pressure rate.

---

## Step 4: sorting the equation into slow and oscillating parts

The left side of (1) splits as $dP/dt = d\bar P/dt + dP'/dt$. Every term on the right now belongs to one of three kinds: purely slow (no oscillating factor), first harmonic (exactly one oscillating factor), or a product of two oscillating factors. A product of two first harmonics is a constant plus a second harmonic,

$$\cos\psi\cos(\psi + \theta) = \tfrac12\cos\theta + \tfrac12\cos(2\psi + \theta) \tag{16}$$

so it contains no first harmonic at all. Terms of the three kinds cannot cancel one another over a cycle — a slow function, a $\cos\psi$ and a $\cos 2\psi$ are orthogonal — so the equation holds kind by kind. Collecting (8), (13) and the constant parts of (10) and (15) gives the slow equation,

$$\frac{d\bar P}{dt} = \frac{\gamma R}{\bar V}\bar S - \frac{\gamma\bar P}{\bar V}\frac{d\bar V}{dt} + \left\langle\text{second-order terms}\right\rangle \tag{17}$$

which is handed to section 7 (there $\bar{\dot m}$ becomes the rectified flow $F_0$ and the average of $-\gamma AP'\dot x'/\bar V$ becomes the acoustic power). Collecting the first-order terms (9), with (11), and (14) gives the oscillating equation before any dropping:

$$\frac{dP'}{dt} = -\frac{\gamma\bar PA}{\bar V}\dot x' - \frac{\gamma R}{\bar V}\left(T_{in}F'_{in} + \bar TF'_{out}\right)P' + \left[-\frac{\gamma R}{\bar V}\bar{\dot m}_{out}T' - \frac{\gamma R}{\bar V^2}\bar S V' - \frac{\gamma}{\bar V}\frac{d\bar V}{dt}P' + \frac{\gamma\bar P}{\bar V^2}\frac{d\bar V}{dt}V'\right] \tag{18}$$

The second-harmonic parts of (10) and (15) are not in (18): they cannot balance a first harmonic, and the resonator's response to them is a bounded ripple of relative size $g/\omega$ (see `higher_harmonics_check.ipynb`).

---

## Step 5: the size of the bracket, and the two terms that survive

Write each term of (18) as a rate multiplying $P'$ (using $\dot x' \approx \omega x'$ in magnitude, $V' = -\bar V P'/(\gamma\bar P)$ from [10], and (6) for $T'$):

$$\text{spring:}\quad \frac{\gamma\bar PA}{\bar V}\dot x' \sim \omega P' \tag{19}$$

$$\text{valve feedback:}\quad \frac{\gamma R}{\bar V}\left(T_{in}F'_{in} + \bar TF'_{out}\right)P' \equiv -gP' \tag{20}$$

$$\text{bracket:}\quad \left[-(\gamma - 1)\frac{\bar{\dot m}_{out}}{\bar M} + \frac{\bar S}{\bar T\bar M} - (\gamma + 1)\frac{1}{\bar V}\frac{d\bar V}{dt}\right]P' \tag{21}$$

(the ideal gas law $R\bar T/\bar V = \bar P/\bar M$ turns each bracket coefficient into a rate; the last two terms of the bracket in (18) combine because $V' \propto P'$). Every bracket rate is a mean-flow rate: $\bar{\dot m}_{out}/\bar M$, $\bar S/(\bar T\bar M) \approx (\bar{\dot m}_{in} - \bar{\dot m}_{out})/\bar M$, and $(d\bar V/dt)/\bar V$. At your numbers these are 0.1 to 0.3 per second. The valve feedback rate, by contrast, is built from the slope of the valve curve rather than the flow itself: using [15], $F' = (c_1\dot m_{max}/\Delta P_{max})f'(u)$ with $f'$ of order one, and $R\bar T/\bar V = \bar P/\bar M$,

$$g \sim \gamma~\frac{\bar P}{\Delta P_{max}}\cdot\frac{\dot m_{max}}{\bar M} \approx 10 \text{ per second} \tag{22}$$

The ratio of a bracket rate to $g$ is therefore about $(\Delta P_{max}/\bar P)\cdot(\bar{\dot m}/\dot m_{max})/\gamma \approx 0.08\times 0.2/1.4 \approx 10^{-2}$. The bracket is a one-percent correction to the valve feedback, and a $10^{-4}$ correction to the spring. This is the "numerically small" statement of section 5.1 made explicit, and it is the same size as the terms the averaging drops elsewhere, so keeping it would not improve the solution consistently. Dropping the bracket and writing $\gamma\bar PA/\bar V = k_{cv}/A$ from the definition $k_{cv} = \gamma\bar PA^2/\bar V$ in [10]:

$$\frac{dP'}{dt} = -\frac{k_{cv}}{A}\frac{dx'}{dt} + gP', \qquad g = -\frac{\gamma R}{\bar V}\left(T_{in}F'_{in} + \bar TF'_{out}\right) \tag{23}$$

which is [14]. The physical reading of the two kept terms: the spring term is the sealed gas answering the piston's motion; the $gP'$ term is the valves answering the pressure fluctuation — the outside world's only handle on the oscillation.

---

## Step 6: the slope of the valve curve

From the valve law [6], $\dot m = c_1\dot m_{max}(1 - u)\sqrt u$ with $u = \Delta P/\Delta P_{max}$, the chain rule gives $d\dot m/d\Delta P = (c_1\dot m_{max}/\Delta P_{max})\cdot d[(1-u)u^{1/2}]/du$, and

$$\frac{d}{du}\left[(1 - u)u^{1/2}\right] = \frac{1}{2}u^{-1/2} - \frac{3}{2}u^{1/2} = \frac{1 - 3u}{2\sqrt u} \qquad\Longrightarrow\qquad F'(\Delta P) = \frac{c_1\dot m_{max}}{\Delta P_{max}}\cdot\frac{1 - 3u}{2\sqrt u} \tag{24}$$

which is [15]. It is positive for $u \lt 1/3$ (the rising branch of the curve), zero at the peak $u = 1/3$, negative on the falling branch $1/3 \lt u \lt 1$, and zero where the valve is closed ($u \le 0$ or $u \ge 1$, where the curve is flat). Through (23), $g \gt 0$ — the valves make the oscillation grow — needs $T_{in}F'_{in} + \bar TF'_{out} \lt 0$, i.e. the two valves sitting, on balance, past the peaks of their curves.

---

## Step 7: what became of the second-order terms

They were not lost. Each product of two fluctuations has, by (16), a constant part and a second-harmonic part. The second-harmonic parts play no role at this order (step 4). The constant parts are cycle averages and belong to the slow equation (17); most of them are small in the same way as the bracket of step 5, and section 7 keeps the one that is not: from (15), the average of $-\gamma AP'\dot x'/\bar V$, which is $-(\gamma/\bar V)\langle P'\dot V'\rangle$ — the acoustic power, the net work per cycle that the gas does on the piston. The neighbouring term $(\gamma\bar PA/\bar V^2)V'\dot x'$ averages to zero, since $V'\dot x' = A x'\dot x'$ is the time derivative of $Ax'^2/2$ and a derivative has zero mean over a cycle.

---

## Two remarks on the packet file's own text

1. In the packet file's bracket of dropped terms the $\bar{\dot m}_{out}T'$ term is written with a plus sign; (9) shows it enters with a minus (an outflow of warmer gas removes more enthalpy), and the fourth term $\gamma\bar P(d\bar V/dt)V'/\bar V^2$ belongs with it. None of this changes the conclusion, since the whole bracket is dropped.
2. The packet file writes [13] with derivatives $d\dot m/dP'$ and [14] with the slopes $F'$. They are the same objects up to the sign bookkeeping of (4) and (5): $d\dot m_{in}/dP' = -F'_{in}$, $d\dot m_{out}/dP' = +F'_{out}$, which is what produces the plus sign inside $g$.
