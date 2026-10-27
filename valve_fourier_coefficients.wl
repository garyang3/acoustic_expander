(* ::Package:: *)

(* valve_fourier_coefficients.wl
   Fourier coefficients of the reed-valve flow over one oscillation cycle,
   companion to valve_fourier_coefficients.ipynb (same definitions, same test points).

   Valve law (dimensionless):  Fhat[u] = (1 - u) Sqrt[u] for 0 < u < 1, else 0,
   with u = ub + dh Cos[psi]   (ub = mean drop / dPmax,  dh = delta / dPmax).
   Physical coefficients are mdotmax * Fhat_n;  the thesis normalization multiplies all of them by c1 = 3 Sqrt[3]/2.

   Coefficients:  Fhat_0 = (1/(2 Pi)) Integrate[Fhat, {psi, 0, 2 Pi}],
                  Fhat_n = (1/Pi) Integrate[Fhat Cos[n psi], {psi, 0, 2 Pi}]   (n >= 1).

   Run:   wolframscript -file valve_fourier_coefficients.wl
   or     Get["valve_fourier_coefficients.wl"] inside a notebook.                                            *)


(* ------------------------------------------------------------------ 1. definitions ------------------------------------------------------------------ *)

Fhat[u_] := Piecewise[{{(1 - u) Sqrt[u], 0 < u < 1}}, 0];

(* numerical coefficient (any ub, dh > 0); the integrand is even about psi = Pi, so integrate over [0, Pi] and double.
   The points where the valve opens or closes are passed to NIntegrate as breakpoints. *)
FnNum[n_Integer, ub_?NumericQ, dh_?NumericQ] := Module[{pref, brk},
  pref = If[n == 0, 1/(2 Pi), 1/Pi];
  brk = Select[N[{ArcCos[-ub/dh], ArcCos[(1 - ub)/dh]}], Im[#] == 0 &];
  brk = DeleteDuplicates[Join[{0.}, Sort[Re /@ brk], {N[Pi]}], Abs[#1 - #2] < 10^-12 &];
  2 pref NIntegrate[Fhat[ub + dh Cos[psi]] Cos[n psi], Evaluate[Prepend[brk, psi]],
    MaxRecursion -> 30, AccuracyGoal -> 12, PrecisionGoal -> 12]
];

clippedQ[ub_, dh_] := (dh > ub) || (ub + dh > 1);


(* ---------------------------------------- 2. closed forms in the unclipped region (dh < ub, ub + dh < 1) ---------------------------------------- *)

(* Icoef[n, a, b] = (1/(2 Pi)) Integrate[Sqrt[a + b Cos[psi]] Cos[n psi], {psi, 0, 2 Pi}]
   Complete elliptic integrals of PARAMETER m = 2 b/(a + b):  EllipticE[m], EllipticK[m]  (same convention as scipy). *)
mpar[a_, b_] := 2 b/(a + b);

Ic0[a_, b_] := 2 Sqrt[a + b] EllipticE[mpar[a, b]]/Pi;
Ic1[a_, b_] := 2 Sqrt[a + b] (a EllipticE[mpar[a, b]] - (a - b) EllipticK[mpar[a, b]])/(3 Pi b);
Ic2[a_, b_] := 2 Sqrt[a + b] ((3 b^2 - 4 a^2) EllipticE[mpar[a, b]] + 4 a (a - b) EllipticK[mpar[a, b]])/(15 Pi b^2);
Ic3[a_, b_] := 2 Sqrt[a + b] ((32 a^3 - 29 a b^2) EllipticE[mpar[a, b]] -
     (32 a^3 - 32 a^2 b - 5 a b^2 + 5 b^3) EllipticK[mpar[a, b]])/(105 Pi b^3);

(* any n: recurrence  (p + 1) J[p + 1] == a (2 p + 1) J[p] - p (a^2 - b^2) J[p - 1]
   for J[p] = Integrate[(a + b Cos[psi])^p, {psi, 0, 2 Pi}], started from J[-1/2] and J[1/2],
   and Cos[n psi] = ChebyshevT[n, Cos[psi]] rewritten as a polynomial in w = a + b Cos[psi]. *)
Icoef[n_Integer, a_, b_] := Module[{J, p, w, coefs},
  J[-1/2] = 4 EllipticK[mpar[a, b]]/Sqrt[a + b];
  J[1/2] = 4 Sqrt[a + b] EllipticE[mpar[a, b]];
  Do[J[p + 1] = (a (2 p + 1) J[p] - p (a^2 - b^2) J[p - 1])/(p + 1), {p, 1/2, n - 1/2}];
  coefs = CoefficientList[Expand[ChebyshevT[n, (w - a)/b]], w];
  Sum[coefs[[k + 1]] J[k + 1/2], {k, 0, n}]/(2 Pi)
];

(* assemble the flow coefficients:  (1 - u) Sqrt[u] = (1 - ub) Sqrt[u] - dh Cos[psi] Sqrt[u],
   and Cos[psi] Cos[n psi] = (Cos[(n - 1) psi] + Cos[(n + 1) psi])/2 *)
FnClosed[0, ub_, dh_] := (1 - ub) Icoef[0, ub, dh] - dh Icoef[1, ub, dh];
FnClosed[n_Integer?Positive, ub_, dh_] := 2 (1 - ub) Icoef[n, ub, dh] - dh (Icoef[n - 1, ub, dh] + Icoef[n + 1, ub, dh]);

(* symbolic closed forms for F0, F1, F2 (valid for dh < ub, ub + dh < 1) *)
simp[expr_] := TimeConstrained[FullSimplify[expr, 0 < dh < ub && ub + dh < 1], 60, Simplify[expr]];
F0sym = simp[FnClosed[0, ub, dh]];
F1sym = simp[FnClosed[1, ub, dh]];
F2sym = simp[FnClosed[2, ub, dh]];
Print["F0 (closed form) = ", F0sym];
Print["F1 (closed form) = ", F1sym];
Print["F2 (closed form) = ", F2sym];

(* small-swing expansions (valid inside the open band):
   F0 ~ Fhat + Fhat''/4 dh^2,  F1 ~ Fhat' dh + Fhat''' dh^3/8,  F2 ~ Fhat''/4 dh^2,  F3 ~ Fhat''' dh^3/24 *)
fsm[u_] := (1 - u) Sqrt[u];
FnTaylor[ub_, dh_] := {fsm[ub] + fsm''[ub] dh^2/4, fsm'[ub] dh + fsm'''[ub] dh^3/8, fsm''[ub] dh^2/4, fsm'''[ub] dh^3/24};


(* ------------------------------------------------ 3. let Mathematica derive the elliptic forms itself ------------------------------------------------ *)

Print["\nMathematica's own Integrate for I0 and I1 (may take a minute; skipped if it exceeds 120 s):"];
symI = TimeConstrained[
  FullSimplify[{Integrate[Sqrt[a + b Cos[psi]], {psi, 0, 2 Pi}, Assumptions -> a > b > 0],
                Integrate[Sqrt[a + b Cos[psi]] Cos[psi], {psi, 0, 2 Pi}, Assumptions -> a > b > 0]}/(2 Pi)],
  120, $Failed];
If[symI === $Failed,
  Print["   timed out"],
  Print["   ", symI];
  Print["   difference from the stated forms at (a, b) = (0.5, 0.3): ",
        Chop[N[symI - {Ic0[a, b], Ic1[a, b]} /. {a -> 0.5, b -> 0.3}, 12]]]
];

(* clipped case: the swing leaves the open band, the integral runs over the psi-interval where 0 < u < 1.
   On [0, Pi] that interval is  ArcCos[Min[1, (1 - ub)/dh]] < psi < ArcCos[Max[-1, -ub/dh]].
   Symbolic evaluation gives incomplete elliptic integrals (EllipticE[phi, m], EllipticF[phi, m]); attempted here
   for the doubly clipped case with a time limit, otherwise use FnNum. *)
Print["\nSymbolic attempt, doubly clipped case (dh > ub, ub + dh > 1):"];
symClip = TimeConstrained[
  FullSimplify[(1/Pi) Integrate[(1 - ub - dh Cos[psi]) Sqrt[ub + dh Cos[psi]],
      {psi, ArcCos[(1 - ub)/dh], ArcCos[-ub/dh]}, Assumptions -> 0 < ub < 1 && dh > ub && ub + dh > 1]],
  120, $Failed];
Print["   F0 (clipped) = ", symClip];
If[symClip =!= $Failed,
  Print["   check at (ub, dh) = (0.6, 0.65): symbolic ", N[symClip /. {ub -> 0.6, dh -> 0.65}, 10],
        "  numeric ", FnNum[0, 0.6, 0.65]]];


(* ------------------------------------------------------------ 4. verification table ------------------------------------------------------------ *)

(* same test points as the Python notebook; its quadrature values are quoted for comparison *)
tests = {{0.5, 0.3}, {0.65, 0.2}, {0.3, 0.25}, {0.6, 0.39}, {0.65, 0.05}, {0.2, 0.02}};
pythonValues = {
  {0.311725, -0.093932, -0.042571, 0.004307},
  {0.267953, -0.115966, -0.014278, 0.000631},
  {0.331464,  0.049023, -0.054704, 0.010271},
  {0.249225, -0.182985, -0.061842, 0.006610},
  {0.281300, -0.029430, -0.000880, 0.000010},
  {0.357320,  0.008970, -0.000450, 0.000010}};
Print["\nunclipped test points:  {ub, dh}  |  F0..F3 numeric  |  max |closed - numeric|  |  max |Taylor - numeric|  |  Python quadrature"];
Do[
  With[{ub = tests[[i, 1]], dh = tests[[i, 2]]},
    Module[{num, cl, ty},
      num = Table[FnNum[n, ub, dh], {n, 0, 3}];
      cl = Table[N[FnClosed[n, ub, dh], 15], {n, 0, 3}];
      ty = N[FnTaylor[ub, dh]];
      Print[{ub, dh}, "  ", NumberForm[num, {7, 6}], "  ", ScientificForm[Max[Abs[num - cl]], 2], "  ",
            ScientificForm[Max[Abs[num - ty]], 2], "  ", pythonValues[[i]]]]],
  {i, Length[tests]}];


(* ------------------------------------------------------------ 5. your operating states ------------------------------------------------------------ *)

dPmax = (6.0*^5 - 5.45*^5)/(2*0.6);   (* Pa *)
mdotmax = 1.0*^-3;                    (* kg/s *)
states = {
  {"start, small swing",    0.545, 0.655, 0.05},
  {"start, after the kick", 0.545, 0.655, 0.15*^5/dPmax},
  {"steady oscillation",    (6.0*^5 - 5.7159*^5)/dPmax, (5.7159*^5 - 5.45*^5)/dPmax, 0.2945*^5/dPmax}};
Print["\nstate | valve | ub | dh | clipped | F0, F1, F2, F3 [g/s] | |F2/F1| | |F3/F1| | N = F1/delta [g/s per bar]"];
Do[
  With[{name = st[[1]], dh = st[[4]]},
    Do[
      With[{ub = st[[k + 1]], valve = {"inlet", "outlet"}[[k]]},
        Module[{c},
          c = 1000 mdotmax Table[FnNum[n, ub, dh], {n, 0, 3}];
          Print[name, " | ", valve, " | ", NumberForm[ub, {4, 3}], " | ", NumberForm[dh, {4, 3}], " | ", clippedQ[ub, dh], " | ",
                NumberForm[c, {6, 4}], " | ", NumberForm[Abs[c[[3]]/c[[2]]], {5, 3}], " | ", NumberForm[Abs[c[[4]]/c[[2]]], {5, 3}], " | ",
                NumberForm[c[[2]]/(dh dPmax/1.0*^5), {5, 3}]]]],
      {k, 1, 2}]],
  {st, states}];


(* ------------------------------------------------------------ 6. plots (exported next to this file) ------------------------------------------------------------ *)

dhs = Range[0.01, 1.2, 0.01];
panels = Table[
  Module[{data},
    data = Table[Table[{dh, FnNum[n, ub, dh]}, {dh, dhs}], {n, 0, 3}];
    ListLinePlot[data, PlotRange -> All, PlotLegends -> {"F0", "F1", "F2", "F3"},
      PlotLabel -> Row[{"ub = ", NumberForm[ub, {4, 3}]}], AxesLabel -> {"dh", "Fn"},
      GridLines -> {Select[{ub, 1 - ub}, 0 < # < 1.2 &], {0}}, ImageSize -> 320]],
  {ub, {0.2, 1/3, 0.55, 0.65, 0.9, 1.15}}];
grid = GraphicsGrid[Partition[panels, 3], ImageSize -> 1000];
Export["valve_fourier_coefficients_mathematica.png", grid];
Print["\nplot written: ", FileNameJoin[{Directory[], "valve_fourier_coefficients_mathematica.png"}]];
