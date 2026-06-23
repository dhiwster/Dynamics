(* ::Package:: *)

(* Executable derivation for the spin-charge-photon coupling note.

   This file replaces the hand-expanded second-order perturbation section
   with a matrix-based Schrieffer-Wolff calculation.  It can be run with

       wolframscript -file Spinspincoupling_clean_executable.wl

   Hilbert-space ordering:

       photon \[CircleTimes] tau1 \[CircleTimes] sigma1 \[CircleTimes] tau2 \[CircleTimes] sigma2

   The oscillator is truncated to {|0>, |1>}.  That is sufficient for the
   photon-vacuum second-order effective Hamiltonian because every second-order
   virtual process starting from |0> visits only |1>.
*)

ClearAll["Global`*"];


(* ---------------------------------------------------------------------- *)
(* Basic operators                                                        *)
(* ---------------------------------------------------------------------- *)

dPhoton = 2;

id[n_Integer] := IdentityMatrix[n];

sx = PauliMatrix[1];
sy = PauliMatrix[2];
sz = PauliMatrix[3];
sp = {{0, 1}, {0, 0}};
sm = {{0, 0}, {1, 0}};

a = SparseArray[Band[{1, 2}] -> Sqrt[Range[dPhoton - 1]], {dPhoton, dPhoton}];
ad = Transpose[a];
nPhoton = ad . a;

opPhoton[op_] := KroneckerProduct[op, id[2], id[2], id[2], id[2]];
opTau1[op_] := KroneckerProduct[id[dPhoton], op, id[2], id[2], id[2]];
opSigma1[op_] := KroneckerProduct[id[dPhoton], id[2], op, id[2], id[2]];
opTau2[op_] := KroneckerProduct[id[dPhoton], id[2], id[2], op, id[2]];
opSigma2[op_] := KroneckerProduct[id[dPhoton], id[2], id[2], id[2], op];

tau1x = opTau1[sx]; tau1y = opTau1[sy]; tau1z = opTau1[sz];
tau2x = opTau2[sx]; tau2y = opTau2[sy]; tau2z = opTau2[sz];
sigma1x = opSigma1[sx]; sigma1y = opSigma1[sy]; sigma1z = opSigma1[sz];
sigma2x = opSigma2[sx]; sigma2y = opSigma2[sy]; sigma2z = opSigma2[sz];

photonX = opPhoton[a + ad];
photonNumber = opPhoton[nPhoton];

comm[a_, b_] := a . b - b . a;


(* ---------------------------------------------------------------------- *)
(* Model                                                                  *)
(* ---------------------------------------------------------------------- *)

assumptions =
  Element[{eTau, eSigma, omegaC, gSigma, gTau}, Reals] &&
  eTau > 0 && eSigma > 0 && omegaC > 0;

siteCoupling[site_Integer] := Switch[
  site,
  1, gSigma tau1z . sigma1x - gTau tau1x,
  2, gSigma tau2z . sigma2x - gTau tau2x
];

h0 =
  eTau (tau1z + tau2z) +
  eSigma (sigma1z + sigma2z) +
  omegaC photonNumber;

hI = (siteCoupling[1] + siteCoupling[2]) . photonX;


(* ---------------------------------------------------------------------- *)
(* Schrieffer-Wolff generator                                             *)
(* ---------------------------------------------------------------------- *)

(* Solve [S, H0] == -HI element-by-element in the eigenbasis of H0.
   Since h0 is diagonal in this basis, S_ij = HI_ij/(E_i - E_j). *)

energies = Diagonal[h0];

swGenerator[v_] := Module[
  {dim, entries},
  dim = Length[v];
  entries = Table[
    If[i == j || v[[i, j]] === 0,
      0,
      Simplify[v[[i, j]]/(energies[[i]] - energies[[j]]), assumptions]
    ],
    {i, dim}, {j, dim}
  ];
  SparseArray[entries]
];

sGenerator = swGenerator[hI];

hEffSecondOrder = 1/2 comm[sGenerator, hI];
hEff = h0 + hEffSecondOrder;


(* ---------------------------------------------------------------------- *)
(* Photon-vacuum projection                                               *)
(* ---------------------------------------------------------------------- *)

projectPhotonVacuum[m_] := Module[
  {p0},
  p0 = KroneckerProduct[{{1, 0}}, id[16]];
  p0 . m . Transpose[p0]
];

hVac0 = projectPhotonVacuum[h0];
hVac2 = projectPhotonVacuum[hEffSecondOrder];
hVac = projectPhotonVacuum[hEff];


(* ---------------------------------------------------------------------- *)
(* Lower-charge projection                                                *)
(* ---------------------------------------------------------------------- *)

(* The paper's effective spin model works in the lower charge manifold
   tau_z^(1) = tau_z^(2) = -1, keeping only the spin degrees of freedom. *)

projectLowerCharge[m_] := Module[
  {pCharge},
  pCharge = KroneckerProduct[{{0, 1}}, id[2], {{0, 1}}, id[2]];
  pCharge . m . Transpose[pCharge]
];

hSpin0 = projectLowerCharge[hVac0];
hSpin2 = projectLowerCharge[hVac2];
hSpin = projectLowerCharge[hVac];


(* ---------------------------------------------------------------------- *)
(* Pauli-string decomposition                                             *)
(* ---------------------------------------------------------------------- *)

paulis = <|
  "I" -> id[2],
  "X" -> sx,
  "Y" -> sy,
  "Z" -> sz
|>;

pauliStrings[n_Integer] := pauliStrings[n] = Module[
  {keys, labels},
  keys = Keys[paulis];
  labels = Tuples[keys, n];
  AssociationThread[
    StringJoin /@ labels,
    KroneckerProduct @@ (paulis /@ #) & /@ labels
  ]
];

pauliDecompose[m_, n_Integer] := Module[
  {basis, coeffs, norm},
  basis = pauliStrings[n];
  norm = 2^n;
  coeffs = Association @ KeyValueMap[
    #1 -> Simplify[Together[Tr[#2 . m]/norm], assumptions] &,
    basis
  ];
  Select[coeffs, ! TrueQ[Simplify[Together[#] == 0, assumptions]] &]
];

pauliBasis4 = pauliStrings[4];
pauliBasis2 = pauliStrings[2];

hVac0Terms = pauliDecompose[hVac0, 4];
hVac2Terms = pauliDecompose[hVac2, 4];
hVacTerms = pauliDecompose[hVac, 4];

hSpin0Terms = pauliDecompose[hSpin0, 2];
hSpin2Terms = pauliDecompose[hSpin2, 2];
hSpinTerms = pauliDecompose[hSpin, 2];


(* ---------------------------------------------------------------------- *)
(* Self and cross pieces                                                  *)
(* ---------------------------------------------------------------------- *)

hI1 = siteCoupling[1] . photonX;
hI2 = siteCoupling[2] . photonX;
s1 = swGenerator[hI1];
s2 = swGenerator[hI2];

self1Vac2 = projectPhotonVacuum[1/2 comm[s1, hI1]];
self2Vac2 = projectPhotonVacuum[1/2 comm[s2, hI2]];
crossVac2 = projectPhotonVacuum[1/2 (comm[s1, hI2] + comm[s2, hI1])];

self1Terms = pauliDecompose[self1Vac2, 4];
self2Terms = pauliDecompose[self2Vac2, 4];
crossTerms = pauliDecompose[crossVac2, 4];

self1Spin2 = projectLowerCharge[self1Vac2];
self2Spin2 = projectLowerCharge[self2Vac2];
crossSpin2 = projectLowerCharge[crossVac2];

self1SpinTerms = pauliDecompose[self1Spin2, 2];
self2SpinTerms = pauliDecompose[self2Spin2, 2];
crossSpinTerms = pauliDecompose[crossSpin2, 2];


(* ---------------------------------------------------------------------- *)
(* Spin-only dispersive Hamiltonian                                       *)
(* ---------------------------------------------------------------------- *)

spinOp1[op_] := KroneckerProduct[op, id[2]];
spinOp2[op_] := KroneckerProduct[id[2], op];
spinExchange = KroneckerProduct[sp, sm] + KroneckerProduct[sm, sp];

spinIICoeff = Lookup[hSpin2Terms, "II", 0];
spinZICoeff = Lookup[hSpin2Terms, "ZI", 0];
spinIZCoeff = Lookup[hSpin2Terms, "IZ", 0];
spinZZCoeff = Lookup[hSpin2Terms, "ZZ", 0];
spinExchangeCoeff = Simplify[
  Lookup[hSpin2Terms, "XX", 0] + Lookup[hSpin2Terms, "YY", 0],
  assumptions
];

hSpinRWA2 = Simplify[
  spinIICoeff id[4] +
  spinZICoeff spinOp1[sz] +
  spinIZCoeff spinOp2[sz] +
  spinZZCoeff KroneckerProduct[sz, sz] +
  spinExchangeCoeff spinExchange,
  assumptions
];

hSpinRWA = Simplify[hSpin0 + hSpinRWA2, assumptions];
hSpinRWA2Terms = pauliDecompose[hSpinRWA2, 2];
hSpinRWATerms = pauliDecompose[hSpinRWA, 2];


(* ---------------------------------------------------------------------- *)
(* Consistency checks                                                     *)
(* ---------------------------------------------------------------------- *)

zeroQ[m_] := Module[
  {entries},
  entries = DeleteCases[DeleteDuplicates[Flatten[Normal[m]]], 0];
  AllTrue[entries, Simplify[Together[#] == 0, assumptions] &]
];

checks = <|
  "SW condition [S,H0] + HI == 0" -> zeroQ[comm[sGenerator, h0] + hI],
  "S is anti-Hermitian" -> zeroQ[sGenerator + ConjugateTranspose[sGenerator]],
  "self + cross reconstructs H2" -> zeroQ[self1Vac2 + self2Vac2 + crossVac2 - hVac2],
  "projected self + cross reconstructs spin H2" ->
    zeroQ[self1Spin2 + self2Spin2 + crossSpin2 - hSpin2],
  "Pauli decomposition reconstructs H2" ->
    zeroQ[
      Total[KeyValueMap[#2 pauliBasis4[#1] &, hVac2Terms]] - hVac2
    ],
  "Pauli decomposition reconstructs spin H2" ->
    zeroQ[
      Total[KeyValueMap[#2 pauliBasis2[#1] &, hSpin2Terms]] - hSpin2
    ]
|>;

If[! And @@ Values[checks],
  Print["One or more checks failed:"];
  Print[checks];
  Exit[1];
];


(* ---------------------------------------------------------------------- *)
(* Report                                                                 *)
(* ---------------------------------------------------------------------- *)

displayRules = {
  eTau -> Subscript[E, \[Tau]],
  eSigma -> Subscript[E, \[Sigma]],
  omegaC -> Subscript[\[Omega], c],
  gSigma -> Subscript[g, \[Sigma]],
  gTau -> Subscript[g, \[Tau]]
};

tex[expr_] := StringReplace[
  ToString[TeXForm[expr /. displayRules]],
  {
    "\n" -> " ",
    "\\text{c}" -> "c",
    "e_{\\tau }" -> "E_{\\tau}",
    "e_{\\sigma }" -> "E_{\\sigma}",
    "\\omega _c" -> "\\omega_c"
  }
];

operatorTex[axis_String, family_String, site_Integer] := Switch[
  axis,
  "I", "",
  "X", "\\" <> family <> "_x^{(" <> ToString[site] <> ")}",
  "Y", "\\" <> family <> "_y^{(" <> ToString[site] <> ")}",
  "Z", "\\" <> family <> "_z^{(" <> ToString[site] <> ")}"
];

labelTexFromSpec[label_String, spec_List] := Module[
  {axes, ops},
  axes = Characters[label];
  ops = MapThread[operatorTex[#1, #2[[1]], #2[[2]]] &, {axes, spec}];
  ops = DeleteCases[ops, ""];
  If[ops === {}, "I", StringRiffle[ops, " "]]
];

pauliLabelTex[label_String] := labelTexFromSpec[
  label,
  {{"tau", 1}, {"sigma", 1}, {"tau", 2}, {"sigma", 2}}
];

spinPauliLabelTex[label_String] := labelTexFromSpec[
  label,
  {{"sigma", 1}, {"sigma", 2}}
];

termTex[label_String, coeff_, labelRenderer_: pauliLabelTex] := Module[
  {op},
  op = labelRenderer[label];
  If[op === "I",
    tex[coeff],
    "\\left(" <> tex[coeff] <> "\\right) " <> op
  ]
];

hamiltonianTex[symbol_String, terms_Association, labelRenderer_: pauliLabelTex] := Module[
  {termList},
  termList = KeyValueMap[termTex[#1, #2, labelRenderer] &, terms];
  If[termList === {},
    symbol <> " &= 0",
    symbol <> " &= " <> First[termList] <>
      StringJoin[("\\\\\n&\\quad + " <> #) & /@ Rest[termList]]
  ]
];

printTerms[title_String, symbol_String, terms_Association, labelRenderer_: pauliLabelTex] := (
  Print[""];
  Print["================================================"];
  Print[title];
  Print["================================================"];
  Print["\\begin{aligned}"];
  Print[hamiltonianTex[symbol, terms, labelRenderer]];
  Print["\\end{aligned}"];
);

Print["Executable spin-spin-coupling derivation"];
Print["All symbolic checks passed."];

printTerms["Zeroth-order photon-vacuum Hamiltonian", "H_{0,\\mathrm{vac}}", hVac0Terms];
printTerms["Second-order photon-vacuum self term: site 1", "H_{2,\\mathrm{self}}^{(1)}", self1Terms];
printTerms["Second-order photon-vacuum self term: site 2", "H_{2,\\mathrm{self}}^{(2)}", self2Terms];
printTerms["Second-order photon-vacuum cross terms", "H_{2,\\mathrm{cross}}", crossTerms];
printTerms["Complete second-order photon-vacuum correction", "H_{2,\\mathrm{vac}}", hVac2Terms];
printTerms["Spin-only lower-charge Hamiltonian", "H_{\\mathrm{spin},0}", hSpin0Terms, spinPauliLabelTex];
printTerms["Spin-only lower-charge second-order correction", "H_{\\mathrm{spin},2}", hSpin2Terms, spinPauliLabelTex];
printTerms["Spin-only lower-charge cross term", "H_{\\mathrm{spin},2}^{\\mathrm{cross}}", crossSpinTerms, spinPauliLabelTex];
printTerms["Spin-only dispersive Hamiltonian after RWA reduction", "H_{\\mathrm{spin},\\mathrm{RWA}}", hSpinRWATerms, spinPauliLabelTex];

Print[""];
Print[
  "Exchange coefficient in paper form: J_disp = " <>
  tex[spinExchangeCoeff] <>
  " multiplying (\\sigma_+^{(1)} \\sigma_-^{(2)} + \\sigma_-^{(1)} \\sigma_+^{(2)})."
];
Print[""];
Print[
  "Use hVac, hVac0, hVac2, self1Vac2, self2Vac2, crossVac2, hSpin, hSpin2, crossSpin2, and hSpinRWA as executable outputs."
];
