(* ::Package:: *)
(* Two-qubit resonator workflow: spin-basis Mathematica checks *)
(* This file mirrors the math steps in the research-library workflow note. *)

ClearAll["Global`*"];

(* ---------- Units and current workflow parameters ---------- *)

hbarUevUs = 6.582119569*10^-4;  (* micro-eV us; E/hbar gives rad/us *)
ghzToUev[fGHz_] := hbarUevUs*2 Pi*1000 fGHz;
mhzToRadPerUs[fMHz_] := 2 Pi fMHz;

parameterRules = {
  tc -> ghzToUev[10.0],
  bx -> 3.0,
  Bz -> 35.0,
  epsPrepare -> ghzToUev[80.0],
  wc -> 2 Pi*5.5*10^3,
  gc -> 2 Pi*50.0,
  hbar -> hbarUevUs
};

(* ---------- One-DQD spin x charge operators ---------- *)

id2 = IdentityMatrix[2];
sx = {{0, 1}, {1, 0}};
sy = {{0, -I}, {I, 0}};
sz = {{1, 0}, {0, -1}};

sigX = KroneckerProduct[sx, id2];
sigY = KroneckerProduct[sy, id2];
sigZ = KroneckerProduct[sz, id2];

tauX = KroneckerProduct[id2, sx];
tauY = KroneckerProduct[id2, sy];
tauZ = KroneckerProduct[id2, sz];

basisLabelsOneDQD = {"plusZ,R", "plusZ,L", "minusZ,R", "minusZ,L"};

hermitianQ[m_] := Chop[m - ConjugateTranspose[m]] ==
  ConstantArray[0, Dimensions[m]];

unitaryQ[u_] := Chop[ConjugateTranspose[u].u - IdentityMatrix[Length[u]]] ==
  ConstantArray[0, Dimensions[u]];

(* ---------- Lab-frame one-DQD Hamiltonian ---------- *)

H0[eps_] := tc tauX + (eps/2) tauZ + (Bz/2) sigZ + (bx/2) sigX.tauZ;

(* ---------- Spin-only diagonalization at finite detuning ---------- *)

Btot = Sqrt[Bz^2 + bx^2];
thetaS = ArcTan[Bz, bx];
Us = MatrixExp[I thetaS sigY.tauZ/2];

Hspin[eps_] := FullSimplify[Us.H0[eps].ConjugateTranspose[Us]];

HspinExpected[eps_] :=
  tc (Cos[thetaS] tauX - Sin[thetaS] sigY.tauY) +
  (eps/2) tauZ +
  (Btot/2) sigZ;

spinBasisHamiltonianCheck =
  FullSimplify[
    Hspin[eps] - HspinExpected[eps],
    {Bz > 0, bx >= 0, tc > 0, Element[eps, Reals]}
  ];

(* ---------- Operator transformations ---------- *)

tauZSpin = FullSimplify[Us.tauZ.ConjugateTranspose[Us]];
spinDriveSpin = FullSimplify[Us.sigX.ConjugateTranspose[Us]];

tauZTransformCheck = FullSimplify[tauZSpin - tauZ];

spinDriveTransformCheck =
  FullSimplify[spinDriveSpin - (Cos[thetaS] sigX + Sin[thetaS] sigZ.tauZ)];

(* ---------- Finite-detuning EDSR frequency ---------- *)

rightPlus = {1, 0, 0, 0};
rightMinus = {0, 0, 1, 0};

overlap[v_, ket_] := Abs[Conjugate[v].ket]^2;

edSrAngularFrequency[epsValue_, rules_: parameterRules] := Module[
  {h, vals, vecs, ord, plusIndex, minusIndex},
  h = N[Hspin[eps] /. eps -> epsValue /. rules];
  {vals, vecs} = Eigensystem[h];
  ord = Ordering[vals];
  vals = vals[[ord]];
  vecs = vecs[[ord]];
  plusIndex = First@Ordering[overlap[#, rightPlus] & /@ vecs, -1];
  minusIndex = First@Ordering[overlap[#, rightMinus] & /@ vecs, -1];
  Abs[vals[[plusIndex]] - vals[[minusIndex]]]/(hbar /. rules)
];

edSrFrequencyGHz[epsValue_, rules_: parameterRules] :=
  edSrAngularFrequency[epsValue, rules]/(2 Pi*1000);

currentOmegaEDSR =
  edSrAngularFrequency[epsPrepare /. parameterRules, parameterRules];

currentFreqEDSRGHz =
  edSrFrequencyGHz[epsPrepare /. parameterRules, parameterRules];

bareOmegaBz = (Bz/hbar) /. parameterRules;

(* ---------- Two-DQD plus resonator embedding ---------- *)

nCav = 10;
idC = IdentityMatrix[nCav, SparseArray];
idDQD = IdentityMatrix[4, SparseArray];

destroy[n_Integer] := SparseArray[
  Table[{i, i + 1} -> Sqrt[i], {i, 1, n - 1}],
  {n, n}
];

a = destroy[nCav];
adag = ConjugateTranspose[a];

opCav[op_] := KroneckerProduct[SparseArray[op], idDQD, idDQD];
opDQD1[op_] := KroneckerProduct[idC, SparseArray[op], idDQD];
opDQD2[op_] := KroneckerProduct[idC, idDQD, SparseArray[op]];

H01Lab[eps1_] := opDQD1[H0[eps1]]/hbar;
H02Lab[eps2_] := opDQD2[H0[eps2]]/hbar;
Hphoton := wc opCav[adag.a];

V1Lab := gc opDQD1[tauZ].opCav[a + adag];
V2Lab := gc opDQD2[tauZ].opCav[a + adag];

HfullLab[eps1_, eps2_] := Hphoton + H01Lab[eps1] + H02Lab[eps2] + V1Lab + V2Lab;

UfullSpinBasis := KroneckerProduct[idC, SparseArray[Us], SparseArray[Us]];

HfullSpinBasis[eps1_, eps2_] :=
  UfullSpinBasis.HfullLab[eps1, eps2].ConjugateTranspose[UfullSpinBasis];

(* ---------- Readout projectors and Bell targets ---------- *)

proj[v_] := Outer[Times, v, Conjugate[v]];

chargeR = {1, 0};
chargeL = {0, 1};
spinPlus = {1, 0};
spinMinus = {0, 1};

oneDQDRightProjector =
  KroneckerProduct[IdentityMatrix[2], proj[chargeR]];

RRProjectorFull =
  KroneckerProduct[
    IdentityMatrix[nCav, SparseArray],
    SparseArray[oneDQDRightProjector],
    SparseArray[oneDQDRightProjector]
  ];

bellPlusMinus =
  (KroneckerProduct[spinPlus, spinMinus] +
     KroneckerProduct[spinMinus, spinPlus])/Sqrt[2];

bellMinusPlus =
  (KroneckerProduct[spinPlus, spinMinus] -
     KroneckerProduct[spinMinus, spinPlus])/Sqrt[2];

bellFidelity[rhoSpin_, bell_: bellPlusMinus] :=
  Chop[Conjugate[bell].rhoSpin.bell];

rrWeight[rhoFullSpinBasis_] := Chop[Tr[RRProjectorFull.rhoFullSpinBasis]];

(* rhoFullSpinBasis must already be in the full spin-basis convention.
   Constructing the normalized RR-conditioned two-spin density matrix requires
   tracing out cavity and charge indices; the Python workflow currently performs
   that numerical partial trace. This file keeps the symbolic operators and
   validation checks explicit. *)

(* ---------- Minimal checks to evaluate in Mathematica ---------- *)

validationChecks := <|
  "H0Hermitian" -> hermitianQ[H0[eps]],
  "UsUnitary" -> unitaryQ[Us],
  "SpinBasisHamiltonian" -> spinBasisHamiltonianCheck,
  "TauZUnchanged" -> tauZTransformCheck,
  "SpinDriveTransform" -> spinDriveTransformCheck,
  "OmegaEDSRRadPerUs" -> currentOmegaEDSR,
  "FreqEDSRGHz" -> currentFreqEDSRGHz,
  "BareOmegaBzRadPerUs" -> bareOmegaBz
|>;

