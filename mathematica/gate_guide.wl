(* ::Package:: *)

(* ================================================================ *)
(*  Executable guide for single- and two-qubit gates                *)
(*  Written in a notebook-style flow to match RWA.nb                *)
(*  Basis: {|00>, |01>, |10>, |11>}                                 *)
(* ================================================================ *)

ClearAll["Global`*"];


(* ================================================================ *)
(*  Basic definitions                                               *)
(* ================================================================ *)

ket0 = {1, 0};
ket1 = {0, 1};

Id2 = IdentityMatrix[2];

\[Sigma]x = {{0, 1}, {1, 0}};
\[Sigma]y = {{0, -I}, {I, 0}};
\[Sigma]z = {{1, 0}, {0, -1}};

Tensor[A_, B_] := KroneckerProduct[A, B];

Op1[A_] := Tensor[A, Id2];
Op2[A_] := Tensor[Id2, A];

State00 = Tensor[ket0, ket0];
State01 = Tensor[ket0, ket1];
State10 = Tensor[ket1, ket0];
State11 = Tensor[ket1, ket1];

StateAssociation[\[Psi]_] := AssociationThread[
  {"00", "01", "10", "11"},
  Chop[\[Psi], 10^-12]
];

Propagator[H_, t_] := MatrixExp[-I H t];

GlobalPhaseDifference[U_, V_] := Module[
  {uFlat, vFlat, idx, phase},
  uFlat = Flatten[U];
  vFlat = Flatten[V];
  idx = FirstCase[
    Range[Length[vFlat]],
    k_ /; Abs[vFlat[[k]]] > 10^-12,
    Missing["NoReference"]
  ];
  If[idx === Missing["NoReference"], Return[U - V]];
  phase = uFlat[[idx]]/vFlat[[idx]];
  Chop[U - phase V, 10^-10]
];

EquivalentUpToGlobalPhaseQ[U_, V_, tol_: 10^-10] :=
  Max[Abs[Flatten[GlobalPhaseDifference[U, V]]]] < tol;


(* ================================================================ *)
(*  Single qubit RWA                                                *)
(* ================================================================ *)

(* From the notebook:

      H1 = {{-\[CapitalDelta]\[Omega]/2, \[CapitalOmega] Exp[I \[Phi]]},
            {\[CapitalOmega] Exp[-I \[Phi]], \[CapitalDelta]\[Omega]/2}}

   On resonance:
   \[Phi] = 0      -> X rotation
   \[Phi] = -Pi/2  -> Y rotation
*)

SingleQubitRWAHamiltonian[\[CapitalDelta]\[Omega]_, \[CapitalOmega]_, \[Phi]_] := {
  {-\[CapitalDelta]\[Omega]/2, \[CapitalOmega] Exp[I \[Phi]]},
  {\[CapitalOmega] Exp[-I \[Phi]], \[CapitalDelta]\[Omega]/2}
};

SingleQubitHamiltonian[
  qubit_Integer,
  \[CapitalDelta]\[Omega]_,
  \[CapitalOmega]_,
  \[Phi]_
] := Switch[
  qubit,
  1, Op1[SingleQubitRWAHamiltonian[\[CapitalDelta]\[Omega], \[CapitalOmega], \[Phi]]],
  2, Op2[SingleQubitRWAHamiltonian[\[CapitalDelta]\[Omega], \[CapitalOmega], \[Phi]]]
];

SingleQubitGate[
  qubit_Integer,
  axis_String,
  \[Theta]_,
  tGate_: 1
] := Module[
  {\[CapitalOmega], \[Phi]},
  Switch[
    ToUpperCase[axis],
    "X",
      \[CapitalOmega] = \[Theta]/(2 tGate);
      \[Phi] = 0;,
    "Y",
      \[CapitalOmega] = \[Theta]/(2 tGate);
      \[Phi] = -Pi/2;,
    "Z",
      Return[MatrixExp[-I (\[Theta]/2) If[qubit == 1, Op1[\[Sigma]z], Op2[\[Sigma]z]]]];
  ];
  Propagator[
    SingleQubitHamiltonian[qubit, 0, \[CapitalOmega], \[Phi]],
    tGate
  ]
];


(* ================================================================ *)
(*  Two qubit effective Hamiltonian                                 *)
(* ================================================================ *)

(* Same static effective Hamiltonian structure as in RWA.nb:

   Hqubit = 1/2 {{ 2 Ez,             0,              0, 0},
                 {    0, -J + \[Delta]Ez,           J, 0},
                 {    0,            J, -J - \[Delta]Ez, 0},
                 {    0,            0,              0, -2 Ez}}
*)

TwoQubitHamiltonian[Ez_, \[Delta]Ez_, J_] := (1/2) {
  {2 Ez, 0, 0, 0},
  {0, -J + \[Delta]Ez, J, 0},
  {0, J, -J - \[Delta]Ez, 0},
  {0, 0, 0, -2 Ez}
};


(* ================================================================ *)
(*  Diagonal form of the exchange block                             *)
(* ================================================================ *)

(* The anti-parallel block has eigenenergies

      (-J +/- Sqrt[J^2 + \[Delta]Ez^2]) / 2

   so the full diagonal Hamiltonian is:
*)

TwoQubitEigenHamiltonian[Ez_, \[Delta]Ez_, J_] := DiagonalMatrix[{
  Ez,
  (-J + Sqrt[J^2 + \[Delta]Ez^2])/2,
  (-J - Sqrt[J^2 + \[Delta]Ez^2])/2,
  -Ez
}];

EntanglingRate[Ez_, \[Delta]Ez_, J_] := Module[
  {E},
  E = Diagonal[TwoQubitEigenHamiltonian[Ez, \[Delta]Ez, J]];
  Simplify[E[[1]] - E[[2]] - E[[3]] + E[[4]]]
];


(* ================================================================ *)
(*  CZ from exchange                                                *)
(* ================================================================ *)

(* The nonlocal part behaves as

      HZZ = (J/4) \[Sigma]z \[CircleTimes] \[Sigma]z

   so for tCZ = Pi/J we get a ZZ phase gate.
   Adding local Z(-Pi/2) on each qubit turns it into CZ
   up to a global phase.
*)

HZZ[J_] := (J/4) Tensor[\[Sigma]z, \[Sigma]z];

UCZ[J_] := Module[
  {tCZ, UZZ, ULocal},
  tCZ = Pi/J;
  UZZ = Propagator[HZZ[J], tCZ];
  ULocal = SingleQubitGate[1, "Z", -Pi/2].SingleQubitGate[2, "Z", -Pi/2];
  Chop[ULocal.UZZ, 10^-12]
];

IdealCZ = DiagonalMatrix[{1, 1, 1, -1}];


(* ================================================================ *)
(*  Display the key Hamiltonians                                    *)
(* ================================================================ *)

Print["================================================"];
Print["Single-qubit rotating-frame Hamiltonian"];
Print["================================================"];

H1 = SingleQubitRWAHamiltonian[\[CapitalDelta]\[Omega], \[CapitalOmega], \[Phi]];
Print[H1 // MatrixForm];

Print[""];
Print["================================================"];
Print["Two-qubit effective Hamiltonian"];
Print["================================================"];

H2 = TwoQubitHamiltonian[Ez, \[Delta]Ez, J];
Print[H2 // MatrixForm];

Print[""];
Print["================================================"];
Print["Two-qubit eigenbasis Hamiltonian"];
Print["================================================"];

H2Diag = TwoQubitEigenHamiltonian[Ez, \[Delta]Ez, J];
Print[H2Diag // MatrixForm];

Print[""];
Print["Entangling rate E00 - E01 - E10 + E11 = ", EntanglingRate[Ez, \[Delta]Ez, J]];
Print["Therefore choose tCZ = Pi/J for the exchange pulse."];


(* ================================================================ *)
(*  Examples of individual gates                                    *)
(* ================================================================ *)

Print[""];
Print["================================================"];
Print["Example single-qubit gates"];
Print["================================================"];

UX90Q1 = SingleQubitGate[1, "X", Pi/2, 1];
UY90Q2 = SingleQubitGate[2, "Y", Pi/2, 1];
UZ90Q1 = SingleQubitGate[1, "Z", Pi/2];

Print["X(pi/2) on qubit 1"];
Print[UX90Q1 // MatrixForm];

Print[""];
Print["Y(pi/2) on qubit 2"];
Print[UY90Q2 // MatrixForm];

Print[""];
Print["Z(pi/2) on qubit 1"];
Print[UZ90Q1 // MatrixForm];


(* ================================================================ *)
(*  Check CZ                                                        *)
(* ================================================================ *)

Print[""];
Print["================================================"];
Print["CZ from exchange"];
Print["================================================"];

JExample = 1;
UCZExample = UCZ[JExample];

Print["CZ gate from exchange pulse"];
Print[UCZExample // MatrixForm];

Print[""];
Print["Equivalent to ideal CZ up to global phase?"];
Print[EquivalentUpToGlobalPhaseQ[UCZExample, IdealCZ]];


(* ================================================================ *)
(*  Full worked example                                             *)
(* ================================================================ *)

(* Sequence:
      1. single qubit gate on qubit 1
      2. two qubit CZ
      3. another single qubit gate

   Here we choose:
      U1 = X(pi/2) on qubit 1
      U2 = CZ
      U3 = Y(-pi/2) on qubit 2
*)

Print[""];
Print["================================================"];
Print["Worked example: single qubit -> CZ -> single qubit"];
Print["================================================"];

U1 = SingleQubitGate[1, "X", Pi/2, 1];
U2 = UCZ[JExample];
U3 = SingleQubitGate[2, "Y", -Pi/2, 1];

UTotal = Chop[U3.U2.U1, 10^-12];

\[Psi]Initial = State00;
\[Psi]After1 = Chop[U1.\[Psi]Initial, 10^-12];
\[Psi]After2 = Chop[U2.\[Psi]After1, 10^-12];
\[Psi]Final = Chop[UTotal.\[Psi]Initial, 10^-12];

Print["Initial state |00>"];
Print[StateAssociation[\[Psi]Initial]];

Print[""];
Print["After first single-qubit gate: X(pi/2) on qubit 1"];
Print[StateAssociation[\[Psi]After1]];

Print[""];
Print["After CZ"];
Print[StateAssociation[\[Psi]After2]];

Print[""];
Print["After final single-qubit gate: Y(-pi/2) on qubit 2"];
Print[StateAssociation[\[Psi]Final]];

Print[""];
Print["Total unitary for the sequence"];
Print[UTotal // MatrixForm];


(* ================================================================ *)
(*  Simple template for your own use                                *)
(* ================================================================ *)

Print[""];
Print["================================================"];
Print["Template"];
Print["================================================"];

Print["Use the following replacements in Mathematica:"];
Print["  1. choose your calibrated J"];
Print["  2. choose gate times for the single-qubit pulses"];
Print["  3. replace X/Y/Z choices as needed"];

(* Example template:

   JCal = 2 Pi 5.0*10^6;
   UPrep = SingleQubitGate[1, \"X\", Pi/2, t1];
   UEnt = UCZ[JCal];
   URead = SingleQubitGate[2, \"Z\", Pi/2];
   USequence = URead.UEnt.UPrep;

*)
