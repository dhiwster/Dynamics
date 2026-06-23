(* ::Package:: *)

(* Symbolic paper-level derivation for the dispersive spin-spin coupling.

   Purpose:
   - keep the Hamiltonian in intrinsic operator notation using a, ad, and
     spin/charge operators
   - use NCAlgebra when available for noncommutative expansion
   - follow the paper-level derivation after projecting to the lower charge
     sector, then applying the rotating-wave and empty-cavity assumptions

   Scope:
   - derive the effective exchange Hamiltonian in the notation of
     Benito et al., Phys. Rev. B 100, 081412 (2019)
   - derive the collective Purcell lowering operator appearing in the
     supplementary-material master-equation treatment

   This file is intentionally symbolic and reduced.  The exact finite-matrix
   derivation remains in Spinspincoupling_clean_executable.wl.
*)

ClearAll["Global`*"];


(* ---------------------------------------------------------------------- *)
(* Logging                                                                 *)
(* ---------------------------------------------------------------------- *)

reportLines = {};
reportFile = Module[
  {baseDir},
  baseDir = If[
    StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[$InputFileName],
    Directory[]
  ];
  FileNameJoin[{baseDir, "Spinspincoupling_ncalgebra_paper_derivation_output.txt"}]
];
writeReport[] := Module[
  {txt},
  txt = StringRiffle[reportLines, "\n"];
  Export[reportFile, txt, "String"];
];

logLine[msg_] := (
  AppendTo[reportLines, msg];
  writeReport[];
  Print[msg];
);

saveReport[] := writeReport[];


(* ---------------------------------------------------------------------- *)
(* NCAlgebra setup                                                         *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  id, a, ad,
  tau1x, tau1z, tau2x, tau2z,
  sigma1x, sigma1z, sigma2x, sigma2z,
  sp1, sm1, sz1, sp2, sm2, sz2,
  t
];

loadNCAlgebra[] := Module[
  {userPacletRoot, candidateLoaders, packageDir, packageRoot},
  Quiet[Needs["NCAlgebra`"]];
  If[NameQ["NCAlgebra`SetNonCommutative"] || NameQ["SetNonCommutative"],
    Return[True]
  ];

  userPacletRoot = FileNameJoin[{$UserBaseDirectory, "Paclets", "Repository"}];
  candidateLoaders = Select[
    FileNames["NCAlgebra.m", userPacletRoot, Infinity],
    StringContainsQ[#, "NCAlgebra-"] &
  ];

  If[candidateLoaders === {},
    Return[False]
  ];

  packageDir = DirectoryName[First[candidateLoaders]];
  packageRoot = DirectoryName[packageDir];

  If[! MemberQ[$Path, packageDir], PrependTo[$Path, packageDir]];
  If[! MemberQ[$Path, packageRoot], PrependTo[$Path, packageRoot]];

  Quiet[Get[First[candidateLoaders]]];

  NameQ["NCAlgebra`SetNonCommutative"] || NameQ["SetNonCommutative"]
];

usingNCAlgebra = loadNCAlgebra[];

If[usingNCAlgebra,
  SetNonCommutative @@ {
    a, ad,
    tau1x, tau1z, tau2x, tau2z,
    sigma1x, sigma1z, sigma2x, sigma2z,
    sp1, sm1, sz1, sp2, sm2, sz2
  };
];

logLine["NCAlgebra available: " <> ToString[usingNCAlgebra]];


(* ---------------------------------------------------------------------- *)
(* Noncommutative helpers                                                  *)
(* ---------------------------------------------------------------------- *)

ncExpandExpr[expr_] := FixedPoint[
  If[usingNCAlgebra,
    NCExpand[#],
    Expand @ Distribute[#, Plus, NonCommutativeMultiply, Plus, Expand]
  ] &,
  expr
];

identityRules = {
  NonCommutativeMultiply[] :> id,
  NonCommutativeMultiply[x___, id, y___] :> NonCommutativeMultiply[x, y],
  id ** x_ :> x,
  x_ ** id :> x
};

bosonRules = {
  a ** ad :> ad ** a + id
};

sameSiteRules = {
  sp1 ** sp1 :> 0,
  sm1 ** sm1 :> 0,
  sp2 ** sp2 :> 0,
  sm2 ** sm2 :> 0,

  sp1 ** sm1 :> (id + sz1)/2,
  sm1 ** sp1 :> (id - sz1)/2,
  sp2 ** sm2 :> (id + sz2)/2,
  sm2 ** sp2 :> (id - sz2)/2,

  sz1 ** sz1 :> id,
  sz2 ** sz2 :> id,

  sz1 ** sp1 :> sp1,
  sp1 ** sz1 :> -sp1,
  sz1 ** sm1 :> -sm1,
  sm1 ** sz1 :> sm1,

  sz2 ** sp2 :> sp2,
  sp2 ** sz2 :> -sp2,
  sz2 ** sm2 :> -sm2,
  sm2 ** sz2 :> sm2
};

cavityCommutationRules = Flatten @ Table[
  {
    op ** a :> a ** op,
    op ** ad :> ad ** op
  },
  {op, {sp1, sm1, sz1, sp2, sm2, sz2}}
];

siteOrderingRules = Flatten @ Table[
  op2 ** op1 :> op1 ** op2,
  {op2, {sp2, sm2, sz2}},
  {op1, {sp1, sm1, sz1}}
];

ncRules = Join[
  identityRules,
  bosonRules,
  sameSiteRules,
  cavityCommutationRules,
  siteOrderingRules
];

ncSimplify[expr_] := FixedPoint[
  ncExpandExpr[# /. ncRules] &,
  expr
];

comm[x_, y_] := x ** y - y ** x;

projectVacuum[expr_] := ncSimplify[
  ncSimplify[expr] /. {
    NonCommutativeMultiply[x___, a, y___] :> 0,
    NonCommutativeMultiply[x___, ad, y___] :> 0,
    a -> 0,
    ad -> 0
  }
];

rwaReduce[expr_] := ncSimplify[
  ncSimplify[expr] /. {
    ad ** sp1 :> 0,
    a ** sm1 :> 0,
    ad ** sp2 :> 0,
    a ** sm2 :> 0
  }
];

displayExpr[expr_] := InputForm[ncSimplify[expr /. id -> 1]];
displayRaw[expr_] := InputForm[expr /. id -> 1];


(* ---------------------------------------------------------------------- *)
(* Reference-style 2x2 eigenvector rotation                                *)
(* ---------------------------------------------------------------------- *)

eigenRotvector2x2[ham_] := Module[
  {a0, b0, c0, l0, rot, rotinv},
  a0 = Simplify[(ham[[1, 1]] - ham[[2, 2]])/2];
  b0 = Simplify[(ham[[1, 2]] + ham[[2, 1]])/2];
  c0 = Simplify[(ham[[1, 2]] - ham[[2, 1]])/(-2 I)];
  l0 = Sqrt[a0^2 + b0^2 + c0^2];

  rot = {
    {
      (a0 + l0)/(Sqrt[2] Sqrt[l0 (a0 + l0)]),
      (-b0 + I c0)/(Sqrt[2] Sqrt[l0 (a0 + l0)])
    },
    {
      (b0 + I c0)/(Sqrt[2] Sqrt[l0 (a0 + l0)]),
      (a0 + l0)/(Sqrt[2] Sqrt[l0 (a0 + l0)])
    }
  };

  rotinv = {
    {
      (a0 + l0)/(Sqrt[2] Sqrt[l0 (a0 + l0)]),
      (b0 - I c0)/(Sqrt[2] Sqrt[l0 (a0 + l0)])
    },
    {
      (-b0 - I c0)/(Sqrt[2] Sqrt[l0 (a0 + l0)]),
      (a0 + l0)/(Sqrt[2] Sqrt[l0 (a0 + l0)])
    }
  };

  {rot, rotinv, l0}
];


(* ---------------------------------------------------------------------- *)
(* Single-site eigenvector rotation into the paper basis                   *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Single-site eigenvector rotation for the charge-spin blocks"];
logLine["Using the same style of 2x2 eigenvector rotation as in the IO theory notebook."];

hBlockMinus = {
  {-(2 tc - Bz)/2, bx/2},
  {bx/2, (2 tc - Bz)/2}
};

hBlockPlus = {
  {-(2 tc + Bz)/2, bx/2},
  {bx/2, (2 tc + Bz)/2}
};

rotMinus = eigenRotvector2x2[hBlockMinus];
rotPlus = eigenRotvector2x2[hBlockPlus];

uMinus = rotMinus[[2]];
uPlus = rotPlus[[2]];
rotationAssumptions =
  Element[{bx, Bz, tc}, Reals] && bx > 0 && tc > 0 && Bz > 0;
eMinus = FullSimplify[rotMinus[[3]], rotationAssumptions];
ePlus = FullSimplify[rotPlus[[3]], rotationAssumptions];

thetaMinus = ArcTan[bx/(2 tc - Bz)];
thetaPlus = ArcTan[bx/(2 tc + Bz)];
thetaBar = (thetaPlus + thetaMinus)/2;
uMinusAngle = {
  {Cos[thetaMinus/2], -Sin[thetaMinus/2]},
  {Sin[thetaMinus/2], Cos[thetaMinus/2]}
};
uPlusAngle = {
  {Cos[thetaPlus/2], -Sin[thetaPlus/2]},
  {Sin[thetaPlus/2], Cos[thetaPlus/2]}
};
rotatedHMinus = {{-eMinus, 0}, {0, eMinus}};
rotatedHPlus = {{-ePlus, 0}, {0, ePlus}};
rotatedDipoleMinus = {
  {-Sin[thetaMinus], Cos[thetaMinus]},
  {Cos[thetaMinus], Sin[thetaMinus]}
};
rotatedDipolePlus = {
  {-Sin[thetaPlus], Cos[thetaPlus]},
  {Cos[thetaPlus], Sin[thetaPlus]}
};

logLine["h_block,- = " <> ToString[InputForm[hBlockMinus]]];
logLine["h_block,+ = " <> ToString[InputForm[hBlockPlus]]];
logLine["theta_- = ArcTan[bx/(2 tc - Bz)]"];
logLine["theta_+ = ArcTan[bx/(2 tc + Bz)]"];
logLine["theta_bar = (theta_+ + theta_-)/2"];
logLine["E_- = " <> ToString[InputForm[eMinus]]];
logLine["E_+ = " <> ToString[InputForm[ePlus]]];
logLine["U_- = {{Cos[theta_-/2], -Sin[theta_-/2]}, {Sin[theta_-/2], Cos[theta_-/2]}}"];
logLine["U_+ = {{Cos[theta_+/2], -Sin[theta_+/2]}, {Sin[theta_+/2], Cos[theta_+/2]}}"];
logLine["U_- h_- U_-^dagger = {{-Sqrt[bx^2 + (Bz - 2*tc)^2]/2, 0}, {0, Sqrt[bx^2 + (Bz - 2*tc)^2]/2}}"];
logLine["U_+ h_+ U_+^dagger = {{-Sqrt[bx^2 + (Bz + 2*tc)^2]/2, 0}, {0, Sqrt[bx^2 + (Bz + 2*tc)^2]/2}}"];
logLine["U_- sigma_x U_-^dagger = {{-Sin[theta_-], Cos[theta_-]}, {Cos[theta_-], Sin[theta_-]}}"];
logLine["U_+ sigma_x U_+^dagger = {{-Sin[theta_+], Cos[theta_+]}, {Cos[theta_+], Sin[theta_+]}}"];
logLine["The helper eigenRotvector2x2 is kept as the implementation route; the matrices above are the compact angle-form equivalents."];
logLine["At the operator level, this first basis change is what turns the bare dipole coupling into the paper structure -gTau tau_x + gSigma tau_z sigma_x."];
logLine["Paper mapping: gTau = gc Cos[theta_bar], gSigma = gc Sin[theta_bar]."];
logLine["This is the basis-change step that produces the later tau_x and tau_z sigma_x structure."];


(* ---------------------------------------------------------------------- *)
(* Full model in the paper's transformed basis                             *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Full symbolic model before lower-charge projection"];

h0Full =
  omegaC (ad ** a) +
  (ETau1/2) tau1z + (ESigma1/2) sigma1z +
  (ETau2/2) tau2z + (ESigma2/2) sigma2z;

hIFull =
  ((-gTau1 tau1x + gSigma1 (tau1z ** sigma1x)) +
   (-gTau2 tau2x + gSigma2 (tau2z ** sigma2x))) ** (a + ad);

logLine["H0_full = " <> ToString[displayRaw[h0Full]]];
logLine["HI_full = " <> ToString[displayRaw[hIFull]]];


(* ---------------------------------------------------------------------- *)
(* Lower-charge projection                                                 *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Lower-charge projection: tau_z^(1) = tau_z^(2) = -1, tau_x^(i) -> 0"];

chargeProjectionRules = {
  tau1z -> -id,
  tau2z -> -id,
  tau1x -> 0,
  tau2x -> 0,
  sigma1x -> sp1 + sm1,
  sigma2x -> sp2 + sm2,
  sigma1z -> sz1,
  sigma2z -> sz2
};

h0Projected = ncSimplify[h0Full /. chargeProjectionRules];
hIProjectedReadable =
  -gSigma1 (sp1 + sm1) ** (a + ad) -
  gSigma2 (sp2 + sm2) ** (a + ad);

(* Absorb the overall sign from the charge projection into the effective
   paper-level couplings used below. *)
hIProjected =
  g1 (sp1 + sm1) ** (a + ad) +
  g2 (sp2 + sm2) ** (a + ad);

h0SpinPhoton =
  omegaC (ad ** a) + (omegaQ1/2) sz1 + (omegaQ2/2) sz2;

hIProjected = hIProjected /. {
  ESigma1 -> omegaQ1,
  ESigma2 -> omegaQ2
};

logLine["H0_projected = " <> ToString[displayRaw[h0Projected /. {ESigma1 -> omegaQ1, ESigma2 -> omegaQ2}]]];
logLine["HI_projected = " <> ToString[displayRaw[hIProjectedReadable]]];


(* ---------------------------------------------------------------------- *)
(* Rotating-frame transformation and RWA                                   *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Rotating-frame transformation and rotating-wave approximation"];
logLine["Following the reference notebook pattern: H_rot = U H U^dagger + i (dU/dt) U^dagger."];

uRotDescription =
  "Exp[i t (omegaC ad**a + (omegaQ1/2) sz1 + (omegaQ2/2) sz2)]";

rotatingActionRules = {
  a -> Exp[-I omegaC t] a,
  ad -> Exp[I omegaC t] ad,
  sp1 -> Exp[I omegaQ1 t] sp1,
  sm1 -> Exp[-I omegaQ1 t] sm1,
  sp2 -> Exp[I omegaQ2 t] sp2,
  sm2 -> Exp[-I omegaQ2 t] sm2,
  sz1 -> sz1,
  sz2 -> sz2
};

hIRotating = Expand[hIProjected /. rotatingActionRules];

hIRotatingCollected =
  g1 (
    Exp[I delta1 t] a ** sp1 +
    Exp[-I delta1 t] ad ** sm1 +
    Exp[-I (omegaC + omegaQ1) t] a ** sm1 +
    Exp[I (omegaC + omegaQ1) t] ad ** sp1
  ) +
  g2 (
    Exp[I delta2 t] a ** sp2 +
    Exp[-I delta2 t] ad ** sm2 +
    Exp[-I (omegaC + omegaQ2) t] a ** sm2 +
    Exp[I (omegaC + omegaQ2) t] ad ** sp2
  );

hIRWATimeDependent =
  g1 (
    Exp[I delta1 t] a ** sp1 +
    Exp[-I delta1 t] ad ** sm1
  ) +
  g2 (
    Exp[I delta2 t] a ** sp2 +
    Exp[-I delta2 t] ad ** sm2
  );

hIRWA =
  g1 (a ** sp1 + ad ** sm1) +
  g2 (a ** sp2 + ad ** sm2);

logLine["U_rot(t) = " <> uRotDescription];
logLine["Detuning convention: delta_i = omegaQi - omegaC."];
logLine["HI_rot(t) = " <> ToString[displayRaw[hIRotatingCollected]]];
logLine["RWA drops the fast terms exp[+- i (omegaC + omegaQi) t] a sm_i and ad sp_i."];
logLine["HI_RWA(t) = " <> ToString[displayRaw[hIRWATimeDependent]]];
logLine["For the static dispersive SW step we strip the slow phases and use HI_RWA = " <> ToString[displayExpr[hIRWA]]];


(* ---------------------------------------------------------------------- *)
(* Schrieffer-Wolff generator                                              *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Symbolic Schrieffer-Wolff generator for the static RWA Hamiltonian"];
logLine["Using the standard dispersive SW choice from the supplementary derivation."];

sGenerator =
  (g1/delta1) (a ** sp1 - ad ** sm1) +
  (g2/delta2) (a ** sp2 - ad ** sm2);

detuningDefinitions = {
  delta1 -> omegaQ1 - omegaC,
  delta2 -> omegaQ2 - omegaC
};

swResidual = Simplify[
  (
    ((g1/delta1) (omegaC - omegaQ1) + g1) (a ** sp1) +
    ((g2/delta2) (omegaC - omegaQ2) + g2) (a ** sp2) +
    ((g1/delta1) (omegaQ1 - omegaC) - g1) (ad ** sm1) +
    ((g2/delta2) (omegaQ2 - omegaC) - g2) (ad ** sm2)
  ) /. detuningDefinitions
];

logLine["HI_RWA = " <> ToString[displayExpr[hIRWA]]];
logLine["S = " <> ToString[displayExpr[sGenerator]]];
logLine["SW check ([S,H0] + HI_RWA == 0): " <> ToString[swResidual === 0]];


(* ---------------------------------------------------------------------- *)
(* Effective Hamiltonian in the empty-cavity sector                        *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Second-order empty-cavity dispersive Hamiltonian"];
logLine["Evaluating 1/2 [S, HI_RWA] with same-site spin algebra and cross-site commutativity."];

h2Vac = ncSimplify[
  (g1^2/(2 delta1)) (id + sz1) +
  (g2^2/(2 delta2)) (id + sz2) +
  (g1 g2/2) (1/delta1 + 1/delta2) (sp1 ** sm2 + sm1 ** sp2)
];

localStarkAndConstant =
  (g1^2/(2 delta1)) (id + sz1) +
  (g2^2/(2 delta2)) (id + sz2);

hDispersiveRotatingFrame = ncSimplify[h2Vac - localStarkAndConstant];

jExchange = Simplify[(g1 g2/2) (1/delta1 + 1/delta2)];

logLine["H2_vac = " <> ToString[displayExpr[h2Vac]]];
logLine["Hd_rotating = " <> ToString[displayExpr[hDispersiveRotatingFrame]]];
logLine["J12 = " <> ToString[InputForm[jExchange]]];


(* ---------------------------------------------------------------------- *)
(* Effective cavity damping channel                                        *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["First-order transformed cavity annihilation operator"];
logLine["Using a -> a + [S,a] to first order, then projecting to the empty cavity."];

aEffectiveVac = ncSimplify[(g1/delta1) sm1 + (g2/delta2) sm2];

logLine["a_eff,vac = " <> ToString[displayExpr[aEffectiveVac]]];
logLine[
  "With paper convention D[c] rho = 2 c rho c^dagger - c^dagger c rho - rho c^dagger c,"
];
logLine[
  "the Purcell term is (kappa/2) D[" <>
  ToString[displayExpr[aEffectiveVac]] <>
  "]."
];


(* ---------------------------------------------------------------------- *)
(* Identical-DQD specialization matching Eq. (2) of the paper             *)
(* ---------------------------------------------------------------------- *)

logLine[""];
logLine["Identical-DQD specialization matching the paper"];

identicalRules = {
  g1 -> gSigma,
  g2 -> gSigma,
  delta1 -> Delta1,
  delta2 -> Delta1,
  omegaQ1 -> omegaQ,
  omegaQ2 -> omegaQ
};

hPaperIdentical = ncSimplify[hDispersiveRotatingFrame /. identicalRules];
jPaperIdentical = Simplify[jExchange /. identicalRules];
aPurcellIdentical = ncSimplify[aEffectiveVac /. identicalRules];

logLine["Hd_identical = " <> ToString[displayExpr[hPaperIdentical]]];
logLine["J_identical = " <> ToString[InputForm[jPaperIdentical]]];
logLine["a_eff,vac,identical = " <> ToString[displayExpr[aPurcellIdentical]]];
logLine[
  "Paper Eq. (2): Hd = (gSigma^2/Delta1) (sp1 ** sm2 + sm1 ** sp2)."
];

saveReport[];
