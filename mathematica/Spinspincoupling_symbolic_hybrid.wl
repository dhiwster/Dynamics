(* ::Package:: *)

(* Readable symbolic front-end for the spin-charge-photon derivation.

   Purpose:
   - present the model in intrinsic operator notation using a, ad, tau, sigma
     and Mathematica's noncommutative product **
   - use NCAlgebra when available for symbolic noncommutative setup
   - keep the exact Schrieffer-Wolff calculation delegated to the verified
     matrix backend in Spinspincoupling_clean_executable.wl

   This is intentionally a wrapper, not a second symbolic algebra engine.
   The exact derivation is already correct and checked in the backend file.

   Run with

       wolframscript -script Spinspincoupling_symbolic_hybrid.wl
*)

ClearAll["Global`*"];


(* ---------------------------------------------------------------------- *)
(* NCAlgebra setup                                                         *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  a, ad,
  tau1x, tau1y, tau1z, tau2x, tau2y, tau2z,
  sigma1x, sigma1y, sigma1z, sigma2x, sigma2y, sigma2z,
  tau1zSigma1x, tau2zSigma2x
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

ncSymbols = {
  a, ad,
  tau1x, tau1y, tau1z, tau2x, tau2y, tau2z,
  sigma1x, sigma1y, sigma1z, sigma2x, sigma2y, sigma2z,
  tau1zSigma1x, tau2zSigma2x
};

If[usingNCAlgebra,
  SetNonCommutative @@ ncSymbols;
];

Print["NCAlgebra available: " <> ToString[usingNCAlgebra]];
If[usingNCAlgebra,
  Print["Using NCAlgebra symbolic layer with ** and SetNonCommutative."],
  Print["NCAlgebra not found. Falling back to Mathematica's built-in NonCommutativeMultiply (**)."]
];
Print[""];


(* ---------------------------------------------------------------------- *)
(* Symbolic presentation layer                                             *)
(* ---------------------------------------------------------------------- *)

symbolicSiteCoupling[site_Integer] := Switch[
  site,
  1, gSigma tau1zSigma1x - gTau tau1x,
  2, gSigma tau2zSigma2x - gTau tau2x
];

h0Sym =
  eTau (tau1z + tau2z) +
  eSigma (sigma1z + sigma2z) +
  omegaC (ad ** a);

hISym = (symbolicSiteCoupling[1] + symbolicSiteCoupling[2]) ** (a + ad);

Print["Readable symbolic Hamiltonians"];
Print["H0 = " <> ToString[InputForm[h0Sym]]];
Print["HI = " <> ToString[InputForm[hISym]]];
Print[""];
Print["Notation map to the paper:"];
Print["  E_tau = 2 eTau,  E_sigma = 2 eSigma,  omega_c = omegaC"];
Print["  tau1zSigma1x and tau2zSigma2x stand for tau_z sigma_x on each site."];
If[usingNCAlgebra,
  Print["  Example symbolic commutator target: a ** ad - ad ** a"];
];
Print[""];
Print["This file is the readable symbolic front-end only."];
Print["It now runs a symbolic paper-level NCAlgebra derivation first,"];
Print["then delegates the exact second-order SW calculation to the verified matrix backend."];
Print[""];


(* ---------------------------------------------------------------------- *)
(* Symbolic paper-level derivation                                         *)
(* ---------------------------------------------------------------------- *)

paperDerivationFile =
  FileNameJoin[{DirectoryName[$InputFileName], "Spinspincoupling_ncalgebra_paper_derivation.wl"}];

If[! FileExistsQ[paperDerivationFile],
  Print["Paper-derivation file not found: " <> paperDerivationFile];
  Exit[1];
];

Print["Running symbolic paper derivation: " <> paperDerivationFile];
Print[""];

Get[paperDerivationFile];


(* ---------------------------------------------------------------------- *)
(* Backend delegation                                                      *)
(* ---------------------------------------------------------------------- *)

backendFile = FileNameJoin[{DirectoryName[$InputFileName], "Spinspincoupling_clean_executable.wl"}];

If[! FileExistsQ[backendFile],
  Print["Backend file not found: " <> backendFile];
  Exit[1];
];

Print["Running backend: " <> backendFile];
Print[""];

Get[backendFile];
