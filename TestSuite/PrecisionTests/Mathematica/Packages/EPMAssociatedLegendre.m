(* ::Package:: *)

BeginPackage["EPMAssociatedLegendre`"]


EPMALegendreGridSize::usage = 
"EPMALegendreGridSize[L] Get the size of the legendre grid.";


EPMALegendreGrid::usage =
"EPMALegendreGrid[L] give the legendre grid points.";


EPMALegendreWeights::usage =
"EPMALegendreWeights[L] give the Legendre quadrature weights.";


EPMALegendreP::usage =
"EPMALegendreP[L,M] get the associated Legendre polynomials on the grid.
EPMALegendreP[L,M,m] get the associated Legendre polynomials on the grid for a given m.";


EPMALegendreIntg::usage = 
"EPMALegendreIntg[L,M] Get the intg integrator . EPMWALegendreIntg[L,M,m] Get the intg integrator for a given m.";


EPMALegendreIntgTh2S::usage = 
"EPMALegendreIntgTh2S[L,M] Get the intgTh2S integrator . EPMWALegendreIntgTh2S[L,M,m] Get the intgTh2S integrator for a given m.";


EPMALegendreIntgTh2T::usage = 
"EPMALegendreIntgTh2T[L,M] Get the intgTh2T integrator . EPMWALegendreIntgTh2T[L,M,m] Get the intgTh2T integrator for a given m.";


EPMALegendreIntgPh2S::usage = 
"EPMALegendreIntgPh2S[L,M] Get the intgPh2S integrator . EPMWALegendreIntgPh2S[L,M,m] Get the intgPh2S integrator for a given m.";


EPMALegendreIntgPh2T::usage = 
"EPMALegendreIntgPh2T[L,M] Get the intgPh2T integrator . EPMWALegendreIntgPh2T[L,M,m] Get the intgPh2T integrator for a given m.";


EPMALegendreProj::usage = 
"EPMALegendreProj[L,M] Get the proj projector. EPMALegendreProj[L,M,m] Get the proj projector for a given m.";


EPMALegendreDProj::usage = 
"EPMALegendreDProj[L,M] Get the dProj projector. EPMALegendreDProj[L,M,m] Get the dProj projector for a given m.";


EPMALegendreProjS2Th::usage = 
"EPMALegendreProjS2Th[L,M] Get the projS2Th projector. EPMALegendreProjS2Th[L,M,m] Get the projS2Th projector for a given m.";


EPMALegendreProjT2Th::usage = 
"EPMALegendreProjT2Th[L,M] Get the projT2Th projector. EPMALegendreProjT2Th[L,M,m] Get the projT2Th projector for a given m.";


EPMALegendreProjS2Ph::usage = 
"EPMALegendreProjS2Ph[L,M] Get the projS2Ph projector. EPMALegendreProjS2Ph[L,M,m] Get the projS2Ph projector for a given m.";


EPMALegendreProjT2Ph::usage = 
"EPMALegendreProjT2Ph[L,M] Get the projT2Ph projector. EPMALegendreProjT2Ph[L,M,m] Get the projT2Ph projector for a given m.";


EPMALegendreProj2GradPh::usage = 
"EPMALegendreProj2GradPh[L,M] Get the proj2GradPh projector. EPMALegendreProj2GradPh[L,M,m] Get the proj2GradPh projector for a given m.";


EPMALegendreProj2GradTh::usage = 
"EPMALegendreProj2GradTh[L,M] Get the proj2GradTh projector. EPMALegendreProj2GradTh[L,M,m] Get the proj2GradTh projector for a given m.";


Begin["`Private`"]


(*****************  Set precisions *******************)
highPrecision=80;
epmPrecision=16;


(*****************  Compute grid size *******************)
alGridN[l_] = Floor[(l+1)3/2];
EPMALegendreGridSize[l_] := alGridN[l];


(*****************  Setup general normalisation function *******************)
schmidtNorm[l_,m_] = If[l==0 && m ==0, 1,If[m==0, Sqrt[(l-m)!/(l+m)!],Sqrt[2] Sqrt[(l-m)!/(l+m)!]]];


weightNorm[l_] = (2 l+1)/4;


(*****************  Create private polynomials and derivatives computations *******************)
alPoly[l_Integer, m_Integer, x_]:= schmidtNorm[l,m]LegendreP[l,m,x];
alDPoly[l_Integer, m_Integer, x_]:= schmidtNorm[l,m]Evaluate[-Sin[ArcCos[x]] D[LegendreP[l,m,t],t]/.{t->x}];
alProj[l_Integer, m_Integer, x_]:= alPoly[l,m,x];
alProjS2Th[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0] alDPoly[l,m,x];
alProjT2Th[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0] 1/Sin[ArcCos[x]] alPoly[l,m,x];
alProjS2Ph[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0] 1/Sin[ArcCos[x]] alPoly[l,m,x];
alProjT2Ph[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0]alDPoly[l,m,x];
alProj2GradTh[l_Integer, m_Integer, x_]:= alDPoly[l,m,x];
alProj2GradPh[l_Integer, m_Integer, x_]:= 1/Sin[ArcCos[x]] alDPoly[l,m,x];
alIntg[l_Integer, m_Integer, x_]:=  weightNorm[l]alPoly[l,m,x];
alIntgTh2S[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0]weightNorm[l]alDPoly[l,m,x];
alIntgTh2T[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0] 1/Sin[ArcCos[x]] weightNorm[l]alPoly[l,m,x];
alIntgPh2S[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0] 1/Sin[ArcCos[x]] weightNorm[l]alPoly[l,m,x];
alIntgPh2T[l_Integer, m_Integer, x_]:= If[l!= 0,1/Sqrt[l (l+1)],0]weightNorm[l]alDPoly[l,m,x];


(*****************  Create high precision grid values *******************)
alGrid[l_] := Sort[N[x/.Solve[LegendreP[alGridN[l],x] ==0, x], highPrecision],Less]


(*****************  Create high precision weights *******************)
alWeights[l_] := Module[{gridN = alGridN[l], grid = alGrid[l]},Table[N[(2(1-grid[[i]]^2))/((gridN+1)^2 LegendreP[gridN+1, grid[[i]]]^2), highPrecision], {i,1,gridN}]]


(*****************  Compute EPM precision grid *******************)
EPMALegendreGrid[l_] := N[alGrid[l], epmPrecision];


(*****************  Compute EPM precision weights *******************)
EPMALegendreWeights[l_] := N[alWeights[l], epmPrecision];


(*****************  Compute EPM precision Associated Legendre polynomials *******************)
EPMALegendreP[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alPoly[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreP[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alPoly[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProj[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProj[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProj[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProj[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreDProj[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alDPoly[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreDProj[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alDPoly[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProjS2Ph[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjS2Ph[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProjS2Ph[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjS2Ph[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProjT2Ph[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjT2Ph[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProjT2Ph[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjT2Ph[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProjS2Th[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjS2Th[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProjS2Th[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjS2Th[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProjT2Th[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjT2Th[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProjT2Th[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProjT2Th[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProj2GradPh[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProj2GradPh[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProj2GradPh[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProj2GradPh[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre projector *******************)
EPMALegendreProj2GradTh[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProj2GradTh[i,k,grid[[j]]], {k, 0, m},{i, k, l},{j, 1, gridN}], epmPrecision]];
EPMALegendreProj2GradTh[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l]},N[Table[alProj2GradTh[i,m,grid[[j]]],{i, m, l},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre integrator *******************)
EPMALegendreIntg[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntg[i,k,grid[[j]]], {k, 0, m},{j, 1, gridN},{i, k, l}], epmPrecision]];
EPMALegendreIntg[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntg[i,m,grid[[j]]],{j, 1, gridN},{i, m, l}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre integrator *******************)
EPMALegendreIntgTh2S[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgTh2S[i,k,grid[[j]]], {k, 0, m},{j, 1, gridN},{i, k, l}], epmPrecision]];
EPMALegendreIntgTh2S[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgTh2S[i,m,grid[[j]]],{j, 1, gridN},{i, m, l}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre integrator *******************)
EPMALegendreIntgTh2T[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgTh2T[i,k,grid[[j]]], {k, 0, m},{j, 1, gridN},{i, k, l}], epmPrecision]];
EPMALegendreIntgTh2T[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgTh2T[i,m,grid[[j]]],{j, 1, gridN},{i, m, l}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre integrator *******************)
EPMALegendreIntgPh2S[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgPh2S[i,k,grid[[j]]], {k, 0, m},{j, 1, gridN},{i, k, l}], epmPrecision]];
EPMALegendreIntgPh2S[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgPh2S[i,m,grid[[j]]],{j, 1, gridN},{i, m, l}], epmPrecision]];


(*****************  Compute EPM precision Associated Legendre integrator *******************)
EPMALegendreIntgPh2T[l_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgPh2T[i,k,grid[[j]]], {k, 0, m},{j, 1, gridN},{i, k, l}], epmPrecision]];
EPMALegendreIntgPh2T[l_,mmax_,m_] := Module[{grid =alGrid[l], gridN = alGridN[l], weights = alWeights[l]},N[Table[weights[[j]]alIntgPh2T[i,m,grid[[j]]],{j, 1, gridN},{i, m, l}], epmPrecision]];


End[]


EndPackage[]
