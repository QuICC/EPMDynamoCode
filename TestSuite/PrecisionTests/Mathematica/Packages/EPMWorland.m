(* ::Package:: *)

BeginPackage["EPMWorland`"]


EPMWorlandGridSize::usage = 
"EPMWorlandGridSize[N,L] Get the size of the radial grid.";


EPMWorlandGrid::usage =
"EPMWorlandGrid[N,L] give the radial worland grid points.";


EPMWorlandWeights::usage =
"EPMWorlandWeights[N,L] give the radial worland quadrature weights.";


EPMWorlandP::usage =
"EPMWorlandP[N,L] get the worland polynomials on the grid.";


EPMWorlandP::usage =
"EPMWorlandP[N,L, l] get the worland polynomials on the grid for a given l.";


EPMWorlandIntg::usage = 
"EPMWorlandIntg[N,L] Get the intg integrator . EPMWorlandIntg[N,L, l] Get the intg integrator for a given l.";


EPMWorlandIntgQ2Pol::usage = 
"EPMWorlandIntgQ2Pol[N,L] Get the intg integrator . EPMWorlandIntgQ2Pol[N,L, l] Get the intg integrator for a given l.";


EPMWorlandIntgT2Tor::usage = 
"EPMWorlandIntgT2Tor[N,L] Get the intg integrator . EPMWorlandIntgT2Tor[N,L, l] Get the intg integrator for a given l.";


EPMWorlandIntgQ2CurlCurlProj::usage = 
"EPMWorlandIntgQ2CurlCurlProj[N,L] Get the intg integrator . EPMWorlandIntgQ2CurlCurlProj[N,L, l] Get the intg integrator for a given l.";


EPMWorlandIntgS2CurlCurlProj::usage = 
"EPMWorlandIntgS2CurlCurlProj[N,L] Get the intg integrator . EPMWorlandIntgS2CurlCurlProj[N,L, l] Get the intg integrator for a given l.";


EPMWorlandIntgT2CurlProj::usage = 
"EPMWorlandIntgT2CurlProj[N,L] Get the intg integrator . EPMWorlandIntgT2CurlProj[N,L, l] Get the intg integrator for a given l.";


EPMWorlandProj::usage = 
"EPMWorlandProj[N,L] Get the proj projector . EPMWorlandProj[N,L, l] Get the proj projector for a given l.";


EPMWorlandDProj::usage = 
"EPMWorlandDProj[N,L] Get the dProj projector . EPMWorlandDProj[N,L, l] Get the dProj prjector for a given l.";


EPMWorlandD2Proj::usage = 
"EPMWorlandD2Proj[N,L] Get the d2Proj projector . EPMWorlandD2Proj[N,L, l] Get the d2Proj projector for a given l.";


EPMWorlandProj2GradTP::usage = 
"EPMWorlandProj2GradTP[N,L] Get the proj2GradTP projector . EPMWorlandProj2GradTP[N,L, l] Get the proj2GradTP projector for a given l.";


EPMWorlandProjPol2Q::usage = 
"EPMWorlandProjPol2Q[N,L] Get the intg integrator . EPMWorlandProjPol2Q[N,L, l] Get the intg integrator for a given l.";


EPMWorlandProjPol2S::usage = 
"EPMWorlandProjPol2S[N,L] Get the intg integrator . EPMWorlandProjPol2S[N,L, l] Get the intg integrator for a given l.";


EPMWorlandProjTor2T::usage = 
"EPMWorlandProjTor2T[N,L] Get the intg integrator . EPMWorlandProjTor2T[N,L, l] Get the intg integrator for a given l.";


EPMWorlandProjTor2CurlQ::usage = 
"EPMWorlandProjTor2CurlQ[N,L] Get the intg integrator . EPMWorlandProjTor2CurlQ[N,L, l] Get the intg integrator for a given l.";


EPMWorlandProjTor2CurlS::usage = 
"EPMWorlandProjTor2CurlS[N,L] Get the intg integrator . EPMWorlandProjTor2CurlS[N,L, l] Get the intg integrator for a given l.";


EPMWorlandProjPol2CurlT::usage = 
"EPMWorlandProjPol2CurlT[N,L] Get the intg integrator . EPMWorlandProjPol2CurlT[N,L, l] Get the intg integrator for a given l.";


Begin["`Private`"]


(*****************  Set precisions *******************)
highPrecision=80;
epmPrecision=16;


(*****************  Compute grid size *******************)
wGridN[n_, l_] = Ceiling[3/2 n] + Ceiling[3/4 l] + 1;
EPMWorlandGridSize[n_,l_] = n + Floor[l/2] + 1;


(*****************  Setup general normalisation function *******************)
jacobiNorm[\[Alpha]_,\[Beta]_,n_]:=If[n==0 && \[Beta]==-1/2,N[\[Pi]/2, highPrecision],N[1/2 1/(2n +\[Alpha]+\[Beta]+1) (Gamma[n+\[Alpha]+1]Gamma[n+\[Beta]+1])/(n!Gamma[n+\[Alpha]+\[Beta]+1]), highPrecision]];


(*****************  Setup Worland normalisation *******************)
wNorm[l_,n_] := jacobiNorm[-1/2, l - 1/2, n]


(*****************  Create private polynomials and derivatives computations *******************)
wPoly[n_Integer, l_Integer, r_Real]:= r^l JacobiP[n,-1/2, l-1/2,2r^2-1];
wDPoly[n_Integer, l_Integer, r_Real]:=Evaluate[D[r^l JacobiP[n,-1/2, l-1/2,2r^2-1], r]];
wDDPoly[n_Integer, l_Integer, r_Real]:= Evaluate[D[r^l JacobiP[n,-1/2, l-1/2,2r^2-1],{r,2}]];
wProj[n_Integer, l_Integer, r_Real]:= r^l JacobiP[n,-1/2, l-1/2,2r^2-1];
wProj2GradTP[n_Integer, l_Integer, r_Real]:= r^(l-1) JacobiP[n,-1/2, l-1/2,2r^2-1];
wProjPol2Q[n_Integer, l_Integer, r_Real]:= l (l+1) r^(l-1) JacobiP[n,-1/2, l-1/2,2r^2-1];
wProjPol2S[n_Integer, l_Integer, r_Real]:=Evaluate[Sqrt[l (l+1)] 1/r D[r^(l+1)JacobiP[n,-1/2, l-1/2,2r^2-1],r]];
wProjTor2T[n_Integer, l_Integer, r_Real]:= -Sqrt[l (l+1)]r^l JacobiP[n,-1/2, l-1/2,2r^2-1];
wProjTor2CurlQ[n_Integer, l_Integer, r_Real]:= l (l+1) r^(l-1) JacobiP[n,-1/2, l-1/2,2r^2-1];
wProjTor2CurlS[n_Integer, l_Integer, r_Real]:=Evaluate[Sqrt[l (l+1)] 1/r D[r^(l+1)JacobiP[n,-1/2, l-1/2,2r^2-1],r]];
wProjPol2CurlT[n_Integer, l_Integer, r_Real]:= Sqrt[l (l+1)](Evaluate[Simplify[D[x^l JacobiP[n,-1/2, l-1/2,2 x^2-1],{x,2}] + 2/x D[x^l JacobiP[n,-1/2, l-1/2,2x^2-1], x] - l (l+1)x^(l-2) JacobiP[n,-1/2, l-1/2,2x^2-1]]])/.{x->r};
wIntg[n_Integer, l_Integer, r_Real]:= r^l JacobiP[n,-1/2, l-1/2,2r^2-1];
wIntgQ2Pol[n_Integer, l_Integer, r_Real]:= If[l != 0,1/(l (l+1)),0]r^(l+1) JacobiP[n,-1/2, l-1/2,2r^2-1];
wIntgT2Tor[n_Integer, l_Integer, r_Real]:= If[l != 0,-1/Sqrt[l (l+1)],0]r^l JacobiP[n,-1/2, l-1/2,2r^2-1];
wIntgQ2CurlCurlProj[n_Integer, l_Integer, r_Real]:= r^(l-1) JacobiP[n,-1/2, l-1/2,2r^2-1];
wIntgS2CurlCurlProj[n_Integer, l_Integer, r_Real]:= If[l != 0,-1/Sqrt[l (l+1)],0](r^(l-2) JacobiP[n,-1/2, l-1/2,2r^2-1] + Evaluate[D[r^(l-1) JacobiP[n,-1/2, l-1/2,2r^2-1], r]]);
wIntgT2CurlProj[n_Integer, l_Integer, r_Real]:= If[l != 0,-1/Sqrt[l (l+1)],0]r^l JacobiP[n,-1/2, l-1/2,2r^2-1];


(*****************  Create high precision grid values *******************)
wGrid[n_,l_] := Sort[N[Sqrt[(Table[Cos[\[Pi]/2 (2i-1)/wGridN[n,l] ],{i,1,wGridN[n,l]}]+1)/2], highPrecision]];


(*****************  Create high precision weights *******************)
wWeights[n_,l_] := Table[N[\[Pi]/(2wGridN[n,l]),highPrecision],{i,wGridN[n,l]}];


(*****************  Compute EPM precision grid *******************)
EPMWorlandGrid[n_,l_] := N[wGrid[n,l], epmPrecision];


(*****************  Compute EPM precision weights *******************)
EPMWorlandWeights[n_,l_] := N[wWeights[n,l], epmPrecision];


(*****************  Compute EPM precision Worland polynomials *******************)
EPMWorlandP[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wPoly[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandP[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wPoly[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision first derivative *******************)
EPMWorlandDP[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wDPoly[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandDP[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wDPoly[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision second derivative *******************)
EPMWorlandDDP[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wDDPoly[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandDDP[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wDDPoly[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProj[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wPoly[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProj[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wPoly[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandDProj[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wDPoly[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandDProj[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wDPoly[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandD2Proj[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wDDPoly[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandD2Proj[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wDDPoly[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProj2GradTP[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProj2GradTP[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProj2GradTP[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProj2GradTP[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProjPol2Q[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProjPol2Q[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProjPol2Q[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProjPol2Q[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProjPol2S[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProjPol2S[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProjPol2S[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProjPol2S[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProjTor2T[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProjTor2T[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProjTor2T[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProjTor2T[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProjTor2CurlQ[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProjTor2CurlQ[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProjTor2CurlQ[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProjTor2CurlQ[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProjTor2CurlS[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProjTor2CurlS[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProjTor2CurlS[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProjTor2CurlS[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland proj projector *******************)
EPMWorlandProjPol2CurlT[n_,l_] := Module[{grid =wGrid[n,l], gridN = wGridN[n,l]},N[Table[wProjPol2CurlT[i,k,grid[[j]]], {k, 0, l},{i, 0, n},{j, 1, gridN}], epmPrecision]];
EPMWorlandProjPol2CurlT[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], gridN = wGridN[n,lmax]},N[Table[wProjPol2CurlT[i,l,grid[[j]]],{i, 0, n},{j, 1, gridN}], epmPrecision]];


(*****************  Compute EPM precision Worland intg integrator *******************)
EPMWorlandIntg[n_,l_] := Module[{grid =wGrid[n,l], weights=wWeights[n,l], gridN = wGridN[n,l]},N[Table[weights[[j]]/wNorm[k,i] wIntg[i,k,grid[[j]]], {k, 0, l},{j, 1, gridN},{i, 0, n}], epmPrecision]];
EPMWorlandIntg[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], weights=wWeights[n,lmax], gridN = wGridN[n,lmax]},N[Table[weights[[j]]/wNorm[l,i] wIntg[i,l,grid[[j]]],{j, 1, gridN},{i, 0, n}], epmPrecision]];


(*****************  Compute EPM precision Worland intg integrator *******************)
EPMWorlandIntgQ2Pol[n_,l_] := Module[{grid =wGrid[n,l], weights=wWeights[n,l], gridN = wGridN[n,l]},N[Table[weights[[j]]/wNorm[k,i] wIntgQ2Pol[i,k,grid[[j]]], {k, 0, l},{j, 1, gridN},{i, 0, n}], epmPrecision]];
EPMWorlandIntgQ2Pol[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], weights=wWeights[n,lmax], gridN = wGridN[n,lmax]},N[Table[weights[[j]]/wNorm[l,i] wIntgQ2Pol[i,l,grid[[j]]],{j, 1, gridN},{i, 0, n}], epmPrecision]];


(*****************  Compute EPM precision Worland intg integrator *******************)
EPMWorlandIntgT2Tor[n_,l_] := Module[{grid =wGrid[n,l], weights=wWeights[n,l], gridN = wGridN[n,l]},N[Table[weights[[j]]/wNorm[k,i] wIntgT2Tor[i,k,grid[[j]]], {k, 0, l},{j, 1, gridN},{i, 0, n}], epmPrecision]];
EPMWorlandIntgT2Tor[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], weights=wWeights[n,lmax], gridN = wGridN[n,lmax]},N[Table[weights[[j]]/wNorm[l,i] wIntgT2Tor[i,l,grid[[j]]],{j, 1, gridN},{i, 0, n}], epmPrecision]];


(*****************  Compute EPM precision Worland intg integrator *******************)
EPMWorlandIntgQ2CurlCurlProj[n_,l_] := Module[{grid =wGrid[n,l], weights=wWeights[n,l], gridN = wGridN[n,l]},N[Table[weights[[j]]/wNorm[k,i]wIntgQ2CurlCurlProj[i,k,grid[[j]]], {k, 0, l},{j, 1, gridN},{i, 0, n}], epmPrecision]];
EPMWorlandIntgQ2CurlCurlProj[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], weights=wWeights[n,lmax], gridN = wGridN[n,lmax]},N[Table[weights[[j]]/wNorm[l,i]wIntgQ2CurlCurlProj[i,l,grid[[j]]],{j, 1, gridN},{i, 0, n}], epmPrecision]];


(*****************  Compute EPM precision Worland intg integrator *******************)
EPMWorlandIntgS2CurlCurlProj[n_,l_] := Module[{grid =wGrid[n,l], weights=wWeights[n,l], gridN = wGridN[n,l]},N[Table[weights[[j]]/wNorm[k,i] wIntgS2CurlCurlProj[i,k,grid[[j]]], {k, 0, l},{j, 1, gridN},{i, 0, n}], epmPrecision]];
EPMWorlandIntgS2CurlCurlProj[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], weights=wWeights[n,lmax], gridN = wGridN[n,lmax]},N[Table[weights[[j]]/wNorm[l,i] wS2CurlCurlProj[i,l,grid[[j]]],{j, 1, gridN},{i, 0, n}], epmPrecision]];


(*****************  Compute EPM precision Worland intg integrator *******************)
EPMWorlandIntgT2CurlProj[n_,l_] := Module[{grid =wGrid[n,l], weights=wWeights[n,l], gridN = wGridN[n,l]},N[Table[weights[[j]]/wNorm[k,i] wIntgT2CurlProj[i,k,grid[[j]]], {k, 0, l},{j, 1, gridN},{i, 0, n}], epmPrecision]];
EPMWorlandIntgT2CurlProj[n_,lmax_,l_] := Module[{grid =wGrid[n,lmax], weights=wWeights[n,lmax], gridN = wGridN[n,lmax]},N[Table[weights[[j]]/wNorm[l,i] wIntgT2CurlProj[i,l,grid[[j]]],{j, 1, gridN},{i, 0, n}], epmPrecision]];


End[]


EndPackage[]
