/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
/*-----------------------------------------------------------------------*\
Class
    Foam::waveModel[surfaceFlowFOAM]

Group
    surfaceFlowFOAM

Description
    Solves the Zero-inertia Model for overland flow.
\*-----------------------------------------------------------------------*/
volVectorField gTerm = alphaO & S0;

h.oldTime() = hN;

fvScalarMatrix hEqn
(
	fvm::ddt(h)
	==
	fvm::laplacian(alphaO, h)
	+
	fvc::div(gTerm)
	+
	P
	-
	I
);
				
hEqn.relax();
hEqn.solve();

h = max(lc*hU, h); 														//Depth-positivity preservation

volScalarField coeff1 = Foam::pow(h, 5.0/3.0)/Foam::pow(hU, 5.0/3.0);
volScalarField coeff2 = Foam::pow(mag(S0), 0.5) ;
coeff2 = max(coeff2, lc);
			
alphaO = (inv(nRGH)*(coeff1/coeff2)) & qUO;
alphaO.correctBoundaryConditions();

U = -alphaO & (fvc::grad(h) + S0);										
qO = fvc::interpolate(U) & (mesh.Sf()/mesh.magSf());					//Updating the face-normal flux
	

	

