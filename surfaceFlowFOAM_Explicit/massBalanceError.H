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
    Foam::massBalanceError[surfaceFlowFOAM_Explicit]

Group
    surfaceFlowFOAM_Explicit

Description
    Calculates the cumulative Mass Balance Error after each time step.
\*-----------------------------------------------------------------------*/

scalar MB1 = 0.0;
scalar MB2_temp = 0.0;
scalar MBE = 0.0;
		
//Calculation of Total additional mass in the domain (MB1)
MB1 = gSum(((h - h0)*mesh.V())());
		
//Calculation of Total net flux into the domain (MB2)
forAll (mesh.boundaryMesh(), i)
{
	forAll(qO.boundaryField()[i], j)
	{
		MB2_temp = MB2_temp + (-qO.boundaryField()[i][j]*mesh.magSf().boundaryField()[i][j]);
	}
}
	
reduce(MB2_temp, sumOp<scalar>());
	
MB2_temp = MB2_temp + gSum(((P - I)*mesh.V())());
	
MB2 = MB2 + MB2_temp*runTime.deltaT().value();
			
MBE = mag(MB2 - MB1);
	
Info << "\nCummulative Mass Balance Error = " << MBE << endl;


