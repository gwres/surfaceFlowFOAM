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
    Foam::specifyOverlandPrecipitation[surfaceFlowFOAM]

Group
    surfaceFlowFOAM

Description
    Specifies the precipitation over the overland surface. (Problem Specific)
\*-----------------------------------------------------------------------*/
if (Pmode.match("steady"))
{
	if (runTime.value() > Tr)
	{
		P = 0.0*P;
		P.correctBoundaryConditions();
	}
}
else if (Pmode.match("unsteady"))
{
	if (runTime.value() <= Tmax)
	{
		forAll(P, i)
		{
			P[i] = Pmax*runTime.value()/Tmax;
		}	
	}
	else if ((runTime.value() > Tmax) && (runTime.value() <= Tr))
	{
		forAll(P, i)
		{
			P[i] = Pmax*(Tr - runTime.value())/(Tr - Tmax);
		}
	}
	else if (runTime.value() > Tr)
	{
		P = 0.0*P;
	}
	P.correctBoundaryConditions();
}
				



