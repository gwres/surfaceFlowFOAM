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
    Foam::outletDischarge[surfaceFlowFOAM_Explicit]

Group
    surfaceFlowFOAM_Explicit

Description
    Calculates the discharge through the outlet boundary
\*-----------------------------------------------------------------------*/
scalar qO_out = 0.0;

label outletID = mesh.boundaryMesh().findPatchID("outlet");

forAll(mesh.boundaryMesh(), patchID)
{
	if ((patchID == outletID)) 
	{
		forAll(qO.boundaryField()[patchID], faceID)
		{
			qO_out = qO_out + qO.boundaryField()[patchID][faceID]*mesh.magSf().boundaryField()[patchID][faceID];
		}
	}
}

reduce(qO_out, sumOp<scalar>());

file2 << runTime.timeName() << "," << qO_out << nl;
