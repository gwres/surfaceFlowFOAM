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
    Foam::readModelParameters[surfaceFlowFOAM_DW]

Group
    surfaceFlowFOAM_DW

Description
    Reads the input parameters for surface flow
\*-----------------------------------------------------------------------*/

word bedSlopeMode(modelDict.lookup("bedSlopeMode"));
word Pmode(modelDict.lookup("precipitationMode"));
const scalar delta_htol = modelDict.lookupOrDefault<scalar>("h_tolerance", 1.0);
const scalar lc = modelDict.lookupOrDefault<scalar>("leastCount", 1.0);
const scalar Tr = modelDict.lookupOrDefault<scalar>("precipitationDuration", 1.0);
const scalar Tmax = modelDict.lookupOrDefault<scalar>("precipitationPeakTime", 1.0);
const scalar Pmax = modelDict.lookupOrDefault<scalar>("precipitationPeak", 1.0);

