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
    Foam::createFields[surfaceFlowFOAM_Explicit]

Group
    surfaceFlowFOAM_Explicit

Description
    Modeling of overland flow [Explicit formulation]
/*------------------------------------------------------------------------------*/
Info<< "Reading SURFACE WATER DEPTH field h\n" << endl;

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading BED ELEVATION field zB\n" << endl;

volScalarField zB
(
    IOobject
    (
        "zB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading PRECIPITATION RATE field P\n" << endl;

volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading INFILTRATION RATE field I\n" << endl;

volScalarField I
(
    IOobject
    (
        "I",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading UNIT FLOW DEPTH field hU\n" << endl;

volScalarField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading BED-SLOPE field S0\n" << endl;

volVectorField S0
(
	IOobject
    (
		"S0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading CELL-CENTRED FLUX field U\n" << endl;

volVectorField U
(
	IOobject
    (
		"U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading FACE-NORMAL FLUX field qO\n" << endl;

surfaceScalarField qO
(
    IOobject
    (
        "qO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading ROUGHNESS TERM field nRGH\n" << endl;

volTensorField nRGH
(
    IOobject
    (
        "nRGH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading ALPHA COEFFICIENT alphaO\n" << endl;

volTensorField alphaO
(
    IOobject
    (
        "alphaO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
Info<< "Reading UNIT DISCHARGE PER UNIT WIDTH field qUO\n" << endl;

volTensorField qUO
(
    IOobject
    (
        "qUO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*------------------------------------------------------------------------------*/
// Reads user-defined data from picardDict and modelDict in system directory

IOdictionary modelDict
(
	IOobject
    (
		"modelDict",
		runTime.system(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    )
);
/*------------------------------------------------------------------------------*/


