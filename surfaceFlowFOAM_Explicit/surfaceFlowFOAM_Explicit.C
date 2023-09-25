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
    Foam::surfaceFlowFOAM_Explicit

Group
    surfaceFlowFOAM_Explicit

Description
    Modeling of overland flow [Explicit formulation]
\*-----------------------------------------------------------------------*/
#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    std::clock_t startT= std::clock(); 									//Start Time
    
    argList::addNote
    (
        "Surface Flow Model"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readTimeControls.H"
    
    //------------------------------------------------------------------------------------------------------------------------------//
    #include "readModelParameters.H" 									//Reading the input parameters for Shallow Water Equation
        
    std::ofstream file2;
	file2.open("outletDischarge.csv");
	
	std::ofstream file3;
	file3.open("massBalanceError.csv");
	
	//------------------------------------------------------------------------------------------------------------------------------//
	label N_ele = mesh.cells().size(); 									//Number of Elements
	
	scalar MB2 = 0.0;			       									//Initialization of total net flux variable
    
    volScalarField h0 = h;    	 	   									//Initial flow depth
	
    //------------------------------------------------------------------------------------------------------------------------------//
    //Initialization of parameters for DW Model
	#include "initializeSurfaceParameters.H"
			
	//--------------------------------------------------------------------------------------------------------------------------//
	//Time-Loop
    while (runTime.loop())
    {
		#include "readTimeControls.H"
		
		Info << "\ndeltaT = " <<  runTime.deltaT().value() << endl;
        
        Info << "\nTime = " << runTime.timeName() << endl;	
        
        #include "specifyPrecipitation.H"
		
		#include "zeroInertiaModel.H"
        
        //--------------------------------------------------------------------------------------------------------------------------//
        #include "massBalanceError.H" 									//Calculation of Mass Balance Error (MBE)
		
		file3 << runTime.timeName() << "," << MB1 << "," << MB2 << "," << MBE << nl;	//Writing the MBE
		
		//--------------------------------------------------------------------------------------------------------------------------//
		#include "outletDischarge.H"									//Calculates the outlet discharge
		
		runTime.write();
		
		Info << "-------------------------------------------------------" << endl;
		
	} 																	//End of Time-Loop
	
	//------------------------------------------------------------------------------------------------------------------------------//
	file2.close();
	file3.close();
	
	std::clock_t endT= std::clock(); 									//End Time
	
	//Writing the Simulation Time
	scalar simTime;
	simTime = (endT - startT)/ (double) CLOCKS_PER_SEC;
	Info<< "\nCPU simulation time = " << simTime << nl << endl;
    
    Info<< "End\n" << endl;

    return 0;
}
