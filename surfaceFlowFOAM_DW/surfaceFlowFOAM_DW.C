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
    Foam::surfaceFlowFOAM_DW

Group
    surfaceFlowFOAM_DW

Description
    Diffusive Wave Modeling of overland flow 
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
        
    #include "readPicardParameters.H" 									//Reading the input parameters for Picard Iteration Method
    
    //Opening files for writing the Number of Picard Iterations, Outlet Discharge, and  MBE for every time-step in a CSV file
    std::ofstream file1;
	file1.open("PicardIterations.csv");
	
	std::ofstream file2;
	file2.open("outletDischarge.csv");
	
	std::ofstream file3;
	file3.open("massBalanceError.csv");
	
	//------------------------------------------------------------------------------------------------------------------------------//
	label N_ele = mesh.cells().size(); 									//Number of Elements
	
	label nMARK = 0;           	       									//Indicator for Convergence Failure
	
	label TrMARK = 0;           	       								//Indicator for end of precipitation 
    	
	scalar err = 0.0;    	    	   									//Initialization of Picard iteration error
        
    label sc = 0;				       									//Initialization of the Stabilization counter
        
    label NP = NPl;           	       									//Initialization of Picard Iteration number
    
    scalar MB2 = 0.0;			       									//Initialization of total net flux variable
    
    scalar Pic_tol = delta_htol;
    
	scalar oldTime = 0.0;	  	       									//Initialization of previous time variable
	label oldTimeIndex = 0;	  	       									//Initialization of previous time-index variable
	
	volScalarField h0 = h;    	 	   									//Initial flow depth
	
    //------------------------------------------------------------------------------------------------------------------------------//
    //Initialization of parameters for DW Model
	#include "initializeSurfaceParameters.H"
			
	//--------------------------------------------------------------------------------------------------------------------------//
	//Time-Loop
    while (runTime.loop())
    {
		#include "readTimeControls.H"
		
		//Reset DeltaT to minDeltaT after precipitation stops
		if ((runTime.value() > Tr) && (TrMARK == 0))
		{
			runTime.setDeltaT
			(
				minDeltaT
			);
			TrMARK = 1;
			sc = 0;
		}
		
		Info << "\ndeltaT = " <<  runTime.deltaT().value() << endl;
        
        Info << "\nTime = " << runTime.timeName() << endl;	
        
        #include "specifyPrecipitation.H"
        
        //--------------------------------------------------------------------------------------------------------------------------//
        //Picard Iteration Loop
		for(label count = 1; count < NPmax + 1; count++)
		{
			#include "waveModel.H"
			
			//----------------------------------------------------------------------------------------------------------------------//	
			//Calculating the Convergence Criterion
			err = Foam::pow((gSumMag(Foam::pow((h - hM), 2.0)())/N_ele), 0.5) ;
			
			//----------------------------------------------------------------------------------------------------------------------//
			//Termination criterion for Picard Iteration Loop
			reduce(count, maxOp<label>());
			
			if ((err < Pic_tol) && (count >= 2))
			{
				NP = count;
				Info << "\nError = " << err << endl;
				Info << "\nNo.of Picard Iterations = " << NP << endl;
				break;
			}
			else if (count == NPmax)
			{
				NP = count;
				Info << "\nConvergence failed!" << endl;
				Info << "\nError = " << err << endl;
				Info << "\nNo.of Picard Iterations = " << NP << endl;
				
				runTime.setTime											//Set Runtime to previous time-level
				(
					oldTime, oldTimeIndex
				);
			 
				h = hN;													//Resetting 'h' and 'hN' field values
				hN = hNm1;
				
				nMARK = 1;				
			}
			else
			{
				NP = count;
			}
	
			//----------------------------------------------------------------------------------------------------------------------//
			hM = h;
			
		} 																//End of Picard Iteration Loop
		
        //--------------------------------------------------------------------------------------------------------------------------//
        #include "massBalanceError.H" 									//Calculation of Mass Balance Error (MBE)
		
		if (nMARK != 1)
		{			
			file1 << runTime.timeName() << "," << NP << nl;				//Writing the Number of Picard Iterations
			file3 << runTime.timeName() << "," << MB1 << "," << MB2 << "," << MBE << nl;	//Writing the MBE
		}
		
		//--------------------------------------------------------------------------------------------------------------------------//
		
		
		#include "outletDischarge.H"									//Calculates the outlet discharge
		
		runTime.write();
		
		oldTime = runTime.value();
		oldTimeIndex = runTime.timeIndex();
		
		hNm1 = hN;
		hN = h;
		
		nMARK = 0;
		
		#include "setDeltaT.H"											//Adjust time-step
		
		Info << "-------------------------------------------------------" << endl;
		
	} 																	//End of Time-Loop
	
	//------------------------------------------------------------------------------------------------------------------------------//
	file1.close();
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
