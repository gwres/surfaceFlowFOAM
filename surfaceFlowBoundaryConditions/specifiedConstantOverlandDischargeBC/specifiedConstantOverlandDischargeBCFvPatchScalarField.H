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
    Foam::specifiedConstantOverlandDischargeBCFvPatchScalarField

Group
    specifiedConstantOverlandDischargeBC/surfaceFlowBoundaryConditions

Description
    This boundary condition calculates the normal gradient of flow depth 
    due to time constant specified discharge at that boundary.
\*-----------------------------------------------------------------------*/
#ifndef specifiedConstantOverlandDischargeBCFvPatchScalarField_H
#define specifiedConstantOverlandDischargeBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{
class specifiedConstantOverlandDischargeBCFvPatchScalarField
:
	public fixedGradientFvPatchScalarField
	{
		//Private Data
		vectorField qC_; //Time Constant Overland Discharge
		
		public:

			//Runtime type information
			TypeName("specifiedConstantOverlandDischargeBC");
			
			//Constructors
				
				//Construct from patch and internal field
				specifiedConstantOverlandDischargeBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);
				
				//Construct from patch,internal field and dictionary
				specifiedConstantOverlandDischargeBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);
				
				//Construct by mapping given specifiedConstantOverlandDischargeBCFvPatchScalarField onto a new patch
				specifiedConstantOverlandDischargeBCFvPatchScalarField
				(
					const specifiedConstantOverlandDischargeBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);
				
				//- Construct as copy
				specifiedConstantOverlandDischargeBCFvPatchScalarField
				(
					const specifiedConstantOverlandDischargeBCFvPatchScalarField&
				);
				
				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new specifiedConstantOverlandDischargeBCFvPatchScalarField(*this)
					);
				}
				
				//Construct as copy setting internal field reference
				specifiedConstantOverlandDischargeBCFvPatchScalarField
				(
					const specifiedConstantOverlandDischargeBCFvPatchScalarField&,
					const DimensionedField<scalar, volMesh>&
				);
				
				//Construct and return a clone setting internal field reference
				virtual tmp<fvPatchScalarField> clone
				(
					const DimensionedField<scalar, volMesh>& iF
				) const
				{
					return tmp<fvPatchScalarField>
					(
						new specifiedConstantOverlandDischargeBCFvPatchScalarField(*this, iF)
					);
				}
				
			//Member Functions
			
				//Update Coefficients
				virtual void updateCoeffs();
				
				//Write
				virtual void write(Ostream&) const;
				
				//Return Time Constant Overland Discharge
				virtual const vectorField& qC() const
				{
					return qC_;
				}

				virtual vectorField& qC()
				{
					return qC_;
				}
				
	};
	
} // End namespace Foam

#endif
