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
    Foam::imperviousSurfaceFlowBCFvPatchScalarField

Group
    imperviousSurfaceFlowBC/surfaceFlowBoundaryConditions

Description
    This boundary condition calculates the normal gradient of the flow depth 
    at an impervious boundary.
\*-----------------------------------------------------------------------*/
#include "imperviousSurfaceFlowBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
	//Constructors
	
	imperviousSurfaceFlowBCFvPatchScalarField::imperviousSurfaceFlowBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF)
		{}
	/*-------------------------------------------------------------------------------------*/
	imperviousSurfaceFlowBCFvPatchScalarField::imperviousSurfaceFlowBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF)
		{
			if (dict.found("gradient"))
			{
				gradient() = scalarField("gradient", dict, p.size());
				fixedGradientFvPatchScalarField::updateCoeffs();
				fixedGradientFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
				gradient() = 0.0;
			}
		}
	/*-------------------------------------------------------------------------------------*/	
	imperviousSurfaceFlowBCFvPatchScalarField::imperviousSurfaceFlowBCFvPatchScalarField
	(
		const imperviousSurfaceFlowBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
		{}
	/*-------------------------------------------------------------------------------------*/
	imperviousSurfaceFlowBCFvPatchScalarField::imperviousSurfaceFlowBCFvPatchScalarField
	(
		const imperviousSurfaceFlowBCFvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf)
		{}
	/*-------------------------------------------------------------------------------------*/
	imperviousSurfaceFlowBCFvPatchScalarField::imperviousSurfaceFlowBCFvPatchScalarField
	(
		const imperviousSurfaceFlowBCFvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF)
		{}
	
	//Member Functions
	
	void imperviousSurfaceFlowBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		const fvPatchField<vector>& S0_B = patch().lookupPatchField<volVectorField, vector>("S0");
                
        gradient() = -(S0_B & patch().nf());
         
        fixedGradientFvPatchScalarField::updateCoeffs();
	}
	
	//Write
	void imperviousSurfaceFlowBCFvPatchScalarField::write(Ostream& os) const
	{
		fixedGradientFvPatchScalarField::write(os);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, imperviousSurfaceFlowBCFvPatchScalarField);
			
} //End namespace Foam
