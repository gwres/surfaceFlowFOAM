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
#include "specifiedConstantOverlandDischargeBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
	//Constructors
	
	specifiedConstantOverlandDischargeBCFvPatchScalarField::specifiedConstantOverlandDischargeBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_(0.0)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedConstantOverlandDischargeBCFvPatchScalarField::specifiedConstantOverlandDischargeBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		qC_("qC", dict, p.size())
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
	specifiedConstantOverlandDischargeBCFvPatchScalarField::specifiedConstantOverlandDischargeBCFvPatchScalarField
	(
		const specifiedConstantOverlandDischargeBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
		qC_(ptf.qC_, mapper)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedConstantOverlandDischargeBCFvPatchScalarField::specifiedConstantOverlandDischargeBCFvPatchScalarField
	(
		const specifiedConstantOverlandDischargeBCFvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf),
		qC_(ptf.qC_)
		{}
	/*-------------------------------------------------------------------------------------*/
	specifiedConstantOverlandDischargeBCFvPatchScalarField::specifiedConstantOverlandDischargeBCFvPatchScalarField
	(
		const specifiedConstantOverlandDischargeBCFvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF),
		qC_(ptf.qC_)
		{}
	
	//Member Functions
	
	void specifiedConstantOverlandDischargeBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		const fvPatchField<tensor>& alphaO_B = patch().lookupPatchField<volTensorField, tensor>("alphaO");
       	const fvPatchField<vector>& S0_B = patch().lookupPatchField<volVectorField, vector>("S0");
                
        gradient() = -((qC_ & inv(alphaO_B)) & patch().nf()) - (S0_B & patch().nf());
         
        fixedGradientFvPatchScalarField::updateCoeffs();
	}
	
	//Write
	void specifiedConstantOverlandDischargeBCFvPatchScalarField::write(Ostream& os) const
	{
		fixedGradientFvPatchScalarField::write(os);
		qC_.writeEntry("qC", os);
		writeEntry("value", os);
	}
	
	makePatchTypeField(fvPatchScalarField, specifiedConstantOverlandDischargeBCFvPatchScalarField);
			
} //End namespace Foam
