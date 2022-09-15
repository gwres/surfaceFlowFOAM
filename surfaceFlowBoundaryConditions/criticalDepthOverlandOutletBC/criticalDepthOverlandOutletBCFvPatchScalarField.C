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
/*---------------------------------------------------------------------------*\
Class
    Foam::criticalDepthOverlandOutletBCFvPatchScalarField
    
Group
	criticalDepthOverlandOutletBC/surfaceFlowBoundaryConditions

Description
    This boundary condition calculates the critical flow depth at the 
    overland outlet boundary.
\*---------------------------------------------------------------------------*/

#include "criticalDepthOverlandOutletBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
	//Constructors
	
	criticalDepthOverlandOutletBCFvPatchScalarField::criticalDepthOverlandOutletBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(p, iF)
		{}
	/*--------------------------------------------------------------------------*/
	criticalDepthOverlandOutletBCFvPatchScalarField::criticalDepthOverlandOutletBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedValueFvPatchScalarField(p, iF)
		{
			if (dict.found("value"))
			{
				fvPatchField<scalar>::operator = (scalarField("value", dict, p.size()));
				fixedValueFvPatchScalarField::updateCoeffs();
				fixedValueFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
			}
		}
	/*-------------------------------------------------------------------------------------*/	
	criticalDepthOverlandOutletBCFvPatchScalarField::criticalDepthOverlandOutletBCFvPatchScalarField
	(
		const criticalDepthOverlandOutletBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedValueFvPatchScalarField(ptf, p, iF, mapper)
		{}
	/*--------------------------------------------------------------------------*/
	criticalDepthOverlandOutletBCFvPatchScalarField::criticalDepthOverlandOutletBCFvPatchScalarField
	(
		const criticalDepthOverlandOutletBCFvPatchScalarField& ptf
	)
	:
		fixedValueFvPatchScalarField(ptf)
		{}
	/*-------------------------------------------------------------------------------------*/
	criticalDepthOverlandOutletBCFvPatchScalarField::criticalDepthOverlandOutletBCFvPatchScalarField
	(
		const criticalDepthOverlandOutletBCFvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedValueFvPatchScalarField(ptf, iF)
		{}

	//Member Functions

	void criticalDepthOverlandOutletBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}

        const fvsPatchField<scalar>& qO_B = patch().lookupPatchField<surfaceScalarField, scalar>("qO");
        
		scalarField hB = pow((pow(qO_B, 2.0)/9.807), (1.0/3.0));
			
		operator== (hB);
		
		fixedValueFvPatchScalarField::updateCoeffs();
	}

	void criticalDepthOverlandOutletBCFvPatchScalarField::write(Ostream& os) const
	{
		fvPatchScalarField::write(os);
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchScalarField, criticalDepthOverlandOutletBCFvPatchScalarField);

} //End namespace Foam

