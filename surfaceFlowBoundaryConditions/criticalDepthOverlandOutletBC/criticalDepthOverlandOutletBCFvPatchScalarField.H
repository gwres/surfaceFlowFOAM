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

#ifndef criticalDepthOverlandOutletBCFvPatchScalarField_H
#define criticalDepthOverlandOutletBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{
class criticalDepthOverlandOutletBCFvPatchScalarField
:
    public fixedValueFvPatchScalarField
	{
    	public:

			//Runtime type information
			TypeName("criticalDepthOverlandOutletBC");

			//Constructors

				//Construct from patch and internal field
				criticalDepthOverlandOutletBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);

				//Construct from patch, internal field and dictionary
				criticalDepthOverlandOutletBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);

				//Construct by mapping given criticalDepthOverlandOutletBCFvPatchScalarField onto a new patch
				criticalDepthOverlandOutletBCFvPatchScalarField
				(
					const criticalDepthOverlandOutletBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);

				//Construct as copy
				criticalDepthOverlandOutletBCFvPatchScalarField
				(
					const criticalDepthOverlandOutletBCFvPatchScalarField&
				);

				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new criticalDepthOverlandOutletBCFvPatchScalarField(*this)
					);
				}

				//Construct as copy setting internal field reference
				criticalDepthOverlandOutletBCFvPatchScalarField
				(
					const criticalDepthOverlandOutletBCFvPatchScalarField&,
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
						new criticalDepthOverlandOutletBCFvPatchScalarField(*this, iF)
					);
				}

			//Member functions

				//Update Coefficients
				virtual void updateCoeffs();

				//Write
				virtual void write(Ostream&) const;
				
	};

} // End namespace Foam

#endif
