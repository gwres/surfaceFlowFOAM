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
    Foam::freeFlowOutletBCFvPatchScalarField
    
Group
	freeFlowOutletBC/surfaceFlowBoundaryConditions

Description
    This boundary condition calculates the free-flow depth at the 
    outlet boundary.
\*---------------------------------------------------------------------------*/

#ifndef freeFlowOutletBCFvPatchScalarField_H
#define freeFlowOutletBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{
class freeFlowOutletBCFvPatchScalarField
:
    public fixedValueFvPatchScalarField
	{
    	public:

			//Runtime type information
			TypeName("freeFlowOutletBC");

			//Constructors

				//Construct from patch and internal field
				freeFlowOutletBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);

				//Construct from patch, internal field and dictionary
				freeFlowOutletBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);

				//Construct by mapping given freeFlowOutletBCFvPatchScalarField onto a new patch
				freeFlowOutletBCFvPatchScalarField
				(
					const freeFlowOutletBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);

				//Construct as copy
				freeFlowOutletBCFvPatchScalarField
				(
					const freeFlowOutletBCFvPatchScalarField&
				);

				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new freeFlowOutletBCFvPatchScalarField(*this)
					);
				}

				//Construct as copy setting internal field reference
				freeFlowOutletBCFvPatchScalarField
				(
					const freeFlowOutletBCFvPatchScalarField&,
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
						new freeFlowOutletBCFvPatchScalarField(*this, iF)
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
