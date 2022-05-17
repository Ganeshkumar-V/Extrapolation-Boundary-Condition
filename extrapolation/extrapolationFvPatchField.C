/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extrapolationFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::extrapolationFvPatchField<Type>::extrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gradient_(p.size(), Zero),
    owner_(iF.mesh().owner()),
    neighbour_(iF.mesh().neighbour()),
    celliF(iF.mesh().C().internalField()),
    inF_(iF.field()),
    bcells_(iF.mesh().boundary()[p.index()].faceCells()),
    icells_(p.size(), -1)
{}


template<class Type>
Foam::extrapolationFvPatchField<Type>::extrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    gradient_(p.size(), Zero),
    owner_(iF.mesh().owner()),
    neighbour_(iF.mesh().neighbour()),
    celliF(iF.mesh().C().internalField()),
    inF_(iF.field()),
    bcells_(iF.mesh().boundary()[p.index()].faceCells()),
    icells_(p.size(), -1)
{
  // Finding the cells next to the boundary cells
  forAll (bcells_, celli)
  {
    forAll (neighbour_, i)
    {
      if (celli == 0)
      {
        if (neighbour_[i] == bcells_[celli])
        {
     	    icells_[celli] = owner_[i];
        }
      }
      else
      {
        if ((neighbour_[i] == bcells_[celli]) && (owner_[i] != bcells_[celli-1]))
        {
           icells_[celli] = owner_[i];
        }
      }
    }
  }
  // Evaluate the boundary Condition
  evaluate();
}


template<class Type>
Foam::extrapolationFvPatchField<Type>::extrapolationFvPatchField
(
    const extrapolationFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper),
    owner_(ptf.owner_),
    neighbour_(ptf.neighbour_),
    celliF(ptf.celliF),
    inF_(ptf.inF_),
    bcells_(ptf.bcells_),
    icells_(ptf.icells_)
{
  Info << "Mapper function is called! " << endl;
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::extrapolationFvPatchField<Type>::extrapolationFvPatchField
(
    const extrapolationFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    gradient_(ptf.gradient_),
    owner_(ptf.owner_),
    neighbour_(ptf.neighbour_),
    celliF(ptf.celliF),
    inF_(ptf.inF_),
    bcells_(ptf.bcells_),
    icells_(ptf.icells_)
{}


template<class Type>
Foam::extrapolationFvPatchField<Type>::extrapolationFvPatchField
(
    const extrapolationFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_),
    owner_(ptf.owner_),
    neighbour_(ptf.neighbour_),
    celliF(ptf.celliF),
    inF_(ptf.inF_),
    bcells_(ptf.bcells_),
    icells_(ptf.icells_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::extrapolationFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    gradient_.autoMap(m);
}


template<class Type>
void Foam::extrapolationFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const extrapolationFvPatchField<Type>& fgptf =
        refCast<const extrapolationFvPatchField<Type>>(ptf);

    gradient_.rmap(fgptf.gradient_, addr);
}

template<class Type>
void Foam::extrapolationFvPatchField<Type>::findGradient()
{
  // Finding Gradient of the field at the boundary cell center
  forAll(bcells_, i)
  {
    gradient_[i] = (inF_[bcells_[i]] - inF_[icells_[i]])
            /Foam::mag(celliF[bcells_[i]] - celliF[icells_[i]]);
  }
}


template<class Type>
void Foam::extrapolationFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    // Finding the Gradients
    findGradient();

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::extrapolationFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
  return tmp<Field<Type>>::New(this->size(), pTraits<Type>::one);
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::extrapolationFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
  return gradient()/this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::extrapolationFvPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::extrapolationFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return gradient();
}


template<class Type>
void Foam::extrapolationFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
}


// ************************************************************************* //
