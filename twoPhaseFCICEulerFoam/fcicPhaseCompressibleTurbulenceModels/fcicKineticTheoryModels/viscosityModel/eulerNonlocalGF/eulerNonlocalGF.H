/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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
Class
    Foam::fcicKineticTheoryModels::viscosityModels::eulerNonlocalGF
Description
SourceFiles
    eulerNonlocalGFViscosity.C
\*---------------------------------------------------------------------------*/

#ifndef eulerNonlocalGF_H
#define eulerNonlocalGF_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fcicKineticTheoryModels
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class eulerNonlocalGF Declaration
\*---------------------------------------------------------------------------*/

class eulerNonlocalGF
:
    public viscosityModel
{
    // Private data

  dictionary eulerNonlocalGFCoeffs_;
  

  dimensionedScalar graind_;
  dimensionedScalar fcs_;
  dimensionedScalar delfc_;
  dimensionedScalar mus_;
  dimensionedScalar in0_;
  dimensionedScalar nu0_;
  dimensionedScalar rho_;
  dimensionedScalar A_;
  dimensionedScalar t0inv_;
  scalar alphaBoundary;
  scalar nNormalCell;
  scalar pres_offset_;
  scalar alphaMax_;

  volScalarField href_nlgf;
  volScalarField g_nlgf;
  volScalarField gsrc_expl;
  volScalarField gdcoeff;
  volScalarField gsrc_impl;
  //  volScalarField srnz;
  //  volScalarField pres_;
  uniformDimensionedVectorField g_;
  //  volScalarField p_rgh;
  //  volScalarField sigNu_;
  volScalarField nu_;
  uniformDimensionedScalarField href_;
  

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("eulerNonlocalGF");


    // Constructors

        //- Construct from components
        eulerNonlocalGF
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
	    const volScalarField& alpha,
	    const volScalarField& rho,
	    const surfaceScalarField& alphaRhoPhi
        );


    //- Destructor
  virtual ~eulerNonlocalGF();
  //  {}


    // Member Functions

          //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }
  
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }
  
        //- Correct the laminar viscosity
          virtual void correct()
        {
            nu_ = calcNu();
        }

  //- Read transportProperties dictionary
        virtual bool read();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace fcicKineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
