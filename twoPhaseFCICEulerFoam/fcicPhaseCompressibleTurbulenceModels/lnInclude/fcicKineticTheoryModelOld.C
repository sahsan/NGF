/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fcicKineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::fcicKineticTheoryModel::fcicKineticTheoryModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<fcicPhaseCompressibleTurbulenceModel>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    phase_(phase),

    viscosityModel_
    (
        fcicKineticTheoryModels::viscosityModel::New
        (
	 type,
	 coeffDict_,
	 U,
	 phi,
	 alpha,
	 rho,
	 alphaRhoPhi
        )
    ),
    conductivityModel_
    (
        fcicKineticTheoryModels::conductivityModel::New
        (
            coeffDict_
        )
    ),
    radialModel_
    (
        fcicKineticTheoryModels::radialModel::New
        (
            coeffDict_
        )
    ),
    granularPressureModel_
    (
        fcicKineticTheoryModels::granularPressureModel::New
        (
            coeffDict_
        )
    ),
    fcicGranularPressureModel_
    (
        fcicKineticTheoryModels::fcicGranularPressureModel::New
	(
	 coeffDict_,
	 U,
	 alpha,
	 rho
        )
	),
    frictionalStressModel_
    (
        fcicKineticTheoryModels::frictionalStressModel::New
        (
            coeffDict_
        )
    ),

    equilibrium_(coeffDict_.lookup("equilibrium")),
    e_("e", dimless, coeffDict_),
    alphaMax_("alphaMax", dimless, coeffDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict_
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffDict_
    ),
//sigNu_(U.mesh().lookupObject<volScalarField>("sigNu")),
    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

  
    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
     )
/*
    sigmaDev_
    (
        IOobject
        (
            IOobject::groupName("sigmaDev", U.group()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        sigmaDev()
     )
*/
{
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::fcicKineticTheoryModel::~fcicKineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::fcicKineticTheoryModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<fcicPhaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        coeffDict().lookup("equilibrium") >> equilibrium_;
        e_.readIfPresent(coeffDict());
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());

        viscosityModel_->read();
        conductivityModel_->read();
        radialModel_->read();
      	fcicGranularPressureModel_->read();
        granularPressureModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::fcicKineticTheoryModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::fcicKineticTheoryModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::fcicKineticTheoryModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
	    - (nut_)*dev(twoSymm(fvc::grad(U_)))//*scalar(0.0)
	    //          - (lambda_*fvc::div(phi_))*symmTensor::I
        )
    );
}


/*
Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::fcicKineticTheoryModel::sigmaDev() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("sigmaDev", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho_*(2.0*nut_*symm(fvc::grad(U_)) - (2.0/3.0)*nut_*tr(symm(fvc::grad(U_)))*symmTensor::I)
        )
    );
}
*/


Foam::tmp<Foam::volScalarField>
Foam::RASModels::fcicKineticTheoryModel::pPrime() const
{
  const volScalarField& rho = phase_.rho();

    tmp<volScalarField> tpPrime
    (
     Theta_
     *granularPressureModel_->granularPressureCoeffPrime
     (
      alpha_,
      radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
      radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
      rho,
      e_
      )
      +
     frictionalStressModel_->frictionalPressurePrime
     (
      phase_,
      alphaMinFriction_,
      alphaMax_
      )
     /*  +
	 fcicGranularPressureModel_->granularPressureCoeffPrime()   */  
     );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }


    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::fcicKineticTheoryModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::fcicKineticTheoryModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (rho_*nut_)
	    *dev(twoSymm(fvc::grad(U_)))//*scalar(0.0)
	    //          - ((rho_*lambda_)*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::fcicKineticTheoryModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
     (- fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            (rho_*nut_)*dev2(T(fvc::grad(U)))
	    //          + ((rho_*lambda_)*fvc::div(phi_))
	    //           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
	 )
      )
    );
}

void Foam::RASModels::fcicKineticTheoryModel::correct()
{
    // Local references
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    const volVectorField& Uc_ =
        refCast<const twoPhaseSystem>(phase_.fluid()).otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    if (!equilibrium_)
    {
        // Particle viscosity (Table 3.2, p.47)

      nut_ = viscosityModel_->calcNu();//nu(alpha, Theta_, gs0_, rho, da, e_);
      //      sigNu_ = viscosityModel_->calcNu();

      //sigNu_.internalField() = nut_.internalField();
      //sigNu_.correctBoundaryConditions();

      volScalarField LunPre(
			    "LunPre",
			    granularPressureModel_->granularPressureCoeffPrime
			    (
			     alpha_,
			     radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
			     radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
			     rho,
			     e_
			     )
			    );
      
      if(U_.time().outputTime())
	{
	  LunPre.write();
	}

      volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
		lambda_ = scalar(0.0)*(4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

		//lambda_ =(4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;
	
        // Stress tensor, Definitions, Table 3.1, p. 43
        volSymmTensorField tau
        (
	 rho*(2.0*nut_*D)// - (2.0/3.0)*nut_*tr(D)*I)// + (lambda_ - (2.0/3.0)*nut_)*tr(D)*I)
        );


	const volScalarField justone
            (
             "justone",nut_/nut_
            );
        double volume=fvc::domainIntegrate(justone).value();
	Info << max(nut_) <<endl;
	Info <<	min(nut_) <<endl;
	Info << "average(nut_):"<<fvc::domainIntegrate(nut_).value()/volume<<"\n";
        // Dissipation (Eq. 3.24, p.50)
        volScalarField gammaCoeff
        (
            "gammaCoeff",
            12.0*(1.0 - sqr(e_))
           *max(sqr(alpha), residualAlpha_)
           *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
        );

        // Drag
        volScalarField beta
        (
            refCast<const twoPhaseSystem>(phase_.fluid()).drag(phase_).K()
        );

        // Eq. 3.25, p. 50 Js = J1 - J2
        volScalarField J1("J1", 3.0*beta);
        volScalarField J2
        (
            "J2",
            0.25*sqr(beta)*da*magSqr(U - Uc_)
           /(
               max(alpha, residualAlpha_)*rho
              *sqrtPi*(ThetaSqrt + ThetaSmallSqrt)
            )
        );

        // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
        volScalarField PsCoeff
        (
            granularPressureModel_->granularPressureCoeff
            (
                alpha,
                gs0_,
                rho,
                e_
            )
        );

        // 'thermal' conductivity (Table 3.3, p. 49)
        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);

        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20:
        //     Ps should be without grad
        //     the laplacian has the wrong sign
        fvScalarMatrix ThetaEqn
        (
            1.5*
            (
                fvm::ddt(alpha, rho, Theta_)
              + fvm::div(alphaRhoPhi, Theta_)
              - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
            )
          - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
         ==
          - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
          + (tau && gradU)
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
          + fvOptions(alpha, rho, Theta_)
        );

        ThetaEqn.relax();
        fvOptions.constrain(ThetaEqn);
        ThetaEqn.solve();
        fvOptions.correct(Theta_);
    }
    else
    {
        // Equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1("K1", 2.0*(1.0 + e_)*rho*gs0_);
        volScalarField K3
        (
            "K3",
            0.5*da*rho*
            (
                (sqrtPi/(3.0*(3.0 - e_)))
               *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha*gs0_)
               +1.6*alpha*gs0_*(1.0 + e_)/sqrtPi
            )
        );

        volScalarField K2
        (
            "K2",
            4.0*da*rho*(1.0 + e_)*alpha*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0
        );

        volScalarField K4("K4", 12.0*(1.0 - sqr(e_))*rho*gs0_/(da*sqrtPi));

        volScalarField trD
        (
            "trD",
            alpha/(alpha + residualAlpha_)
           *fvc::div(phi_)
        );
        volScalarField tr2D("tr2D", sqr(trD));
        volScalarField trD2("trD2", tr(D & D));

        volScalarField t1("t1", K1*alpha + rho);
        volScalarField l1("l1", -t1*trD);
        volScalarField l2("l2", sqr(t1)*tr2D);
        volScalarField l3
        (
            "l3",
            4.0
           *K4
           *alpha
           *(2.0*K3*trD2 + K2*tr2D)
        );

        Theta_ = sqr
        (
            (l1 + sqrt(l2 + l3))
           /(2.0*max(alpha, residualAlpha_)*K4)
        );

        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);
    }

    Theta_.max(0);
    Theta_.min(100);

    {
        // particle viscosity (Table 3.2, p.47)
      /*      volScalarField pres_regK
      	(
      	 "pres_regK",
      	 fcicGranularPressureModel_->granularPressureCoeffPrime()
      	 );
      */
      nut_ = viscosityModel_->calcNu();
      //sigNu_ = viscosityModel_->calcNu();
      //
      //sigNu_.internalField() = nut_.internalField();
      //sigNu_.correctBoundaryConditions();
      //nu(alpha, Theta_, gs0_, rho, da, e_);

        volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
	        lambda_ = scalar(0.0)*(4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;
	//lambda_ =(4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;
        // Frictional pressure
        volScalarField pf
        (
            frictionalStressModel_->frictionalPressure
            (
                phase_,
                alphaMinFriction_,
                alphaMax_
            )
        );

        nuFric_ = frictionalStressModel_->nu
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            pf/rho,
            D
        );

        // Limit viscosity and add frictional viscosity
        //nut_.min(maxNut_);
        //nuFric_ = min(nuFric_, maxNut_ - nut_);
        //nut_ += nuFric_;
    }

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //
