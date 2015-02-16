/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "OceanBasin.H"
#include "numericalBeachFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void OceanBasin::SelectZone()
{
    volScalarField CoodX(U_.mesh().C()&X_);
    volScalarField CoodZ(U_.mesh().C()&Z_);
    
    sourceZone_ = pos(
                            (deltaL_*waveL_/2)-mag(CoodX-zoneCenterXcoor_)
                      )
                  *pos(
                            ((h1_-h2_)*waveD_/2)-mag(CoodZ-(waveD_-zoneCenterZcoor_))
                       );
         
}
    
dimensionedScalar OceanBasin::sourceZoneVolume() const
{
    return sum(sourceZone_*U_.mesh().V()); //NOT parallel
}

void OceanBasin::CalculatePhaseShift() const
{
     // do nothing?
}

scalar OceanBasin::curTsoft() const
{
    return scalar(1.0);
}

scalar OceanBasin::waveTsoft() const
{
    return scalar(1.0);
}
    
dimensionedScalar OceanBasin::scaleFactor() const
{
    return OrgVolume_/sourceZoneVolume();
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

OceanBasin::OceanBasin(volVectorField& U)
:
    U_(U),
    IOdictionary
    (
        IOobject
        (
            "waveProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
         )
     ),
    waveH_(lookup("waveH")),
    waveT_(lookup("waveT")),
    waveL_(lookup("waveL")),
    waveD_(lookup("waveD")),
    waveAmp_(0.5*waveH_),
    waveOmega_(2*mathematicalConstant::pi/waveT_),
    waveK_(2*mathematicalConstant::pi/waveL_),
    waveC_(waveL_/waveT_),
    waveTsoft_(lookupOrDefault<scalar>("waveTsoft",scalar(0.0))),
    massSource_(lookupOrDefault<Switch>("massSource",true)), //false
    momentumSource_(lookupOrDefault<Switch>("momentumSource",false)), //false
    secondOrder_(lookupOrDefault<Switch>("SecondOrder",false)),
    X_(lookupOrDefault<vector>("streamDirection",vector(1, 0, 0))),
    Z_(lookupOrDefault<vector>("streamDirection",vector(0, 1, 0))),
    Tdamp_(lookupOrDefault<tensor>("beachDirection",tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)-X_*X_)),
    h1_(lookupOrDefault<scalar>("h1",scalar(2.0))), // h1=0
    h2_(lookupOrDefault<scalar>("h2",scalar(0.0))), // h2> D+H/2 to include all cells in the vertical direction
    deltaL_(lookupOrDefault<scalar>("detL",scalar(0.016666))), // default value based on 60 cells/wave length
    Lx_(deltaL_* waveL_), //according to Gauss Theory, Source =2U/Lx,
    zoneCenterXcoor_(lookupOrDefault<dimensionedScalar>("zoneCenterXcoor",dimensionedScalar("Xc",dimLength,0.0))),
    zoneCenterZcoor_((h1_+h2_)/2*waveD_),
    currentV_(lookupOrDefault<dimensionedScalar>("currentV",dimensionedScalar("cv",dimVelocity,0))),
    curTsoft_(lookupOrDefault<scalar>("curTsoft",0.0)),
    sourceZone_(
                 IOobject
                 (
                    "sourceZone",
                    U_.db().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                  ),
                    U_.mesh(),
                    dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0)
                ),
    OrgVolume_(Lx_*waveD_*Lx_)
{
    
    Info << "git test" <<endl;
    //- wave parameters
   /* waveH_ = waveProperties.lookup("waveH");
    waveT_ = waveProperties.lookup("waveT");
    waveL_ = waveProperties.lookup("waveL");
    waveD_ = waveProperties.lookup("waveD"); //water depth, h: in the papers;
 
    waveAmp_    = 0.5*waveH_;
    waveOmega_  = 2*mathematicalConstant::pi/waveT_;
    waveK_      = 2*mathematicalConstant::pi/waveL_;
    waveC_      = waveL_/waveT_;
 
    //wave maker paramters from O.S. Madesen, On the Generation of Long waves, December 20, 1971 ?? no need for private ddata?
    //dimensionedScalar n1_; //eqn (13) in Madesen' spaper
    //dimensionedScalar zeta0_;//eqn (13) in Madesen' spaper
    waveTsoft_ = waveProperties.lookupOrDefault<scalar>("waveTsoft",scalar(0.0));
    
    massSource_ = waveProperties.lookupOrDefault<Switch>("massSource",true)); //false
    
    if (massSource_){
        momentumSource_ = false;
    }
    else {
        momentumSource_ = true;
    }
    
    secondOrder = waveProperties.lookupOrDefault<Switch>("SecondOrder",false);
   
    //- Geometry parameters
    X_     = waveProperties.lookupOrDefault<vector>("streamDirection",vector(1, 0, 0));
    Z_     = waveProperties.lookupOrDefault<vector>("streamDirection",vector(0, 1, 0));
    Tdamp_ = waveProperties.lookupOrDefault<tensor>("beachDirection",tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)-X*X);
    
    h1_ ＝ waveProperties.lookupOrDefault<scalar>("h1",scalar(2.0)); // h1=0
    h2_ ＝ waveProperties.lookupOrDefault<scalar>("h2",scalar(0.0)); // h2> D+H/2 to include all cells in the vertical direction
    deltaL_ ＝ waveProperties.lookupOrDefault<scalar>("detL",scalar(0.016666)); // default value based on 60 cells/wave length
    Lx_ = deltaL_* waveL_; //according to Gauss Theory, Source =2U/Lx,
    //Ly_ = waveProperties.lookupOrDefault<dimensionedScalar>("DomainThickness",dimensionedScalar("Ly",dimLength,0.1));  // spanwise direction size of domain ( thickness of domainsize)
    zoneCenterXcoor_ ＝ waveProperties.lookupOrDefault<dimensionedScalar>("zoneCenterXcoor",dimensionedScalar("Xc",dimLength,0.0));
    zoneCenterZcoor_ = (h1_+h2_)/2*waveD_;
   
    
    // for Beach direction.
    currentV_ ＝ waveProperties.lookupOrDefault<dimensionedScalar>("currentV",dimensionedScalar("cv",dimVelocity,0));
    curTsoft_ = waveProperties.lookupOrDefault<scalar>("currentTsoft",scalar(0.0));
  */
    Info << "Initialise OceanBasin" << endl;
    this->SelectZone();
    OrgVolume_ = this->sourceZoneVolume();
    Info<<"OceanBasin is prepared ... \n"
        <<"wave amptitude is   :" << waveAmp_ << "\n"
        <<"wave number is      :" << waveK_ << "\n"
        <<"wave frequecy is    :" << waveOmega_ <<"\n"
        <<"water depth is      :" << waveD_ << "\n"
        <<"sourceZoneVolume is :" << OrgVolume_<< endl;

}


OceanBasin::~OceanBasin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> OceanBasin::BeachDamping() const
{
    tmp<volScalarField> tDamping
    (
        new volScalarField
        (
            IOobject
            (
                "damping",
                U_.db().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0)
        )
    );

    volScalarField& damping = tDamping();

    const volVectorField::GeometricBoundaryField& bvf = U_.boundaryField();
    forAll(bvf, patchi)
    {
        if (isType<numericalBeachFvPatchField>(bvf[patchi]))
        {
            const numericalBeachFvPatchField& beach =
                dynamic_cast<const numericalBeachFvPatchField&> (bvf[patchi]);

            damping = max(damping, beach.internalDamping());
        }
    }
    
    return tDamping;
}
    
tmp<volScalarField> OceanBasin::currentSource() const
{
    tmp<volScalarField> tCurrentSource
    (
        new volScalarField
        (
            IOobject
            (
                "CurrrentSource",
                U_.db().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
            U_.mesh(),
            2*currentV_/Lx_
         )
     );
    
    volScalarField& CurSource = tCurrentSource();

    return this->scaleFactor()*CurSource*sourceZone_;
    
}
    
tmp<volScalarField> OceanBasin::massPaddle() const
{
    
    //wave maker paramters from O.S. Madesen, On the Generation of Long waves, December 20, 1971
    dimensionedScalar n1(0.5*(1+2*waveK_*waveD_/Foam::sinh(2*waveK_*waveD_))); //eqn (13) in Madesen' spaper
    dimensionedScalar zeta0(waveAmp_*n1/Foam::tanh(waveK_*waveD_));//eqn (13) in Madesen' spaper
    tmp<volScalarField> tmassPaddle
    (
        new volScalarField
        (
            IOobject
            (
                "massPaddle",
                U_.db().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
             U_.mesh(),
            2*zeta0*waveOmega_/Lx_
        )
    );
        
        volScalarField& massPad   = tmassPaddle();

    
    if(!secondOrder_){
        return scaleFactor()*massPad*sourceZone_;
    }
    else{
        
        dimensionedScalar SecOrderTerm = 2*zeta0*waveOmega_*(waveAmp_/(waveD_*n1))
                                        *(0.75/Foam::sqr(Foam::sinh(waveK_*waveD_))-n1/2)/Lx_;
        return scaleFactor()*(massPad*sourceZone_+ SecOrderTerm*sourceZone_);
    }
    
}

tmp<volVectorField> OceanBasin::forcePaddle() const
{
        // haven't implemented yet
    tmp<volVectorField> tforcePaddle
    (
        new volVectorField
        (
            IOobject
            (
                "forcePaddle",
                U_.db().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
            U_.mesh(),
            dimensionedVector("test",dimless,vector(1.0,0,0))
      )
     );
    
    volVectorField& m_forcePad   = tforcePaddle();
    return m_forcePad;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
