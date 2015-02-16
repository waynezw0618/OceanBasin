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

Description
    numerical damping in the internalfield for selected boundaries

Class
    Foam::OceanBasin

Author
    Wei Zhang, XiaoDong Bai, 
    National University of Singapore.
    All rights reserved.

SourceFiles
    OceanBasin.C

\*---------------------------------------------------------------------------*/

#ifndef OceanBasin_H
#define OceanBasin_H

#include "volFields.H"
#include "IOdictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class numericalBeach Declaration
\*---------------------------------------------------------------------------*/

class OceanBasin
:
    public IOdictionary
{

    // Private data

        //- Velocity field to be damped at boundaries
        volVectorField& U_;
    
        //- wave parameters
        dimensionedScalar waveH_;
        dimensionedScalar waveT_;
        dimensionedScalar waveL_;
        dimensionedScalar waveD_; //water depth, h: in the papers;
    
        dimensionedScalar waveAmp_;
        dimensionedScalar waveOmega_;
        dimensionedScalar waveK_;
        dimensionedScalar waveC_;
    
        //wave maker paramters from O.S. Madesen, On the Generation of Long waves, December 20, 1971 ?? no need for private ddata?

        scalar waveTsoft_;
    
        Switch massSource_;
        Switch momentumSource_;
        Switch secondOrder_;
    
        // Read geometry parameters
        vector X_;
        vector Z_;
        tensor Tdamp_;
    
        scalar h1_; // h1=0
        scalar h2_; // h2> D+H/2 to include all cells in the vertical direction
        scalar deltaL_;
        dimensionedScalar Lx_;
        dimensionedScalar zoneCenterXcoor_;
        dimensionedScalar zoneCenterZcoor_;
    
        // for Beach current
        dimensionedScalar currentV_;
        scalar curTsoft_;
    
        volScalarField sourceZone_;
        dimensionedScalar OrgVolume_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OceanBasin(const OceanBasin&);

        //- Disallow default bitwise assignment
        void operator=(const OceanBasin&);
    
        void SelectZone();
    
        dimensionedScalar sourceZoneVolume() const;
    
        //- phase shift for second order paddle mass source
        void CalculatePhaseShift() const;
    
        scalar curTsoft() const;
    
        scalar waveTsoft() const;
    
        dimensionedScalar scaleFactor() const;


public:

    // Constructors

        //- Construct from components
        OceanBasin(volVectorField& U);


    // Destructor

        ~OceanBasin();


    // Member Functions

        // Access

        //- Return the source term for the momentum equation for beach
        tmp<volScalarField> BeachDamping() const;
    
        //- Return the source term for the momentum equation for beach
        tmp<volScalarField> currentSource() const;
    
        //- Return the source term for the momentum equation for beach
        tmp<volScalarField> massPaddle() const;
    
        //- Return the source term for the momentum equation for beach
        tmp<volVectorField> forcePaddle() const;
    
        Switch massSource() const
        {
            return massSource_;
        }
    
        Switch momentumSource() const
        {
            return momentumSource_;
        }
    
        tensor Tdamp() const
        {
            return Tdamp_;
        }



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
