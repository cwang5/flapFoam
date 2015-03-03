Info << "Reading the Extra Dict File" << endl;
IOdictionary newDict
(
 IOobject
 (
  "flapDict",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
 )
);
//
// New addition to output flap angle
OFstream str(runTime.path()/"angleTime.txt");

// for information on the use of "wordList" to read in patch names,
// please refer to openFoam function forces::read()

labelHashSet patchList_=mesh.boundaryMesh().patchSet(wordReList(newDict.lookup("patches")));

const wordList motionPatches(newDict.lookup("patches"));

//const surfaceVectorField::GeometricBoundaryField& Sfb = mesh.Sf().boundaryField();

scalar pArea = 0;

forAllConstIter(labelHashSet, patchList_, iter)
{
	label patchi = iter.key();
	pArea += gSum(mag(mesh.Sf().boundaryField()[patchi]));
//	pArea += patchArea.value();
}
Info << '\n' << "The area of the moving patch is " << pArea/2 << '\n';
// for information on the use of "readScalar" please consult openFoam function "readTimeControls.H"
vector AoR = newDict.lookup("AxisOfRotation");
AoR = AoR/mag(AoR);

vector CoR = newDict.lookup("CenterOfRotation");

const scalar flapRho
(
 readScalar(newDict.lookup("flapDensity"))
);

const scalar flapL
(
 readScalar(newDict.lookup("flapLength"))
);

Info << '\n' << "The flap length is " << flapL << '\n';

const scalar rhoFluid
(
 readScalar(newDict.lookup("FluidDensity"))
);

// thickness of the flap and density of the flap are used to calculate the moment of inertia
// All measurement and values are scaled up to match the chord length of 1 meter, the flap
// material is assumed to be a 0.5cm thick carbon-fiber composite plate
const scalar flapThickness = 0.005;

// For a rectangular board the moment of inertia along an edge is I = (1/3)Ma^2, where a is the
// length and M is the mass

const scalar flapMass = pArea*flapRho*flapThickness;

Info << '\n'<< "The flap mass is " << flapMass << '\n';

const scalar flapInertia = 0.33333  * flapMass * flapL * flapL;

Info << '\n' << "The flap inertia is " << flapInertia << '\n';

// Initialize the flap angle, the initial angular velocity, the initial gravitational moment,
// and the initial fluid induced moment

scalar beta_0 
(
 readScalar(newDict.lookup("initialBeta"))
);
scalar beta_dot_0 = 0;
scalar momentG = 0;
scalar momentF = 0;

// Set angle limit based on Dict input

scalar beta_limit
(
 readScalar(newDict.lookup("betaLimit"))
)
;



// Initialization for the cellDisplacement and pointDisplacement portion
volVectorField& cellDisplacement = const_cast<volVectorField&>
(
 mesh.objectRegistry::lookupObject<volVectorField>
 (
  "cellDisplacement"
 )
);
pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
 mesh.objectRegistry::lookupObject<pointVectorField>
 (
  "pointDisplacement"
 )
);

//scalar counter = 0;
