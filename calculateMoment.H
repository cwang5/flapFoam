// This section of the code calculates the fluid-induced moment acting on the
// flap as well as (if applicable) the gravitation induced moment.
// The moment will be calculated under similar formulation as the code used
// for the openFoam force.C function.
Info << '\n' << "Initializing Moment Calculation..." << '\n';


vector sumMoment = vector(0,0,0);
/*if (obr_.foundObject<incompressible::RASModel>("RASProperties"))
{
	const incompressible::RASModel& ras = obr_.lookupObject<incompressible::RASModel>("RASProperties");
} */
//tmp<volSymmTensorField> tdevRhoReff = Foam::forces::devRhoReff();
//

tmp<volSymmTensorField> tdevRhoReff = rhoFluid* turbulence->devReff();
const volSymmTensorField::GeometricBoundaryField& devRhoReffb = tdevRhoReff().boundaryField();


Info << '\n' << "Calculating Moment..." << '\n';

forAllConstIter(labelHashSet, patchList_, iter)
{
	// replace rho(p) and devRhoReffb[]
	label patchi = iter.key();
	vectorField Md = mesh.C().boundaryField()[patchi] - CoR;
	vectorField pf = mesh.Sf().boundaryField()[patchi]*p.boundaryField()[patchi];
	sumMoment += gSum(Md ^ pf);
//	vector vm = sum (Md ^ (-U.boundaryField()[patchi].snGrad()*Sfb[patchi]*1.55e-5));
//	Tensor<double> vf = (-U.boundaryField()[patchi].snGrad()*Sfb[patchi]* 1.55e-5);
	vectorField vf = mesh.Sf().boundaryField()[patchi] & devRhoReffb[patchi];
	sumMoment += gSum(Md ^ vf);
//	sumMoment += vm;
}

momentF = sumMoment & AoR;

Info << '\n' << "The force moment acting on the flap is " << momentF << '\n';
