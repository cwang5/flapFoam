/****************************************************************************\ 
 * Rotating a block by -10 degrees to 10 degrees in oscillating motion
 * This code segment is used to test the mesh motion files
\****************************************************************************/
// To Do List:
// Oscillate the block based on the value of runtime.timeOutputValue()
floatScalar realDeltaT = runTime.deltaT().value();
floatScalar currentTime = runTime.timeOutputValue();
//cout << "the current run time is " << currentTime << '\n';

// cellposition for a patch: mesh.C().boundaryField()[patchID]
// pointposition ofr a patch: mesh.boundaryMesh()[patchID].localPoints()

// forAllConstIter (labelHashSet, patchList_, iter)
// {
// 	label patchi = iter.key();
// 	vectorField oldCell = mesh.C().boundaryField()[patchi];
// 	pointField oldPoint = mesh.boundaryMesh()[patchi].localPoints();
// }

// oscillatory motion on the flap below, insert new line and comment out this one for forced motion

//  scalar theta = 45*(Foam::sin(currentTime*10)-Foam::sin((currentTime-realDeltaT)*10))*M_PI/180;
//  scalar theta = 45*(Foam::sin(currentTime*10))*M_PI/180;


// Insert forced motion solver below:


scalar Moment_Inertia = -(momentF + momentG)/flapInertia;
scalar beta_move = 0.5 * Moment_Inertia * realDeltaT * realDeltaT + beta_dot_0 * realDeltaT;
beta_0 += beta_move;
beta_dot_0 += Moment_Inertia * realDeltaT;

Info << '\n' << "The angle of rotation is " << beta_move << '\n';
Info << '\n' << "The new flap angle is " << beta_0 << '\n';


// Obtaining the rotational tensor:
// Calculating the RotTen
/*
scalar w = Foam::cos (0.5 * beta_move);
vector v = (Foam::sin (0.5 * beta_move)) * AoR;

Info << '\n' << "-1 * AoR =  " << -1*AoR << '\n';
Info << '\n' << "Scalar w is " << w << " vector v is " << v << '\n';

tensor RotTen_1 
(	sqr(w)+sqr(v.x())-sqr(v.y())-sqr(v.z()),
	2*(v.x()*v.y()-w*v.z()),
	2*(v.x()*v.z()+w*v.y()),
	2*(v.x()*v.y()+w*v.z()),
	sqr(w)-sqr(v.x())+sqr(v.y())-sqr(v.z()),
	2*(v.y()*v.z()-w*v.x()),
	2*(v.x()*v.z()-w*v.y()),
	2*(v.y()*v.z()+w*v.x()),
	sqr(w)-sqr(v.x())-sqr(v.y())+sqr(v.z()));
	
Info << '\n' << v.x()*v.y() <<"  " << -w*v.z() << "  " << w*v.z() << '\n';
*/

//  quaternion RotationMatrix(AoR,theta); // For oscillation
quaternion RotationMatrix(AoR,beta_0); // For fluid-induced motion
tensor RotTen = RotationMatrix.R();
Info << '\n' << "The Z component of the Axis of Rotation is " << AoR.z() << '\n' << endl;
Info << '\n' << "The rotational Tensor is " << RotTen << '\n' << endl;


forAllConstIter(labelHashSet, patchList_, iter)
{
	label patchi = iter.key();
	vectorField oldCell = mesh.C().boundaryField()[patchi];
	pointField oldPoint = mesh.boundaryMesh()[patchi].localPoints();
// Rotate the patch by RotTem
	vectorField newCell = ((oldCell-CoR) & RotTen) + CoR;
	pointField newPoint = ((oldPoint-CoR) & RotTen) + CoR;
// Cummulative point displacement
	pointDisplacement.boundaryField()[patchi] == (newPoint - oldPoint);
	cellDisplacement.boundaryField()[patchi] == (newCell - oldCell);
}

