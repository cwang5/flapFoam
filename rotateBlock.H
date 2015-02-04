/****************************************************************************\ 
 * Rotating a block by -10 degrees to 10 degrees in oscillating motion
 * This code segment is used to test the mesh motion files
\****************************************************************************/
// To Do List:
// Oscillate the block based on the value of runtime.timeOutputValue()
floatScalar realDeltaT = runTime.deltaT().value();


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
// beta_dot_0 += Moment_Inertia * realDeltaT;
beta_dot_0 = beta_move/realDeltaT;
if (beta_dot_0 > 8e-4)
{
    beta_dot_0 = 8e-4;
    Info << '\n' << "Angular velocity too high, capped at 1e-3 deg/s." << '\n';
}

//Info << '\n' << "The angle of rotation is " << beta_move << "with rotation speed" << beta_dot_0 << "degree per seconds." << '\n';
Info << '\n' << "The new flap angle is " << beta_0 << '\n';


// Obtaining the rotational tensor:
// Calculating the RotTen

//  quaternion RotationMatrix(AoR,theta); // For oscillation
quaternion RotationMatrix(AoR,beta_0); // For fluid-induced motion
tensor RotTen = RotationMatrix.R();
// Info << '\n' << "The Z component of the Axis of Rotation is " << AoR.z() << '\n' << endl;
// Info << '\n' << "The rotational Tensor is " << RotTen << '\n' << endl;


forAllConstIter(labelHashSet, patchList_, iter)
{
	label patchi = iter.key();
        vectorField oldCell = mesh.C().boundaryField()[patchi];
	pointField oldPoint = mesh.boundaryMesh()[patchi].localPoints();
// Rotate the patch by RotTem
        vectorField newCell = ((oldCell-CoR) & RotTen) + CoR;
	pointField newPoint = ((oldPoint-CoR) & RotTen) + CoR;
// Cummulative point displacement
	pointDisplacement.boundaryField()[patchi] = (newPoint - oldPoint);
        cellDisplacement.boundaryField()[patchi] == (newCell - oldCell);
}

