//                                  Isotropic Elasticity Parallel

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "../../mfemplus/custombilininteg.hpp"
#include "../../mfemplus/customfunctions.hpp"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{

  // Initialize MPI and HYPRE.
  Mpi::Init(argc, argv);
  int num_procs = Mpi::WorldSize();
  int myid = Mpi::WorldRank();
  Hypre::Init();

  // Parse command-line options.
  const char *mesh_file = "../mesh/AnnularCylinder-tet.mesh";
  int order = 1;
  bool visualization = 0;
  int ref_levels = 0;
  int iterations = 1000;
  const char *device_config = "cpu";

  OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
                "Finite element order (polynomial degree).");
  args.AddOption(&ref_levels, "-r", "--ref_levels",
                "Number of uniform mesh refinements.");
  args.AddOption(&iterations, "-it", "--iterations",
                "Number of solver iterations.");
  args.Parse();

  // Enable hardware devices such as GPUs, and programming models such as
  // CUDA, OCCA, RAJA and OpenMP based on command line options.
  Device device(device_config);
  if (myid == 0) { device.Print(); }

  // Read the (serial) mesh from the given mesh file on all processors. 
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  for (int l = 0; l < ref_levels; l++)
  {
    mesh->UniformRefinement();
  }
  
  // Define a parallel mesh by a partitioning of the serial mesh. Refine this mesh further in 
  // parallel to increase the resolution. Once the parallel mesh is defined, the serial mesh can be deleted.

  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;

  // Define a parallel and serial vector finite element space on the mesh.
  // Serial fespace is for saving data at the end.
  FiniteElementCollection *fec;
  ParFiniteElementSpace *fespace;

  fec = new H1_FECollection(order, dim);
  fespace = new ParFiniteElementSpace(pmesh, fec, dim);

  // Determine the essential (Dirichlet) boundary dofs.

  Array<int> ess_tdof_listx, ess_tdof_listy, ess_tdof_listz, 
  ess_tdof_listx_ser, ess_tdof_listy_ser, ess_tdof_listz_ser,
  ess_bdr_x(pmesh->bdr_attributes.Max()), ess_bdr_y(pmesh->bdr_attributes.Max()), ess_bdr_z(pmesh->bdr_attributes.Max());

  ess_bdr_x = 0;
  ess_bdr_x[0] = 1;

  ess_bdr_y = 0; 
  ess_bdr_y[0] = 1;

  ess_bdr_z = 0; 
  ess_bdr_z[0] = 1;
  ess_bdr_z[1] = 1;

  fespace->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listx, 0);
  fespace->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listy, 1);
  fespace->GetEssentialTrueDofs(ess_bdr_z, ess_tdof_listz, 2);

  Array<int> ess_tdof_list, ess_tdof_list_ser;

  ess_tdof_list.Append(ess_tdof_listx); 
  ess_tdof_list.Append(ess_tdof_listy); 
  ess_tdof_list.Append(ess_tdof_listz);


  //   Set up the linear form b(.) which corresponds to the right-hand side of the FEM linear system.

  // Surface traction
  // VectorArrayCoefficient f(dim);
  // for (int i = 0; i < dim; i++)
  // {
  //    f.Set(i, new ConstantCoefficient(0.0));
  // }
  
  // Vector pull_force(mesh->bdr_attributes.Size());
  // pull_force = 0.0;
  // pull_force(1) = 250000.0;

  // f.Set(2, new PWConstCoefficient(pull_force));

  ParLinearForm* b = new ParLinearForm(fespace);
  // b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));

  // VectorArrayCoefficient bodyforce(dim);
  // bodyforce.Set(0, new ConstantCoefficient(0.0));
  // bodyforce.Set(1, new ConstantCoefficient(0.0));
  // bodyforce.Set(2, new ConstantCoefficient(-10.0));

  // b->AddDomainIntegrator(new VectorDomainLFIntegrator(bodyforce));

  b->UseFastAssembly(true);
  b->Assemble();

  cout << "RHS assembled" << endl;

  // Define the solution vector x as a finite element grid function corresponding to fespace. 
  // Initialize x with an initial guess of zero, which satisfies the boundary conditions.

  ParGridFunction u(fespace);
  u = 0.0;

  VectorArrayCoefficient dirichlet(dim);
  Vector dispbcz(pmesh->bdr_attributes.Size());
  dispbcz = 0.0;
  dispbcz(1) = -0.1;

  dirichlet.Set(0, new ConstantCoefficient(0.0));
  dirichlet.Set(1, new ConstantCoefficient(0.0));
  dirichlet.Set(2, new PWConstCoefficient(dispbcz));

  u.ProjectBdrCoefficient(dirichlet, ess_bdr_z);

  // Set up the bilinear form a(.,.) on the finite element space. 
  // Use the ThreeDIsotropicElasticityIntegrator from mfemplus.

  ParBilinearForm* a = new ParBilinearForm(fespace);

  float E = 100.0;
  float NU = 0.4;
  ConstantCoefficient E_func(E);
  ConstantCoefficient NU_func(NU);

  a->AddDomainIntegrator(new mfemplus::ThreeDIsotropicElasticityIntegrator(E_func, NU_func));
  
  //  Assemble the bilinear form and the corresponding linear system.

  a->Assemble();  
  
  // cout << "Bilinear form assembled" << endl;

  HypreParMatrix A;
  Vector B, X;

  a->FormLinearSystem(ess_tdof_list, u, *b, A, X, B);

  // Use a Hypre's Conjugate Gradient solver to solve the linear system with Hypre's Algebraic Multigrid preconditioner.

  HypreBoomerAMG *amg = new HypreBoomerAMG(A);
  amg->SetSystemsOptions(dim);
  amg->SetPrintLevel(0);

  HyprePCG *pcg = new HyprePCG(MPI_COMM_WORLD);
  pcg->SetOperator(A);
  pcg->SetTol(1e-8);
  pcg->SetMaxIter(500);
  pcg->SetPrintLevel(1);
  pcg->SetPreconditioner(*amg);
  pcg->Mult(B, X);

  a->RecoverFEMSolution(X, *b, u);


//  Optional. Save displacement data.
//  ofstream dispdata("../results/IsotropicDisp_prl.dat");
//  dispdata.precision(8);
//  u.Save(dispdata);

// Create ParGridFunctions for stress and strain. First create fespaces.

// Set up L2 FESpace to obtain one dof per element for each stress and strain component.
int numels = fespace->GetNE();
int str_comp = dim * 2;

mfem::L2_FECollection* L2fec;
mfem::ParFiniteElementSpace* L2fespace;

L2fec = new mfem::L2_FECollection(0, dim);
L2fespace = new mfem::ParFiniteElementSpace(pmesh, L2fec, str_comp);

ParGridFunction strain(L2fespace), stress(L2fespace);
strain = stress = 0.0;

// Create the GlobalStressStrain object and compute strains and stresses.
mfemplus::GlobalStressStrain StressStrain(pmesh, fespace);
StressStrain.GlobalStrain(u, strain);
StressStrain.GlobalStress(strain, E_func, NU_func, stress);

// Output max strain and stress from each parallel process.
cout << "Maximum strain: " << strain.Max() << endl;
cout << "Maximum stress: " << stress.Max() << endl;


  //   Free the used memory.
  delete a;
  delete b;
  delete amg;
  delete pcg;
  delete fespace;
  delete fec;
  delete pmesh;

  return 0;
}

