//                                  Isotropic Elasticity

#include "mfem.hpp"
#include <filesystem> // C++17
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip> // at the top of your file if not already included
#include <cmath>   // Include cmath for pow function
#include <map>
#include <initializer_list>
#include <stdexcept>
#include "../../mfemplus/mfemplus.hpp"

using namespace std;
using namespace mfem;
using json = nlohmann::json;
constexpr double π = M_PI;

int main(int argc, char *argv[])
{
   // Parse command-line options.
   const char *mesh_file = "../mesh/AnnularCylinder-tet.mesh";  
   int order = 1;
   bool visualization = 0;
   int iterations = 1000;
   int ref_levels = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&ref_levels, "-r", "--ref_levels",
                  "Number of uniform mesh refinements.");
   args.AddOption(&iterations, "-it", "--iterations",
                     "Number of solver iterations.");
   args.Parse();

   // Read the mesh from the given mesh file. 
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // Refine the mesh.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   //  Define a finite element space on the mesh. Here we use a vector H1 finite space. L2 space for stress and strain.
   FiniteElementCollection *fec, *l2fec;
   FiniteElementSpace *fespace, *l2fespace;
  
   fec = new H1_FECollection(order, dim);
   fespace = new FiniteElementSpace(mesh, fec, dim);

   int str_comp = (dim == 2) ? 3 : 6;

   l2fec = new L2_FECollection(0, dim);
   l2fespace = new FiniteElementSpace(mesh, l2fec, str_comp);

   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl ;

   // Determine the list of essential (Dirichlet) boundary dofs in each vector dimension.
   // Boundary attributes marked with 1 correspond to essential boundaries.
   // bdr_attribute 1 is bottom, 2 is top, and 3 are radial boundaries.

   Array<int> ess_tdof_listx, ess_tdof_listy, ess_tdof_listz, 
   ess_bdr_x(mesh->bdr_attributes.Max()), ess_bdr_y(mesh->bdr_attributes.Max()), ess_bdr_z(mesh->bdr_attributes.Max());

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

   Array<int> ess_tdof_list;

   ess_tdof_list.Append(ess_tdof_listx); ess_tdof_list.Append(ess_tdof_listy); ess_tdof_list.Append(ess_tdof_listz);

   // Set up the linear form b(.) which corresponds to the right-hand side of the FEM linear system. 

   // Surface traction   
   // VectorArrayCoefficient f(dim);
   // for (int i = 0; i < dim; i++)
   // {
   //    f.Set(i, new ConstantCoefficient(0.0));
   // }
   
   // Vector pull_force(mesh->bdr_attributes.Size());
   // pull_force = 0.0;
   // pull_force(1) = 25.0;

   // f.Set(0, new PWConstCoefficient(pull_force));

   LinearForm *b = new LinearForm(fespace);
   // b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));

   b->Assemble();
   cout << "RHS assembled" << endl;

   // Define the solution vector x as a finite element grid function corresponding to fespace. 
   // Initialize x with an initial guess of zero, which satisfies the boundary conditions.

   GridFunction u(fespace);
   u = 0.0; 

   // Displacement condition.
   VectorArrayCoefficient dirichlet(dim);
   Vector zdisp (mesh->bdr_attributes.Size());
   zdisp = 0.0;
   zdisp(1) = -0.1;
   
   dirichlet.Set(0, new ConstantCoefficient(0.0));
   dirichlet.Set(1, new ConstantCoefficient(0.0));
   dirichlet.Set(2, new PWConstCoefficient(zdisp)); // Non zero displacement in z direction on the top surface. 

   u.ProjectBdrCoefficient(dirichlet, ess_bdr_z);

   // Set up the bilinear form a(.,.) on the finite element space. 
   // Use the ThreeDIsotropicElasticityIntegrator from mfemplus.

   BilinearForm *a = new BilinearForm(fespace);

   float young_mod = 100.0;
   float poisson_ratio = 0.4;
   ConstantCoefficient young_mod_func(young_mod);
   ConstantCoefficient poisson_ratio_func(poisson_ratio);

   a->AddDomainIntegrator(new mfemplus::IsotropicElasticityIntegrator(young_mod_func,poisson_ratio_func));

   // Assemble the bilinear form and the corresponding linear system. Elimate boundary conditions.

   a->Assemble();

   cout << "Bilinear forms assembled" << endl;

   SparseMatrix A;
   Vector B, X;

   a->FormLinearSystem(ess_tdof_list, u, *b, A, X, B);

   // Use a Conjugate Gradient solver to solve the linear system with Gauss-Seidel preconditioner.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, iterations, 1e-8, 0.0);

   cout << "Linear system solved" << endl;

   // Recover the solution as a finite element grid function and save it.
  
   a->RecoverFEMSolution(X,*b,u);
   
   //  Optional. Save displacement data.
   //  ofstream dispdata("../results/IsotropicDisp.dat");
   //  dispdata.precision(8);
   //  u.Save(dispdata);

   // Create GlobalStressStrain object to calculate strains and stresses.

   mfemplus::GlobalStressStrain StressStrain(mesh, fespace);
   GridFunction strain(l2fespace), stress(l2fespace);
   strain = 0.0; stress = 0.0;
   StressStrain.GlobalStrain(u, strain);
   StressStrain.GlobalStress(strain, young_mod_func, poisson_ratio_func, stress);

   cout << "Number of elements: " << fespace->GetNE() << endl;
   cout << "Maximum strain: " << strain.Max() << endl;
   cout << "Maximum stress: " << stress.Max() << endl;

   // Free the used memory.
   delete a;
   delete b;
   delete fespace;
   delete fec;
   delete mesh;

   return 0;
}
