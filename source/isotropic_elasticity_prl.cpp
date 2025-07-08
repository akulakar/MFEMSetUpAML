//                                  Isotropic Elasticity Parallel

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "custombilininteg.hpp"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    //    Initialize MPI and HYPRE.
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

   //     Parse command-line options.
   const char *mesh_file = "../mesh/AnnularCylinder-tet.mesh";
   int order = 1;
   bool visualization = 0;
   int ref_levels = 0;
   int iterations = 1000;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&ref_levels, "-r", "--ref_levels",
                  "Number of uniform mesh refinements.");
   args.AddOption(&iterations, "-it", "--iterations",
                  "Number of solver iterations.");
   args.Parse();

   //    Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   //    Read the (serial) mesh from the given mesh file on all processors. 
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   for (int l = 0; l < ref_levels; l++)
   {
     mesh->UniformRefinement();
   }
   
   //    Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.

   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  
   //    Define a parallel and serial vector finite element space on the mesh.
   //    Serial fespace is for saving data at the end.
   FiniteElementCollection *fec;
   ParFiniteElementSpace *fespace;
   FiniteElementSpace *serfespace;
  
   fec = new H1_FECollection(order, dim);
   fespace = new ParFiniteElementSpace(pmesh, fec, dim);
   serfespace = new FiniteElementSpace(mesh, fec, dim);

   // Outputs from each processor that it is run on.
   // Finite element unknowns DOES NOT double count the degrees of freedom from parallel fespace.
   // Nodes and degrees of freedom DOES double count the degrees of freedom.

   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl ;
   cout << "Number of degrees of freedom: " << fespace->GetVSize() << endl ;
   cout << "Number of nodes: " << fespace->GetNDofs() << endl ;


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

   // serfespace->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listx_ser, 0);
   // serfespace->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listy_ser, 1);
   // serfespace->GetEssentialTrueDofs(ess_bdr_z, ess_tdof_listz_ser, 2);

   Array<int> ess_tdof_list, ess_tdof_list_ser;

   ess_tdof_list.Append(ess_tdof_listx); 
   ess_tdof_list.Append(ess_tdof_listy); 
   ess_tdof_list.Append(ess_tdof_listz);

   // ess_tdof_list_ser.Append(ess_tdof_listx_ser); 
   // ess_tdof_list_ser.Append(ess_tdof_listy_ser); 
   // ess_tdof_list_ser.Append(ess_tdof_listz_ser);

   cout << "Num of ess dofs: " << ess_tdof_list.Size() << endl;
   // cout << "Num of ess dofs (serial): " << ess_tdof_list_ser.Size() << endl;

   // Array<int> ess_tdof_list, ess_bdr(pmesh->bdr_attributes.Max());
   // ess_bdr = 0;
   // ess_bdr[0] = 1;
   // ess_bdr[1] = 1;
   // fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   //   Set up the linear form b(.) which corresponds to the right-hand side of
   //   the FEM linear system.

   // // VectorArrayCoefficient f(dim);
   // // for (int i = 0; i < dim; i++)
   // // {
   // //    f.Set(i, new ConstantCoefficient(0.0));
   // // }
   
   // // Vector pull_force(mesh->bdr_attributes.Size());
   // // pull_force = 0.0;
   // // pull_force(1) = 250000.0;

   // // f.Set(2, new PWConstCoefficient(pull_force));

   ParLinearForm* b = new ParLinearForm(fespace);
   // // b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));

   // // VectorArrayCoefficient bodyforce(dim);
   // // bodyforce.Set(0, new ConstantCoefficient(0.0));
   // // bodyforce.Set(1, new ConstantCoefficient(0.0));
   // // bodyforce.Set(2, new ConstantCoefficient(-10.0));

   // // b->AddDomainIntegrator(new VectorDomainLFIntegrator(bodyforce));

   b->Assemble();
   // cout << "RHS assembled" << endl;

   //    Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.

   ParGridFunction u(fespace);
   u = 0.0;
   cout << "ParGridFunction size: " << u.Size() << endl;

   VectorArrayCoefficient dirichlet(dim);
   Vector dispbcx(pmesh->bdr_attributes.Size()), dispbcy(pmesh->bdr_attributes.Size()), dispbcz(pmesh->bdr_attributes.Size());
   dispbcx = dispbcy = dispbcz = 0.0;
   dispbcz(1) = 1.0;

   dirichlet.Set(0, new PWConstCoefficient(dispbcx));
   dirichlet.Set(1, new PWConstCoefficient(dispbcy));
   dirichlet.Set(2, new PWConstCoefficient(dispbcz));

   u.ProjectBdrCoefficient(dirichlet, ess_bdr_x);
   u.ProjectBdrCoefficient(dirichlet, ess_bdr_y);
   u.ProjectBdrCoefficient(dirichlet, ess_bdr_z);

   //    Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.

   ParBilinearForm* a = new ParBilinearForm(fespace);

   float E = 100.0;
   float NU = 0.4;
   ConstantCoefficient E_func(E);
   ConstantCoefficient NU_func(NU);

   // float lambda = 142.857;
   // float mu = 35.7143;
   // ConstantCoefficient lambda_coeff(lambda);
   // ConstantCoefficient mu_coeff(mu);


   a->AddDomainIntegrator(new mfemplus::ThreeDIsotropicElasticityIntegrator(E_func, NU_func));
   // // a->AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff,mu_coeff));

   //  Assemble the bilinear form and the corresponding linear system.

   a->Assemble();  
   // cout << "Bilinear form assembled" << endl;

   HypreParMatrix A;
   Vector B, X;

   a->FormLinearSystem(ess_tdof_list, u, *b, A, X, B);

   HypreBoomerAMG *amg = new HypreBoomerAMG(A);
   amg->SetSystemsOptions(dim);
   amg->SetPrintLevel(0);

   HyprePCG *pcg = new HyprePCG(A);
   pcg->SetTol(1e-8);
   pcg->SetMaxIter(500);
   pcg->SetPrintLevel(0);
   pcg->SetPreconditioner(*amg);
   pcg->Mult(B, X);

   a->RecoverFEMSolution(X, *b, u);

   cout << "Linear systems solved" << endl;

   //  Saving data. Create serial gridfunction.
   ofstream dispdata("../results/IsotropicDisp_prl.dat");
   dispdata.precision(8);
   u.Save(dispdata);

   // // GridFunction *u_ser = pmesh->GetNodes();
   // // *u_ser += u_ani;
   // // // u_ani.ParallelProject(u_ser);
   // // ofstream outputfile("output/IsotropicDisp.dat");
   // // outputfile << u_ser << endl;


   // ParaViewDataCollection paraview_elasticity("isotropic_output", mesh);
   // paraview_elasticity.SetPrefixPath("results");
   // paraview_elasticity.SetLevelsOfDetail(order);
   // paraview_elasticity.SetDataFormat(VTKFormat::ASCII);
   // paraview_elasticity.SetHighOrderOutput(true);
   // paraview_elasticity.RegisterField("displacement",&u);
   // paraview_elasticity.Save();
  
   //   Free the used memory.
   delete a;
   delete b;
   delete amg;
   delete pcg;

   delete fespace;
   delete serfespace;
   delete fec;
   delete pmesh;

   return 0;
}

