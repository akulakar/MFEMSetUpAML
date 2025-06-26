//                                  Isotropic Elasticity

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
   // 1. Parse command-line options.
   const char *mesh_file = "../mesh/CylindricalRod-hex.mesh";  
   int order = 1;
   bool static_cond = false;
   bool visualization = 0;
   int ref_levels = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&ref_levels, "-r", "--ref_levels",
                  "Number of uniform mesh refinements.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();


   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   //    Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.
   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
  
   fec = new H1_FECollection(order, dim);
   fespace = new FiniteElementSpace(mesh, fec, dim);

   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl ;

   Array<int> dofs;
   cout << "Element dofs: " << endl;
   fespace->GetElementVDofs(0, dofs);
   for (int i = 0; i < dofs.Size(); i++){
      cout << dofs[i] << " ";
   };
   cout << endl;

   // int numels = fespace->GetNE();
   // cout << "Mesh dimension: " << mesh->Dimension() << endl;

   // cout << " Number of elements: " << numels << endl;

   // fespace->GetElementVDofs(1,dofs);
   // cout << dofs[0] << ", " << dofs[1] << ", " << dofs[2] << ", " << dofs[3] << endl;
   // cout << dofs[4] << ", " << dofs[5] << ", " << dofs[6] << ", " << dofs[7] << endl;

   //    Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking only
   //    boundary attribute 1 from the mesh as essential and converting it to a
   //    list of true dofs.

   Array<int> ess_tdof_listx, ess_tdof_listy, ess_tdof_listz, 
   ess_bdr_x(mesh->bdr_attributes.Max()), ess_bdr_y(mesh->bdr_attributes.Max()), ess_bdr_z(mesh->bdr_attributes.Max());

   ess_bdr_x = 0;
   ess_bdr_x[0] = 1;
   // ess_bdr_x[2] = 1;

   ess_bdr_y = 0;  // 1 for Constrain the other surface too for the manufactured solution.
   ess_bdr_y[0] = 1;
   // ess_bdr_y[2] = 1;

   ess_bdr_z = 0; // 1 for Constrain the other surface too for the manufactured solution.
   ess_bdr_z[0] = 1;
   ess_bdr_z[1] = 1;

   fespace->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listx, 0);
   fespace->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listy, 1);
   fespace->GetEssentialTrueDofs(ess_bdr_z, ess_tdof_listz, 2);

   Array<int> ess_tdof_list;

   ess_tdof_list.Append(ess_tdof_listx); ess_tdof_list.Append(ess_tdof_listy); ess_tdof_list.Append(ess_tdof_listz);

   // Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
   // ess_bdr = 0;
   // ess_bdr[0] = 1;
   // ess_bdr[1] = 1;

   // fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);


   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system. In this case, b_i equals the boundary integral
   //    of f*phi_i where f represents a "pull down" force on the Neumann part
   //    of the boundary and phi_i are the basis functions in the finite element
   //    fespace. The force is defined by the VectorArrayCoefficient object f,
   //    which is a vector of Coefficient objects. The fact that f is non-zero
   //    on boundary attribute 2 is indicated by the use of piece-wise constants
   //    coefficient for its last component.

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

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.

   // GridFunction u_mfem(fespace);
   GridFunction u_mfemplus(fespace);
   // u_mfem = 0.0; 
   u_mfemplus = 0.0; 

   VectorArrayCoefficient dirichlet(dim);
   Vector zdisp (mesh->bdr_attributes.Size());
   zdisp(1) = -0.1;
   
   dirichlet.Set(2, new PWConstCoefficient(zdisp));
   dirichlet.Set(1, new ConstantCoefficient(0.0));
   dirichlet.Set(0, new ConstantCoefficient(0.0));

   // dirichlet.Set(0, new FunctionCoefficient ([](const Vector &xcoord) {
   //    return -1*xcoord(0)/30; 
   // }));
   // dirichlet.Set(1, 
   //    new FunctionCoefficient ([](const Vector &xcoord) {
   //    return (xcoord(0))/5;  })
   //    ); 
   //    dirichlet.Set(2, 
   //       new FunctionCoefficient ([](const Vector &xcoord) {
   //       return (xcoord(0))/5; })
   //    );
   

   // u_mfem.ProjectBdrCoefficient(dirichlet, ess_bdr_z);
   u_mfemplus.ProjectBdrCoefficient(dirichlet, ess_bdr_z);

      
   // u.ProjectBdrCoefficient(dirichlet, ess_bdr_x);
   // u.ProjectBdrCoefficient(dirichlet, ess_bdr_y);
   // u.ProjectBdrCoefficient(dirichlet, ess_bdr_z);


   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.

   // BilinearForm *a_mfem = new BilinearForm(fespace);

   // float lambda = 142.857;
   // float mu = 35.7143;
   // ConstantCoefficient lambda_coeff(lambda);
   // ConstantCoefficient mu_coeff(mu);

   // a_mfem->AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff,mu_coeff));

   BilinearForm *a_mfemplus = new BilinearForm(fespace);

   float young_mod = 100.0;
   float poisson_ratio = 0.4;
   ConstantCoefficient young_mod_func(young_mod);
   ConstantCoefficient poisson_ratio_func(poisson_ratio);

   a_mfemplus->AddDomainIntegrator(new mfemplus::ThreeDIsotropicElasticityIntegrator(young_mod_func,poisson_ratio_func));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.

   // a_mfem->Assemble(); 
   a_mfemplus->Assemble();

   cout << "Bilinear forms assembled" << endl;

   // SparseMatrix A_mfem;
   SparseMatrix A_mfemplus;
   Vector B_mfem, X_mfem;
   Vector B_mfemplus, X_mfemplus;

   // a_mfem->FormLinearSystem(ess_tdof_list, u_mfem, *b, A_mfem, X_mfem, B_mfem);
   a_mfemplus->FormLinearSystem(ess_tdof_list, u_mfemplus, *b, A_mfemplus, X_mfemplus, B_mfemplus);

   // GSSmoother M_mfem(A_mfem);
   GSSmoother M_mfemplus(A_mfemplus);
   // PCG(A_mfem, M_mfem, B_mfem, X_mfem, 0, 500, 1e-8, 0.0);
   PCG(A_mfemplus, M_mfemplus, B_mfemplus, X_mfemplus, 0, 500, 1e-8, 0.0);

   cout << "Linear system solved" << endl;

   // 12. Recover the solution as a finite element grid function.
   // a_mfem->RecoverFEMSolution(X_mfem, *b, u_mfem);
   a_mfemplus->RecoverFEMSolution(X_mfemplus,*b,u_mfemplus);

   
   
   // ofstream mfemoutput("output/mfemoutput.dat");
   // mfemoutput << u_mfemplus << endl;
   // mfemoutput.close();

   // // Now calculate strains and stress with the customfunctions.

   mfemplus::GlobalStressStrain StressStrain(mesh, fespace);
   GridFunction strain, stress;
   StressStrain.GlobalStrain(u_mfemplus, strain);
   StressStrain.GlobalStress(strain, young_mod_func, poisson_ratio_func, stress);

   cout << "Number of elements: " << fespace->GetNE() << endl;
   cout << "Size of strain GridFunction: " << strain.Size() << endl;
   cout << "Maximum strain: " << strain.Max() << endl;
   cout << "Size of stress GridFunction: " << stress.Size() << endl;
   cout << "Maximum stress: " << stress.Max() << endl;

   // // ofstream sig3output("output/sig33.dat");
   // // sig3output << stress << endl;
   // // sig3output.close();


   {
      double top_force = 0.0;
      double top_area = 0.0;
      array<double,6> top_strain, top_stress;
      int num_bdr_els = fespace->GetNBE();
      int num_els = fespace->GetNE();
      for (int i = 0; i < num_bdr_els; i++)
      {
         int bdr_attr = fespace->GetBdrAttribute(i);
         if (bdr_attr == 2)
         {
            FaceElementTransformations *ftr = fespace->GetMesh()->GetBdrFaceTransformations(i);
            int volume_elem_index = ftr->Elem1No;
            double element_sig = stress(volume_elem_index + (2 * num_els));
            double bdr_element_force = 0.0;
            real_t bdr_element_area = 0.0;
            const FiniteElement &surf_el = *(fespace->GetBE(i));
            ElementTransformation &el_trans = *(fespace->GetBdrElementTransformation(i));
            mfemplus::ElementStressStrain ElementArea(surf_el, el_trans);
            ElementArea.ComputeBoundaryElementArea(bdr_element_area);
            bdr_element_force = bdr_element_area * element_sig;
            top_force += bdr_element_force;
            top_area += bdr_element_area;

            for (int comp = 0; comp < 6; comp ++){
            top_strain[comp] = strain(volume_elem_index + (comp * num_els));
            top_stress[comp] = stress(volume_elem_index + (comp * num_els));
            }

         }
      }
      cout << "Force on top surface is: " << top_force << endl;
      cout << "The area of the top surface is: " << top_area << endl;
      cout << "Strain: " << top_strain[0] << ", " << top_strain[1] << ", " << top_strain[2] << ", " << top_strain[3] << ", " << top_strain[4] << ", " << top_strain[5] << ", "  << endl;
      cout << "Stress: " << top_stress[0] << ", " << top_stress[1] << ", " << top_stress[2] << ", " << top_stress[3] << ", " << top_stress[4] << ", " << top_stress[5] << ", "  << endl;
   }
      


   // // 15. Send the above data by socket to a GLVis server. Use the "n" and "b"
   // //     keys in GLVis to visualize the displacements.
   // if (visualization)
   // {
   //    char vishost[] = "localhost";
   //    int  visport   = 19916;
   //    socketstream sol_sock(vishost, visport);
   //    sol_sock.precision(8);
   //    sol_sock << "solution\n" << *mesh << u_mfemplus  << "window_title 'Displacement Field'" << flush;
   // }

   // //   Free the used memory.
   // delete a_mfem;
   delete a_mfemplus;
   delete b;
   delete fespace;
   delete fec;
   delete mesh;

   return 0;
}
