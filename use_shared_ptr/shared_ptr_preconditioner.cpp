// This is a demonstration based upon my own implementation

// **********************************************************
// Assume Petsc AMG preconditioner is used;
// Denote the shared_ptr by precond_ptr
// Denote the global matrix by system_matrix
// **********************************************************

using namespace dealii;

template <int dim>
Problem<dim>::preprocess_preconditioner ()
{
  // reset the pointer
  precond_ptr.reset ();
  
  // assign the pointer
  precond_ptr = std_cxx11::shared_ptr<PETScWrappers::MPI::PreconditionAMG> (new PETScWrappers::MPI::PreconditionAMG);
  
  // The following are how we construct a preconditioner in a usual way
  PETScWrappers::MPI::PreconditionAMG::AdditionalData data;
  
  // for symmetric system:
  data.symmetric_operator = true;
  
  // initialize the preconditioner:
  precond_ptr->initialize (system_matrix, data);
}
