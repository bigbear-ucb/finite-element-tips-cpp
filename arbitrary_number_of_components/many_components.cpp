// ***************************************************
// Assume n_tot_comp is the total number of components,
// n_tot_comp > 100, and denote the finite element order
// by fe_order
// ***************************************************

using namespace dealii;

template <int dim>
class Problem ()
{
public:
  Problem ();
  ~Problem ();
  
  void run ();
private:
  void initialize_components ();
  void extract_to_local (FEValue<dim> &fv,
                         Vector<double> &global_solution,
                         std::vector<double> &local_solution,
                         unsigned int &desired_component);
  
  FESystem<dim> fe;
  
  // We need a container to contain all the component masks.
  // Each mask is a scalar.
  std::vector<FEValuesExtractors::Scalar> comp;
}

// We initialize the FESystem
template <int dim>
Problem<dim>::Problem ()
:
fe(FE_Q<dim>(fe_order), n_tot_comp),
dof_handler(triangulation)
{}

template <int dim>
void Problem<dim>::initialize_components ()
{
  // **************************************************
  // Here is how we set up FESystem and initialize the
  // component masks.
  for (unsigned int i=0; i<n_tot_comp; ++i)
  {
    // Assign a temporary mask to component with number i
    FEValuesExtractors::Scalar tmp(i);
    comp.push_back (tmp);
  }
  // **************************************************
}

// Here's how we use the component mask: we show how you
// extract solution to local cell for a specific component i
template <int dim>
void Problem<dim>::extract_to_local (FEValue<dim> &fv,
                                     Vector<double> &global_solution,
                                     std::vector<double> &local_solution,
                                     unsigned int &desired_component)
{
  fv[comp[desired_component]].get_function_values (global_solution,
                                                   local_solution);
}
