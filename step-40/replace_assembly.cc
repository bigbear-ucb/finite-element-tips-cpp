template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  TimerOutput::Scope t(computing_timer, "assembly");
  
  const QGauss<dim>  quadrature_formula(3);
  
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points |
                           update_JxW_values);
  
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  
  std::vector<FullMatrix<double> >
  cell_mat_at_qp (n_q_points, FullMatrix<double> (dofs_per_cell, dofs_per_cell));
  
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  bool is_preassembly_finished = false;
  
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
    {
      cell_matrix = 0;
      cell_rhs = 0;
      
      fe_values.reinit (cell);
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      // This is what I have added
      if (!is_preassembly_finished)
      {
        for (unsigned int qi=0; qi<n_q_points; ++qi)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_mat_at_qp[qi](i,j) += (fe_values.shape_grad(i,qi) *
                                          fe_values.shape_grad(j,qi));
        is_preassembly_finished = true;
      }
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        const double
        rhs_value
        = (fe_values.quadrature_point(q_point)[1]
           >
           0.5+0.25*std::sin(4.0 * numbers::PI *
                             fe_values.quadrature_point(q_point)[0])
           ? 1 : -1);
        
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            cell_matrix(i,j) += (cell_mat_at_qp[q_point](i,j) *
                                 fe_values.JxW(q_point));
          /*
           This is the original assembly
             cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                 fe_values.shape_grad(j,q_point) *
                                 fe_values.JxW(q_point));
          */
          
          cell_rhs(i) += (rhs_value *
                          fe_values.shape_value(i,q_point) *
                          fe_values.JxW(q_point));
        }
      }
      
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix,
                                              cell_rhs,
                                              local_dof_indices,
                                              system_matrix,
                                              system_rhs);
    }
  
  // Notice that the assembling above is just a local operation. So, to
  // form the "global" linear system, a synchronization between all
  // processors is needed. This could be done by invoking the function
  // compress(). See @ref GlossCompress  "Compressing distributed objects"
  // for more information on what is compress() designed to do.
  system_matrix.compress (VectorOperation::add);
  system_rhs.compress (VectorOperation::add);
}
