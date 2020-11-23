#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "boost/multi_array.hpp"

using std::cout;
using std::endl;

// Diagonalize real symmetric matrix, which is represented as an
// array.  Note that here we are using gsl - GNU Scientific Library -
// for diagonalisation.

void diagonalize(double* data, const int dim, double* res, double* eig)
{
   gsl_matrix_view m = gsl_matrix_view_array (data, dim, dim);
   gsl_vector *eval = gsl_vector_alloc (dim);
   gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
   gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dim);
   gsl_eigen_symmv (&m.matrix, eval, evec, w);
   gsl_eigen_symmv_free (w);
   gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
   
   for (int i = 0; i != dim; ++i) {
      double eval_i = gsl_vector_get (eval, i);
      gsl_vector_view evec_i = gsl_matrix_column (evec, i);

      // res stores eigenvectors as columns
      for (int j = 0; j < dim; ++j) {
         res[i + dim * j] = gsl_vector_get(&evec_i.vector, j);
      }
      // eig stores eigenvalues
      eig[i] = eval_i;
   }

   gsl_vector_free (eval);
   gsl_matrix_free (evec);
}

// Print eigenvalue information.

void printeig(const double* res, const double* eig, const int dim)
{
   for (int i = 0; i != dim; ++i) {
      cout << "eigenvalue = " << eig[i] << endl;
      cout << "eigenvector = " << endl;
      for (int j = 0; j != dim; ++j) {
         cout << res[i + dim * j] << endl;
      }
   }
}
