#include "util.hpp"

/*
 * helper function to print a matrix
 */
void printMatrix(gsl_matrix* mat) {
  printMatrix(stdout, mat);
}
void printMatrix(FILE* stream, gsl_matrix* mat) {
  for(unsigned int i = 0; i < mat->size1; i++) {
    for(unsigned int j = 0; j < mat->size2; j++) {
      fprintf(stream, "%0.15g \t", gsl_matrix_get(mat, i, j));
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}

/*
 * helper function to print a row vector
 */
void printRowVector(gsl_vector* vec) {
  printRowVector(stdout, vec);
}
void printRowVector(FILE* stream, gsl_vector* vec) {
  for(unsigned int i = 0; i < vec->size; i++) {
    fprintf(stream, "%.40f \t", gsl_vector_get(vec, i));
  }
  fprintf(stream, "\n");
}

/*
 * helper function to print a column vector
 */
void printColVector(gsl_vector* vec) {
  printColVector(stdout, vec);
}
void printColVector(FILE* stream, gsl_vector* vec) {
  for(unsigned int i = 0; i < vec->size; i++) {
    fprintf(stream, "%.40f\n", gsl_vector_get(vec, i));
  }
  fprintf(stream, "\n");
}

/*
 * helper function to take in a src matrix and col number,
 * and turn the specified column vector into
 * a matrix, with the src elements on the diagonals and
 * zeros in the off-diagonals, stored in dest.
 * Assumes src and dest are the same dimensions
 */
void createDiagFromCol(gsl_matrix* dest, int colNum, gsl_matrix* src) {
  // extract column
  gsl_vector_view col = gsl_matrix_column(src, colNum);

  // reset destination matrix
  gsl_matrix_set_zero(dest);

  // extract diagonal elements and set them to be the src col elements
  gsl_vector_view diagElements = gsl_matrix_diagonal(dest);
  gsl_vector_memcpy(&diagElements.vector, &col.vector);
}

/*
 * helper function to raise src matrix to power, stores in dest.
 * assumes src and dest are the same dimension, and assumes
 * src (and dest) are multipliable
 */
void raiseMatrixToPower(gsl_matrix* dest, int power, gsl_matrix* src) {
  // set up dest as a copy of src
  gsl_matrix_memcpy(dest, src);
  gsl_matrix* currRes = gsl_matrix_alloc(src->size1, src->size2);

  for(int i = 0; i < power; i++) {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, src, dest, 0.0, currRes); // Double GEneral Matrix Multiplication
    gsl_matrix_memcpy(dest, currRes);
  }
  gsl_matrix_free(currRes);
}

/*
 * helper function to calculate the matrix exponential, P = e^(Qt),
 * via eigendecomposition. Q must be diagonalizable and square (this is assumed, no checks are done)
 *
 * P = e^(Qt)
 *   = X e^(Dt) X^-1, where Q = X D X^(-1), and X is the eigenvectors of Q, D is the diagonal eigenvalues, X^(-1) is the inverse of X
 *   and each entry of e^(Dt) is just e^(D[i,i]*t)
 *
 * Matrix exponential is stored in destP, srcQ is unchanged
 * eigenWorkspace, eval, evec, realEvecMat, diagMat, LUdecompMat, perm, inverseMat are all intermediates that can and should be reused/not realloc every time this is called. Contents are not guaranteed
 * returns GSL_NAN if an error is thrown in any matrix operation
 */
double matrixExponential(gsl_matrix* destP, gsl_matrix* srcQ, double t, gsl_eigen_nonsymmv_workspace* eigenWorkspace, gsl_vector_complex* eval, gsl_matrix_complex* evec, gsl_matrix* realEvecMat, gsl_matrix* diagMat, gsl_matrix* LUdecompMat, gsl_permutation* perm, gsl_matrix* inverseMat) {
  // turn off error handler to get rid of "gsl: lu.c:266: ERROR: matrix is singular" message that stops program execution
  // see https://www.gnu.org/software/gsl/doc/html/err.html#c.gsl_set_error_handler_off
  gsl_error_handler_t* errHandler = gsl_set_error_handler_off();

  // first find eigenvalues and eigenvectors: https://www.gnu.org/software/gsl/doc/html/eigen.html#c.gsl_eigen_nonsymmv
  // need an intermediate matrix to keep srcQ from changing in the eigen call, use LUdecompMat for convenience
  gsl_matrix_memcpy(LUdecompMat, srcQ);
  double status = gsl_eigen_nonsymmv(LUdecompMat, eval, evec, eigenWorkspace);
  if(status != GSL_SUCCESS) {
    std::cerr << "ERROR: could not calculate eigenvalues/eigenvectors in matrixExponential" << std::endl;

    // restore error handler
    gsl_set_error_handler(errHandler);
    return GSL_NAN;
  }

  // convert eval into diagMat
  gsl_matrix_set_zero(diagMat);
  gsl_vector_view evalReal = gsl_vector_complex_real(eval);
  for(unsigned int i = 0; i < eval->size; i++) {
    gsl_matrix_set(diagMat, i, i, gsl_vector_get(&evalReal.vector, i));
  }

  // extract real portions of evec
  for(unsigned int row = 0; row < evec->size1; row++) {
    for(unsigned int col = 0; col < evec->size2; col++) {
      gsl_matrix_set(realEvecMat, row, col, GSL_REAL(gsl_matrix_complex_get(evec, row, col)));
    }
  }

  // find inverse of eigenvectors matrix. from https://gist.github.com/bjd2385/7f4685e703f7437e513608f41c65bbd7
  int s;
  gsl_matrix_memcpy(LUdecompMat, realEvecMat);

  // Compute the LU decomposition of this matrix: https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_LU_decomp
  status = gsl_linalg_LU_decomp(LUdecompMat, perm, &s);
  if(status != GSL_SUCCESS) {
    std::cerr << "ERROR: could not calculate LU decomposition in matrixExponential" << std::endl;

    // restore error handler
    gsl_set_error_handler(errHandler);
    return GSL_NAN;
  }

  // Compute the  inverse of the LU decomposition: https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_LU_invert
  status = gsl_linalg_LU_invert(LUdecompMat, perm, inverseMat);
  if(status != GSL_SUCCESS) {
    std::cerr << "ERROR: could not invert LU decomposition in matrixExponential" << std::endl;

    // restore error handler
    gsl_set_error_handler(errHandler);
    return GSL_NAN;
  }

  // raise diagMat to the correct power
  double expValue = 0;
  for(unsigned int i = 0; i < diagMat->size1; i++) {
    expValue = gsl_matrix_get(diagMat, i, i);
    expValue = exp(expValue * t);
    gsl_matrix_set(diagMat, i, i, expValue);
  }
  
  // multiply everything back together: https://www.gnu.org/software/gsl/doc/html/blas.html#c.gsl_blas_dgemm
  // need an intermediate matrix for multiplication storage , use LUdecompMat for convenience
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, realEvecMat, diagMat, 0.0, LUdecompMat); // Double GEneral Matrix Multiplication
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, LUdecompMat, inverseMat, 0.0, destP); // Double GEneral Matrix Multiplication
  if(status != GSL_SUCCESS) {
    std::cerr << "ERROR: could not multiply matrices in matrixExponential" << std::endl;

    // restore error handler
    gsl_set_error_handler(errHandler);
    return GSL_NAN;
  }

  // restore error handler
  gsl_set_error_handler(errHandler);

  return status;
}

/*
 * Helper function to create state space, given number of cells and maxPloidy
 * from https://stackoverflow.com/a/28712605
 */
std::vector<std::string>* createStateSpace(int numCells, int maxPloidy) {
  // there are (maxPloidy + 1)^numCells states (+1 for 0 ploidy state)
  std::vector<std::string>* states = new std::vector<std::string>(std::pow(maxPloidy + 1, numCells));
  std::vector<int> ploidies;
  // for each ploidy, add it numCells times
  for(int i = 0; i <= maxPloidy; i++) {
    for(int j = 0; j < numCells; j++) {
      ploidies.push_back(i);
    }
  }
  std::sort(ploidies.begin(), ploidies.end());

  // create all permutations
  int stateCtr = 0;
  do
  {
    std::ostringstream os;
    for (int i = 0; i < numCells - 1; i++)
    {
      os << ploidies[i] << ",";
    }
    os << ploidies[numCells - 1];
    (*states)[stateCtr] = os.str();
    stateCtr++;
    std::reverse(ploidies.begin() + numCells, ploidies.end());
  } while (next_permutation(ploidies.begin(), ploidies.end()));
  return states;
}

/*
 * helper function to parse a state string (ex "0,0") into
 * a vector of int states stored in dest (ex [0, 0]).
 * based on https://stackoverflow.com/a/236976
 */
std::vector<int>* parseStateToInt(std::string state) {
  std::vector<int>* dest = new std::vector<int>();
  std::vector<std::string> parsedStrs;
  boost::split(parsedStrs, state, boost::is_any_of(","));
  for(std::vector<std::string>::iterator itr = parsedStrs.begin(); itr != parsedStrs.end(); ++itr) {
    dest->push_back(atoi((*itr).c_str()));
  }
  return dest;
}

/*
 * helper function to take the absolute value of a vector, in place
 */
void vectorAbsoluteValue(gsl_vector* vec) {
  for(unsigned int i = 0; i < vec->size; i++) {
    gsl_vector_set(vec, i, std::abs(gsl_vector_get(vec, i)));
  }
}

/*
 * helper function to compare depths, stored as doubles. Simplistic comparison is ok here
 */
bool compareDoubles(double depth1, double depth2, double epsilon) {
  return std::abs(depth1 - depth2) < epsilon;
}

/*
 * helper function to extract the sample name from a filename.
 * If the filename contains an SRA number, then returns the SRA number (if multiple SRA numbers exist,
 * only the first is returned).
 * Otherwise returns the basename without the file extension.
 * file parsing from: https://stackoverflow.com/a/21827261
 * regex from: http://www.cplusplus.com/reference/regex/regex_search/
 */
std::string parseSampleName(std::string filename) {
  boost::filesystem::path p(filename);
  std::string sampleName = p.stem().string();
  std::regex rgx("SRR[0-9]+");
  std::smatch match;
  while(std::regex_search(sampleName, match, rgx) ) {
    sampleName = match[0];
    break;
  }
  return sampleName;
}

/*
 * function to calculate element wise chi squared value of matrices. Chi-squared tests
 * should only be used if each expected value > ~5, but no such assumptions are applied here.
 * Assumes matrices are the same size
 * Returns sum_ij (expected_ij - observed_ij)^2 / expected_ij
 */
double calcChiSqOfMatrix(gsl_matrix* expected, gsl_matrix* observed) {
  double chiSq = 0;
  double expVal = 0;
  double obsVal = 0;
  for(unsigned int i = 0; i < expected->size1; i++) {
    for(unsigned int j = 0; j < expected->size2; j++) {
      expVal = gsl_matrix_get(expected, i, j);
      obsVal = gsl_matrix_get(observed, i, j);

      // if expect to see 0, skip this entry to avoid dividing by 0 ==> nan
      if(compareDoubles(expVal, 0)) {
        continue;
      }
      chiSq += std::pow(expVal - obsVal, 2) / expVal;
      //chiSq += std::pow(expVal - obsVal, 2) ;/// expVal; // TODO Wed 07 Oct 2020 08:32:47 PM PDT change back
    }
  }
  return chiSq;
}

/*
 * function to calculate element wise squared difference of matrices
 * Same as calcChiSqOfMatrix but doesn't divide by the expVal
 * Assumes matrices are the same size
 * Returns sum_ij (mat1_ij - mat2_ij)^2
 */
double calcSqDiffOfMatrix(gsl_matrix* mat1, gsl_matrix* mat2) {
  double sqDiff = 0;
  double val1 = 0;
  double val2 = 0;
  for(unsigned int i = 0; i < mat1->size1; i++) {
    for(unsigned int j = 0; j < mat1->size2; j++) {
      val1 = gsl_matrix_get(mat1, i, j);
      val2 = gsl_matrix_get(mat2, i, j);

      sqDiff += std::pow(val1 - val2, 2) ;
    }
  }
  return sqDiff;
}

/*
 * function to check if all entries in a matrix are the same.
 * for example, if the least squares transition parameter finding step
 * fails, the transition matrix can become uniform.
 * returns true if matrix is uniform, false otherwise
 */
bool checkUniformMatrix(gsl_matrix* mat, double epsilon) {
  double prevVal = gsl_matrix_get(mat, 0, 0);
  double currVal = 0;
  for(unsigned int row = 0; row < mat->size1; row++) {
    for(unsigned int col = 0; col < mat->size2; col++) {
      currVal = gsl_matrix_get(mat, row, col);
      // return false as soon a value is different from the 0,0 position
      if(std::abs(prevVal - currVal) > epsilon) {
        return false;
      }
    }
  }
  // if haven't previously returned, all entries must be the same as the 0,0 position, ergo every entry is the same and the matrix is uniform
  return true;
}

/*
 * function to check if val1 is a near double of val2 (within epsilon).
 * That is, checks if if val1 ~= 2 * val2. This is useful for
 * checking if average ploidies are approximately doubles of each other
 * (ie the lib scaling is wrong).
 *
 * returns true if val1 / val2 - 2 < epsilon, false otherwise
 */
bool checkNearDouble(double val1, double val2, double epsilon) {
  return (std::abs((val1 / val2) - 2.0) < epsilon);
}

