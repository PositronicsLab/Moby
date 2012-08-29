#include <ctype.h>
#include <cmath>
#include <assert.h> 
#include <stdio.h>
#include <Moby/Constants.h>
#include <Moby/PathLCPSolver.h>

extern "C"
{
#include "Options.h"
#include "Error.h"
#include "Macros.h"
#include "Memory.h"
#include "Output.h"
#include "Output_Interface.h"
#include "Path.h"
#include "PathOptions.h" 
#include "Types.h" 
#include "MCP_Interface.h"
}

using namespace Moby;
using std::vector;
using boost::shared_array;
using boost::shared_ptr;

#if SAFESTATIC != static
#error Path LCP solver likely not reentrant!
#endif

// for fixing ctype link errors
__const unsigned short int *__ctype_b = *(__ctype_b_loc());
__const __int32_t *__ctype_tolower = *(__ctype_toupper_loc());
__const __int32_t *__ctype_toupper = *(__ctype_tolower_loc());

typedef struct
{
  int variables;

  int n;
  int nnz;

  double *z;
  double *lb;
  double *ub;

  int *m_start;
  int *m_len;
  int *m_row;
  double *m_data;

  double *q;
} Problem;

static Problem problem;
static int filled;

static CB_FUNC(void) start(void *v)
{
  filled = 0;
  return;
}

static CB_FUNC(void) problem_size(void *v, int *n, int *nnz)
{
  *n = problem.n;
  *nnz = problem.nnz + 1;
  return;
}

static CB_FUNC(void) bounds(void *v, int n, double *z, double *lb, double *ub)
{
  int i;

  for (i = 0; i < n; i++) {
    z[i] = problem.z[i];
    lb[i] = problem.lb[i];
    ub[i] = problem.ub[i];
  }
  return;
}

static CB_FUNC(int) function_evaluation(void *v, int n, double *z, double *f)
{
  int col, colStart, colEnd, row;
  double value;

  for (col = 0; col < n; col++) {
    f[col] = problem.q[col];
  }

  for (col = 0; col < n; col++) {
    value = z[col];

    if (value != 0) {
      colStart = problem.m_start[col] - 1;
      colEnd = colStart + problem.m_len[col];

      while(colStart < colEnd) {
	row = problem.m_row[colStart] - 1;
	f[row] += problem.m_data[colStart]*value;
	colStart++;
      }
    }
  }

  return 0;
}

static CB_FUNC(int) jacobian_evaluation(void *v, int n, double *z, int wantf, 
			                double *f, int *nnz,
                                        int *col_start, int *col_len, 
			                int *row, double *data)
{
  int element;

  if (wantf) {
    function_evaluation(v, n, z, f);
  }

  if (!filled) {
    for (element = 0; element < problem.n; element++) {
      col_start[element] = problem.m_start[element];
      col_len[element] = problem.m_len[element];
    }

    for (element = 0; element < problem.nnz; element++) {
      row[element] = problem.m_row[element];
      data[element] = problem.m_data[element];
    }

    filled = 1;
  }

  *nnz = problem.nnz;
  return 0;
}

static CB_FUNC(void) mcp_typ(void *d, int nnz, int *typ)
{
  int i;

  for (i = 0; i < nnz; i++) {
    typ[i] = PRESOLVE_LINEAR;
  }
  return;
}

static MCP_Interface mcp_interface =
{
  NULL,
  problem_size, bounds,
  function_evaluation, jacobian_evaluation,
  start, NULL,
  NULL, NULL,
  NULL
};

static Presolve_Interface mcp_presolve =
{
  NULL,
  NULL, NULL,
  NULL, NULL,
  mcp_typ,
  NULL
};

static void install_interface(MCP *m)
{
  MCP_SetInterface(m, &mcp_interface);
  MCP_SetPresolveInterface(m, &mcp_presolve);
  return;
}

static void sort(int rows, int cols, int elements,
		 int *row, int *col, double *data)
{
  double *m_data;
  int *m_start;
  int *m_len;
  int *m_row;

  int i, cs, ce;

  m_start = (int *)Memory_Allocate(sizeof(int)*(cols + 1));
  m_len = (int *)Memory_Allocate(sizeof(int)*(cols + 1));
  m_row = (int *)Memory_Allocate(sizeof(int)*(elements + 1));
  m_data = (double *)Memory_Allocate(sizeof(double)*(elements + 1));

  for(i = 0; i < cols; i++) {
    m_len[i] = 0;
  }

  for(i = 0; i < elements; i++) {
    if ((col[i] < 1) || (col[i] > cols)) {
      Error("column incorrect.\n");
    }

    if ((row[i] < 1) || (row[i] > rows)) {
      Error("column incorrect.\n");
    }

    m_len[col[i] - 1]++;
  }

  m_start[0] = 0;
  for (i = 1; i < cols; i++) {
    m_start[i] = m_start[i - 1] + m_len[i - 1];
    m_len[i - 1] = 0;
  }
  m_len[i - 1] = 0;

  for(i = 0; i < elements; i++) {
    cs = col[i] - 1;
    ce = m_start[cs] + m_len[cs];
    m_row[ce] = row[i];
    m_data[ce] = data[i];
    m_len[cs]++;
  }

  elements = 0;
  for (i = 0; i < cols; i++) {
    cs = m_start[i];
    ce = cs + m_len[i];

    while (cs < ce) {
      row[elements] = m_row[cs];
      col[elements] = i + 1;
      data[elements] = m_data[cs];
      elements++;
      cs++;
    }
  }

  Memory_Free(m_data);
  Memory_Free(m_row);
  Memory_Free(m_len);
  Memory_Free(m_start);
  return;
}

static void create(int variables,
	           int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
	           double *z, double *lb, double *ub)
{
  double inf;
  int m_index;
  int m_count;
  int i;

  inf = 1e20;

  problem.n = variables;
  problem.nnz = m_nnz;

  problem.z = (double *)Memory_Allocate(sizeof(double)*problem.n);
  problem.lb = (double *)Memory_Allocate(sizeof(double)*problem.n);
  problem.ub = (double *)Memory_Allocate(sizeof(double)*problem.n);

  problem.m_start = (int *)Memory_Allocate(sizeof(int)*problem.n);
  problem.m_len = (int *)Memory_Allocate(sizeof(int)*problem.n);
  problem.m_row = (int *)Memory_Allocate(sizeof(int)*problem.nnz+1);
  problem.m_data = (double *)Memory_Allocate(sizeof(double)*problem.nnz+1);

  problem.q = (double *)Memory_Allocate(sizeof(double)*problem.n);

  sort(variables, variables, m_nnz, m_i, m_j, m_ij);

  for (i = 0; i < variables; i++) {
    problem.z[i] = z[i];

    problem.q[i] = q[i];
    problem.lb[i] = lb[i];
    problem.ub[i] = ub[i];
  }

  m_index = 0;
  m_count = 0;
  for (i = 0; i < variables; i++) {
    problem.m_start[i] = m_count + 1;
    problem.m_len[i] = 0;

    while ((m_index < m_nnz) && (m_j[m_index] <= i + 1)) {
      if (m_ij[m_index] != 0) {
	problem.m_len[i]++;
	problem.m_row[m_count] = m_i[m_index];
	problem.m_data[m_count] = m_ij[m_index];
	m_count++;
      }
      m_index++;
    }
  }
  problem.nnz = m_count;
  return;
}

static void destroy(void)
{
  Memory_Free(problem.z);
  Memory_Free(problem.lb);
  Memory_Free(problem.ub);
  Memory_Free(problem.m_start);
  Memory_Free(problem.m_len);
  Memory_Free(problem.m_row);
  Memory_Free(problem.m_data);
  Memory_Free(problem.q);
  return;
}

void SimpleLCP(int variables,
               int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
	       double *lb, double *ub,
               MCP_Termination *status, double *z, double tol)
{
  Options_Interface *o;
  MCP *m;
  Information info;

  double *x;
  double dnnz;
  int i;

  o = Options_Create();
  Path_AddOptions(o);
  Options_Default(o);
  Options_SetDouble(o, "convergence_tolerance", tol);
  Options_SetBoolean(o, "output", False);

//  Output_Printf(Output_Log | Output_Status | Output_Listing,
//		"%s: LCP Link\n", Path_Version());

  create(variables, m_nnz, m_i, m_j, m_ij, q, z, lb, ub);
     
  if (problem.n == 0) {
//    Output_Printf(Output_Log | Output_Status,
//		  "\n ** EXIT - solution found (degenerate model).\n");
    (*status) = MCP_Solved;
    Options_Destroy(o);
    return;
  }

  dnnz = MIN(1.0*problem.nnz, 1.0*problem.n*problem.n);
  if (dnnz > std::numeric_limits<int>::max()) {
//    Output_Printf(Output_Log | Output_Status,
//		  "\n ** EXIT - model too large.\n");
    (*status) = MCP_Error;
    Options_Destroy(o);
    return;
  }
  problem.nnz = (int) dnnz;

//  Output_Printf(Output_Log | Output_Status | Output_Listing,
//		"%d row/cols, %d non-zeros, %3.2f%% dense.\n\n",
//		problem.n, problem.nnz, 
//		100.0*problem.nnz/(1.0*problem.n*problem.n));

  m = MCP_Create(problem.n, problem.nnz+1);
  MCP_Jacobian_Structure_Constant(m, True);
  install_interface(m);

//  Options_Read(o, "path.opt");
//  Options_Display(o);

//  info.generate_output = Output_Log | Output_Status | Output_Listing;
  info.generate_output = 0;
  info.use_start = True;
  info.use_basics = True;

  (*status) = Path_Solve(m, &info);

  x = MCP_GetX(m);
     
  for (i = 0; i < variables; i++) {
    z[i] = x[i];
  }

  MCP_Destroy(m);
  destroy();

  Options_Destroy(o);
  return;
}

bool PathLCPSolver::solve_lcp(const MatrixNN& MM, const VectorN& q, VectorN& z, double tol)
{
  const Real INF = std::numeric_limits<Real>::max();
  int i, j, k;
  int variables;
  shared_array<int> m_i, m_j;
  shared_array<double> dq, lb, ub, m_ij, x_end;
  MCP_Termination termination;
   
  Output_Printf(Output_Log | Output_Status | Output_Listing, "%s: Standalone-C Link\n", Path_Version());

  // We need a sparse representation of the matrix
  assert(q.size() == MM.size());
  assert(q.size() == z.size());

  // Copy q to dq
  variables = q.size();
  dq = shared_array<double>(new double[variables+1]);
  for (i = 0;i < variables; i++)
    dq[i] = q[i];

  // Copy q to dq
  variables = q.size();
  dq = shared_array<double>(new double[variables+1]);
  for (i = 0;i < variables; i++)
    dq[i] = q[i];

  int num_non_zeroes = 0;
  for (i = 0;i < variables; i++)
    for (j = 0;j < variables; j++)
      if (std::fabs(MM(i,j)) > NEAR_ZERO)
        num_non_zeroes++;

  m_i = shared_array<int>(new int[num_non_zeroes + 1]);
  m_j = shared_array<int>(new int[num_non_zeroes + 1]);
  m_ij = shared_array<double>(new double[num_non_zeroes + 1]);

  k = 0;
  for (i = 0;i < variables; i++)
  {
    for (j = 0;j < variables; j++)
    {
      if (std::fabs(MM(i,j)) > NEAR_ZERO)
      {
        m_i[k] = i+1;
        m_j[k] = j+1;
        m_ij[k] = MM(i,j);
        k++;
      }
    }
  }
    
  lb = shared_array<double>(new double[variables+1]);
  ub = shared_array<double>(new double[variables+1]);
  for (i = 0;i < variables; i++)
  {
    lb[i] = 0.0;
    ub[i] = INF;
  }

  x_end = shared_array<double>(new double[variables+1]);

  SimpleLCP(variables, num_non_zeroes, m_i.get(), m_j.get(), m_ij.get(), dq.get(), lb.get(), ub.get(), &termination, x_end.get(), tol);

  // We now copy x_end into z
  z.resize(variables); 
  for (i = 0;i < variables; i++)
    z[i] = x_end[i];

  return (termination == MCP_Solved);
}



