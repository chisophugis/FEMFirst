// I'm mostly looking at MFEM for guidance.
// Also Bathe's lecture: https://www.youtube.com/watch?v=oNqSzzycRhw
// (and the accompanying material on OCW)

#include <iomanip>
#include <iostream>

#include <algorithm>
#include <tuple>
#include <vector>

#include <assert.h>
#include <stdint.h>
#include <string.h>

// Following MFEM, we'll try to keep separate the two concepts:
// - general array datastructure for programmatic use
// - vector in the mathematical sense of the FEM formulation
using Vector = std::vector<double>;
template <typename T>
using Array = std::vector<T>;

void print(const Vector &V, std::ostream &OS) {
  for (double D : V)
    OS << D << " ";
  OS << "\n";
}

//=====================
// Roller spring system
//=====================

struct Element {
  int DoF1, DoF2;
  double K; // Spring constant.
};

struct RollerSpringSystem {
  int NumDoFs;
  Array<Element> Elements;
};

void readSystem(RollerSpringSystem &RSS, std::istream &IS) {
  std::string Identifier; // We read into this as a dummy variable.
  IS >> Identifier;
  if (Identifier != "dofs") {
    std::cerr << "Expecting `dofs`";
    return;
  }

  IS >> RSS.NumDoFs;

  IS >> Identifier;
  if (Identifier != "elements") {
    std::cerr << "Expecting `elements`";
    return;
  }

  int NumElements;
  IS >> NumElements;
  RSS.Elements.resize(NumElements);
  for (int i = 0, e = RSS.Elements.size(); i != e; ++i) {
    Element *E = &RSS.Elements[i];
    IS >> E->DoF1;
    IS >> E->DoF2;
    IS >> E->K;
    // TODO: write some useful macros:
    // ERROR(expr, "message: " << foo << "info");
    // ASSERT(expr, "message: " << foo << "info");
    assert((E->DoF1 < RSS.NumDoFs && E->DoF2 < RSS.NumDoFs) &&
           "Invalid DoF number");
  }
}

void printSystem(const RollerSpringSystem &RSS, std::ostream &OS) {
  OS << "dofs\n";
  OS << RSS.NumDoFs << "\n";

  OS << "elements\n";
  OS << RSS.Elements.size() << "\n";
  //OS << std::setprecision(18);
  for (int i = 0, e = RSS.Elements.size(); i != e; ++i) {
    const Element *E = &RSS.Elements[i];
    OS << E->DoF1 << " " << E->DoF2 << " " << E->K << "\n";
  }
}

//==============
// Sparse matrix
//==============

class SparseMatrixBuilder;
// Simple CSR (compressed sparse row) sparse matrix.
// Square.
class SparseMatrix {
public:
  int NumRows = -1;
  int NumCols = -1;

  // StoredElements and StoredElementsPositionInRow go hand-in-hand.
  // StoredElements[k] is the entry in column
  // StoredElementsPositionInRow[k] of row i where
  // k lies in the index range [RowSeparators[i],RowSeparators[i+1])
  // The number of elements stored is RowSeparators[Size].
  // RowSeparators[0] is always 0.
  // Since RowSeparators are "separators", there is actually one more
  // element than the number of rows.
  double *StoredElements = nullptr;
  int *StoredElementsPositionInRow = nullptr;
  int *RowSeparators = nullptr;

  // Think of StoredElements as a bunch of books on a bookshelf. We
  // then insert separators between each group of books that are part
  // of the same row (RowSeparators). RowSeparators[0] is on the left
  // of the first book (at position 0; hence RowSeparators[0] == 0
  // always).
  // RowSeparators[Size] is to the right of the last book, notionally
  // to the left of the "one past the end" book (the book at position
  // "number of books total"; hence RowSeparators[Size] is the number
  // of books).

  SparseMatrix() {}
  SparseMatrix(SparseMatrix &&RHS) {
    NumRows = RHS.NumRows;
    NumCols = RHS.NumCols;
    StoredElements = RHS.StoredElements;
    StoredElementsPositionInRow = RHS.StoredElementsPositionInRow;
    RowSeparators = RHS.RowSeparators;
  }
  ~SparseMatrix() {
    delete[] StoredElements;
    delete[] StoredElementsPositionInRow;
    delete[] RowSeparators;
  }

  double elem(int i, int j) const {
    if (double *Ep = elemPtr(i, j))
      return *Ep;
    return 0;
  }
  // I really don't like this function as a public interface; it's a
  // weird compromise between exposing internal implementation details
  // (i.e. which elements are present) and avoiding directly exposing
  // the underlying CSR.
  // Usually routines that need this could be significantly sped up by
  // utilizing the internal implementation details directly, and
  // should just be members.
  double *elemPtr(int i, int j) const {
    int Begin = RowSeparators[i];
    int End = RowSeparators[i + 1];
    int *B = StoredElementsPositionInRow + Begin;
    int *E = StoredElementsPositionInRow + End;
    int *El = std::lower_bound(B, E, j);
    if (El == E || *El != j)
      return nullptr;
    return &StoredElements[El - StoredElementsPositionInRow];
  }

  void print(std::ostream &OS) const {
    //OS << std::setprecision(4); // Make output more compact.
    for (int i = 0; i != NumRows; ++i) {
      for (int j = 0; j != NumCols; ++j)
        std::cout << elem(i, j) << "\t";
      std::cout << "\n";
    }
  }

  void mult(const Vector &X, Vector &Y) const {
    assert(X.size() == NumCols && "Invalid size for X");
    Y.resize(NumRows);
    int *RS = RowSeparators;
    for (int i = 0, e = NumRows; i != e; ++i) {
      double Sum = 0;
      for (int jj = RS[i], jje = RS[i + 1]; jj != jje; ++jj) {
        Sum += StoredElements[jj] * X[StoredElementsPositionInRow[jj]];
      }
      Y[i] = Sum;
    }
  }

private:
  friend class SparseMatrixBuilder;
};

// Could be made more generic.
// Basically like std::unique except that you can apply a custom
// function to "compress" elements that compare equal.
template <typename T, typename F1, typename F2>
T *compress(T *Begin, T *End, F1 Equal, F2 Compress) {
  if (Begin == End)
    return Begin;

  // Loop invariant:
  // [Begin,Cur+1) is the correct result if the input were just [Begin,Next)
  T *Cur = Begin;
  T *Next = Cur + 1;
  for (; Next != End; ++Next) {
    if (Equal(*Cur, *Next)) {
      Compress(*Cur, *Next);
      continue;
    }
    ++Cur;
    *Cur = *Next;
  }
  return Cur + 1;
}

// Notionally represents a bunch of individual entries summed
// together.
class SparseMatrixBuilder {
public:
  SparseMatrixBuilder(int R, int C) : NumRows(R), NumCols(C) {}
  void addToElementAt(int Row, int Col, double Val) {
    Entries.push_back({Row, Col, Val});
  }
  SparseMatrix finalize() {
    assertMatrixConformsToAssumptions();

    SparseMatrix Ret;
    Ret.NumRows = NumRows;
    Ret.NumCols = NumCols;
    std::sort(Entries.begin(), Entries.end(), [](Entry &LHS, Entry &RHS) {
      return std::tie(LHS.Row, LHS.Col) < std::tie(RHS.Row, RHS.Col);
    });

    Entry *NewEnd = compress(
        &*Entries.begin(), &*Entries.end(),
        [](Entry &LHS, Entry &RHS) {
          return std::tie(LHS.Row, LHS.Col) == std::tie(RHS.Row, RHS.Col);
        },
        [](Entry &LHS, Entry &RHS) { LHS.Val += RHS.Val; });
    Entries.resize(NewEnd - &*Entries.begin());

    Ret.StoredElements = new double[Entries.size()];
    Ret.StoredElementsPositionInRow = new int[Entries.size()];
    Ret.RowSeparators = new int[NumRows + 1];
    int CurSep = 0;
    Ret.RowSeparators[0] = 0;
    for (int i = 0, e = Entries.size(); i != e; ++i) {
      Ret.StoredElements[i] = Entries[i].Val;
      Ret.StoredElementsPositionInRow[i] = Entries[i].Col;
      if (Entries[i].Row != CurSep) {
        CurSep = Entries[i].Row;
        // We assume at least one entry in each row, so this visits
        // sequential indices.
        Ret.RowSeparators[CurSep] = i;
      }
    }
    assert(CurSep == (NumRows - 1));
    Ret.RowSeparators[NumRows] = Entries.size();
    return Ret;
  };

private:
  int NumRows;
  int NumCols;
  struct Entry {
    int Row, Col;
    double Val;
  };
  Array<Entry> Entries;

  void assertMatrixConformsToAssumptions() const {
    // We assume that the matrix has an entry in at least every row.
    Array<int8_t> Buf(NumRows, 0);
    for (const Entry &E : Entries)
      Buf[E.Row] = 1;
    assert(std::all_of(Buf.begin(), Buf.end(), [](bool B) { return B; }) &&
           "Matrix has rows that vanish. Disconnected graph?");
  }
};

//=========
// Assembly
//=========

SparseMatrix assemble(const RollerSpringSystem &RSS) {
  SparseMatrixBuilder SMB(RSS.NumDoFs, RSS.NumDoFs);
  for (int i = 0, e = RSS.Elements.size(); i != e; ++i) {
    const Element *E = &RSS.Elements[i];
    SMB.addToElementAt(E->DoF1, E->DoF1, E->K);
    SMB.addToElementAt(E->DoF1, E->DoF2, -E->K);
    SMB.addToElementAt(E->DoF2, E->DoF1, -E->K);
    SMB.addToElementAt(E->DoF2, E->DoF2, E->K);
  }
  return SMB.finalize();
}

//====================
// Boundary conditions
//====================

void eliminateBoundaryConditions(SparseMatrix &Stiffness, Vector &Load,
                                 double Sol, int EssentialDoF) {
  // This could exploit the sparsity, but for now, just be naive.

  const int Ess = EssentialDoF;

  // For the purposes of solving the system, we don't care about the
  // row with the essential DoF, except that it is consistent.
  // Just zero it out (leaving the diagonal equal to 1), and hardcode
  // the load vector to be consistent with it.
  // With the matrix representing the physical system, this actually
  // changes the physical meaning.
  // I.e. in the context of the roller spring system, once we do this
  // modification, multiplying the solution by the stiffness matrix
  // will *not* correctly give the force on the cart with the
  // essential boundary condition.
  // MFEM has an overload that fills in a SparseMatrix from which the
  // removed entries can be recovered (notionally it just needs to be
  // added back in, but in order to avoid summing sparse matrices you
  // would likely just use distributivity and do two separate matvecs
  // and then add the two resulting vectors).

  for (int j = 0, je = Stiffness.NumCols; j != je; ++j)
    if (double *D = Stiffness.elemPtr(Ess, j))
      if (j != Ess)
        *D = 0;
  // Force Load[Ess] to be consistent.
  Load[Ess] = Sol;
  assert(Stiffness.elem(Ess, Ess) != 0 && "DoF does not play into solution?");
  *Stiffness.elemPtr(Ess, Ess) = 1;


  // Note that the step above ruins the symmetry of the stiffness
  // matrix. We can fix that by zeroing the corresponding column
  // (except for the diagonal element, which we already set above).
  // The way we do that is by exploiting that we already know the
  // solution value Sol for EssentialDoF, so we can incorporate the
  // forces due to EssentialDoF into the load vector directly, and
  // remove the coupling from the stiffness matrix.
  //
  // If you are wondering how we can add in a force contribution from
  // EssentialDoF without knowing the offset of the other DoF's, the
  // answer is that you are probably thinking about the force as:
  //     F1 = k (u1 - u2)
  // but from the perspective of the linear system, it is
  //     F1 = k u1 + (-k) u2
  // I.e. the contribution from each DoF can be computed independently
  // if we know the solution value for it already.
  // If, say, u2 is the essential DoF, we can then rewrite this as:
  //     F1' = (F1 - (-k) u2) = k u1
  for (int i = 0, e = Stiffness.NumRows; i != e; ++i) {
    if (i != Ess) {
      if (double *D = Stiffness.elemPtr(i, Ess)) {
        Load[i] -= Sol * (*D);
        *D = 0;
      }
    }
  }

  // The big idea here is that the forces are "supposed to be" the
  // givens in the roller-spring system; except for givens in the form
  // of essential boundary conditions. We are dealing with them now in
  // order to reduce the problem to a problem where only forces are
  // the givens.
  //
  //>So to summarize what we did in this function, we first took the
  // system
  // F1 = k (u1 - u2)
  // F2 = k (u2 - u1)
  // (a single spring with spring constant k)
  // which in matrix form is:
  // [  k  -k ] [ u1 ] = [ F1 ]
  // [ -k   k ] [ u2 ] = [ F2 ]
  //>and, say, knowing u2 = 5, we first say "let's force the load
  // vector to be consistent with u2 = 5; we know u2, so ignore the
  // effect on u2 from u1", getting us:
  // [  k  -k ] [ u1 ] = [ F1 ]
  // [  0   1 ] [ u2 ] = [ (the numerical value of u2, as a force) ]
  // We could also have done:
  // [  k  -k ] [ u1 ] = [   F1 ]
  // [  0   k ] [ u2 ] = [ k u2 ]
  // and MFEM actually has an option for EliminateRow and related
  // functions which allows choosing between them. For no particular
  // reason, I have done the former.
  //>Finally, we incorporate the column into the RHS, which makes the
  // matrix symmetric again (i.e. for the roller cart system, this
  // means that it corresponds to a proper physical system).
  // [  k   0 ] [ u1 ] = [ F1 - (-k) u2 ]
  // [  0   1 ] [ u2 ] = [ (the numerical value of u2, as a force) ]
  //
  // Note that this physical situation is completely analogous to many
  // others. Basically all the ones that are the lumped version of
  // Laplace's equation, i.e. a diffusion problem.
  // E.g the carts can be nodes in a resistive network, with the
  // forces being currents injected at each node; the solution is just
  // the steady-state of the diffusion given the boundary conditions.
  // Or the springs could be capacitors, and the forces be an amount
  // of charge injected at each node.
  // I.e. the force represents a flux across one of the boundaries.
}

//=======================
// Linear system solution
//=======================

// Z = X - Y
void subtract(const Vector &X, const Vector &Y, Vector &Z) {
  // Note: one of X or Y can be the destination.
  for (int i = 0, e = Z.size(); i != e; ++i)
    Z[i] = X[i] - Y[i];
}

double dot(const Vector &X, const Vector &Y) {
  double Sum = 0;
  for (int i = 0, e = X.size(); i != e; ++i)
    Sum += X[i] * Y[i];
  return Sum;
}

// Z = X + a Y
void add(const Vector &X, double a, const Vector &Y, Vector &Z) {
  for (int i = 0, e = X.size(); i != e; ++i)
    Z[i] = X[i] + a * Y[i];
}

void cg(const SparseMatrix &A, const Vector &B, Vector &X, int max_iters,
        double rel_tol, double abs_tol) {
  assert(A.NumRows == B.size() && "Invalid size for B");
  X.assign(B.size(), 0);
  Vector D(A.NumRows);
  Vector R(A.NumRows);
  Vector AD(A.NumRows);
  // d = r = b - A x
  A.mult(X, R);
  subtract(B, R, R);
  D = R;
  double BetaNumerator = dot(R, R);
  double InitialRNormSquared = BetaNumerator;
  double TargetRNormSquared =
      std::max(rel_tol * rel_tol * BetaNumerator, abs_tol * abs_tol);
  std::cout << "TargetRNormSquared = " << TargetRNormSquared << "\n";
  std::cout << "(Rel = " << (rel_tol * rel_tol * BetaNumerator)
            << ", Abs = " << (abs_tol * abs_tol) << ")\n";
  for (int i = 0;; ++i) {
    std::cout << "Iter: " << std::setw(3) << i << "   (r,r) = " << BetaNumerator
              << " (compared to " << InitialRNormSquared << ")\n";
    if (BetaNumerator <= TargetRNormSquared) {
      std::cout << "Converged.\n";
      return;
    }
    if (i == max_iters) {
      std::cout << "No Convergence.\n";
      return;
    }
    double AlphaNumerator = BetaNumerator;
    A.mult(D, AD);
    double AlphaDenominator = dot(D, AD);
    double Alpha = AlphaNumerator / AlphaDenominator;
    add(X, Alpha, D, X);
    add(R, -Alpha, AD, R);
    BetaNumerator = dot(R, R);
    double Beta = BetaNumerator / AlphaNumerator;
    add(R, Beta, D, D);
  }
}

int main(int, char **) {
  std::cout << std::setprecision(18);
  std::cout.setf(std::ios_base::scientific);
  RollerSpringSystem RSS;
  readSystem(RSS, std::cin);
  //printSystem(RSS, std::cout);
  //const int kDim = 5;
  //SparseMatrixBuilder SMB(kDim, kDim);
  //for (int i = 0; i != kDim; ++i) {
  //  SMB.addToElementAt(i, i, 1.25);
  //  if (i + 1 < kDim)
  //    SMB.addToElementAt(i, i + 1, 2.25);
  //}
  //SparseMatrix SM = SMB.finalize();
  //SM.print(std::cout);

  SparseMatrix Stiffness = assemble(RSS);
  //std::cout << "Stiffness (initially assembled):\n";
  //Stiffness.print(std::cout);

  // Essential/Dirichlet boundary conditions...
  Vector Load(Stiffness.NumRows, 0);
  // TODO: Initialize in a more interesting way.
  for (int i = 0, e = Load.size(); i != e; ++i)
    Load[i] = i;
  Vector X(Stiffness.NumCols, 0);
  // FIXME: hardcoded.
  // Should be part of the problem description.
  Vector EssentialDoFs = {0};
  // Reference cart is at "0" position. Actual displacement doesn't
  // matter though for only one essential boundary condition though.
  X[EssentialDoFs[0]] = 1;

  for (int Ess : EssentialDoFs)
    eliminateBoundaryConditions(Stiffness, Load, X[Ess], Ess);
  //std::cout << "Stiffness (ess. bdr. cond. eliminated):\n";
  //Stiffness.print(std::cout);
  //std::cout << "Load: ";
  //print(Load, std::cout);
  //std::cout << "X: ";
  //print(X, std::cout);

  std::cerr << "Starting cg.\n";
  cg(Stiffness, Load, X, Stiffness.NumRows, 1e-14, 0);

  std::cout << "X: ";
  print(X, std::cout);

  // TODO: need better way to control what is printed.
  // Maybe just command line options, one for each thing?
  // There is a tradeoff between human-understandability and
  // machine-readability. Much of what I want to do with the output is
  // just feed it into Mathematica, so options should be e.g. "print
  // the initially assembled stiffness matrix" and it prints *just
  // that*.

  return 0;
}
