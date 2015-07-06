#include <iomanip>
#include <iostream>

#include <algorithm>
#include <tuple>
#include <vector>

#include <assert.h>
#include <string.h>


//=====================
// Roller spring system
//=====================

struct Element {
  int DoF1, DoF2;
  double K; // Spring constant.
};

struct RollerSpringSystem {
  int NumDoFs;
  std::vector<Element> Elements;
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
  }
}

void printSystem(const RollerSpringSystem &RSS, std::ostream &OS) {
  OS << "dofs\n";
  OS << RSS.NumDoFs << "\n";

  OS << "elements\n";
  OS << RSS.Elements.size() << "\n";
  OS << std::setprecision(18);
  for (int i = 0, e = RSS.Elements.size(); i != e; ++i) {
    const Element *E = &RSS.Elements[i];
    OS << E->DoF1 << " " << E->DoF2 << " " << E->K << "\n";
  }
}

//=============
// Sparse array
//=============

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
    StoredElements = RHS.StoredElements;
    StoredElementsPositionInRow = RHS.StoredElementsPositionInRow;
    RowSeparators = RHS.RowSeparators;
  }
  ~SparseMatrix() {
    delete[] StoredElements;
    delete[] StoredElementsPositionInRow;
    delete[] RowSeparators;
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
    SparseMatrix Ret;
    Ret.NumRows = NumRows;
    Ret.NumCols = NumCols;
    // We don't know the size of the other members of SparseMatrix
    // until later.
    Ret.RowSeparators = new int[NumRows];
    memset(Ret.RowSeparators, 0, NumRows * sizeof(int));
    if (Entries.size() == 0)
      return Ret;
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
    for (int i = 0, e = Entries.size(); i != e; ++i) {
      Ret.StoredElements[i] = Entries[i].Val;
      Ret.StoredElementsPositionInRow[i] = Entries[i].Col;
    }
    return Ret;
  };

private:
  int NumRows;
  int NumCols;
  struct Entry {
    int Row, Col;
    double Val;
  };
  std::vector<Entry> Entries;
};

//=========
// Assembly
//=========

SparseMatrix assemble(const RollerSpringSystem &RSS) {
  SparseMatrixBuilder SMB;
  for (int i = 0, e = RSS.Elements.size(); i != e; ++i) {
    const Element *E = &RSS.Elements[i];
    SMB.addToElementAt(E->DoF1, E->DoF1, -E->K);
    SMB.addToElementAt(E->DoF1, E->DoF2, E->K);
    SMB.addToElementAt(E->DoF2, E->DoF1, -E->K);
    SMB.addToElementAt(E->DoF2, E->DoF2, E->K);
  }
}

int main(int, char **) {
  RollerSpringSystem RSS;
  readSystem(RSS, std::cin);
  printSystem(RSS, std::cout);
  SparseMatrixBuilder SMB(100, 100);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SMB.addToElementAt(3, 3, 1.8);
  SparseMatrix SM = SMB.finalize();
  //SM = SMB.finalize();

  // Essential/Dirichlet boundary conditions...
  return 0;
}
